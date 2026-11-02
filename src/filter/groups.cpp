#include "arise/filter/groups.hpp"

#include "arise/core/error.hpp"
#include "arise/io/number.hpp"
#include "arise/io/window_csv.hpp"
#include "group_shifts_data.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <utility>

namespace arise::filter {

namespace {

constexpr std::array kAllTypes{
    nmr1d::CarbonType::C,
    nmr1d::CarbonType::CH,
    nmr1d::CarbonType::CH2,
    nmr1d::CarbonType::CH3,
};

// Columns of the measured table (data/group_shifts.tsv).
constexpr std::size_t kColumns = 12;
constexpr std::size_t kNColumn = 3;
constexpr std::size_t kMoleculesColumn = 4;
constexpr std::size_t kMinColumn = 5;
constexpr std::size_t kP1Column = 6;
constexpr std::size_t kP99Column = 10;
constexpr std::size_t kMaxColumn = 11;

std::vector<std::string_view> split(std::string_view text, char separator) {
    std::vector<std::string_view> parts;
    while (true) {
        const auto end = text.find(separator);
        parts.push_back(text.substr(0, end));
        if (end == std::string_view::npos) {
            return parts;
        }
        text = text.substr(end + 1);
    }
}

core::Error parse_error(std::size_t line, const std::string& message) {
    return {
        .code = core::ErrorCode::ParseError,
        .message = "hallmark windows, line " + std::to_string(line) + ": " + message,
    };
}

std::optional<chem::Group> group_named(std::string_view name) {
    for (const auto group : chem::kAllGroups) {
        if (chem::to_string(group) == name) {
            return group;
        }
    }
    return std::nullopt;
}

// The inverse of nmr1d::to_string(CarbonTypeSet) for non-empty sets: "{C, CH}".
std::optional<nmr1d::CarbonTypeSet> parse_types(std::string_view text) {
    if (text.size() < 2 || text.front() != '{' || text.back() != '}') {
        return std::nullopt;
    }
    nmr1d::CarbonTypeSet types;
    for (const auto name : split(text.substr(1, text.size() - 2), ',')) {
        const auto trimmed = name.starts_with(' ') ? name.substr(1) : name;
        bool known = false;
        for (const auto type : kAllTypes) {
            if (nmr1d::to_string(type) == trimmed) {
                types = types | nmr1d::CarbonTypeSet{type};
                known = true;
            }
        }
        if (!known) {
            return std::nullopt;
        }
    }
    return types;
}

std::optional<std::size_t> parse_count(std::string_view text) {
    if (text.empty() ||
        !std::all_of(text.begin(), text.end(), [](char c) { return c >= '0' && c <= '9'; })) {
        return std::nullopt;
    }
    std::size_t value = 0;
    for (const char c : text) {
        value = (value * 10) + static_cast<std::size_t>(c - '0');
    }
    return value;
}

core::Result<std::size_t> parse_positive_count(std::size_t line,
                                               std::span<const std::string_view> fields,
                                               std::size_t column, std::string_view name) {
    const auto count = parse_count(fields[column]);
    if (!count || *count == 0) {
        return parse_error(line, std::string{name} + " must be a positive whole number, got '" +
                                     std::string{fields[column]} + "'");
    }
    return *count;
}

// The shift columns (min to max), parsed; index i holds column kMinColumn + i.
core::Result<std::vector<double>> parse_shifts(std::size_t line,
                                               std::span<const std::string_view> fields) {
    std::vector<double> values;
    for (std::size_t column = kMinColumn; column < kColumns; ++column) {
        const auto value = io::parse_decimal(fields[column]);
        if (!value || !std::isfinite(value.value())) {
            return parse_error(line, "not a finite number: '" + std::string{fields[column]} + "'");
        }
        values.push_back(value.value());
    }
    const auto at = [&](std::size_t column) { return values[column - kMinColumn]; };
    if (at(kP1Column) > at(kP99Column)) {
        return parse_error(line, "1st percentile " + std::string{fields[kP1Column]} +
                                     " is above the 99th " + std::string{fields[kP99Column]});
    }
    if (at(kMinColumn) > at(kMaxColumn)) {
        return parse_error(line, "min " + std::string{fields[kMinColumn]} + " is above the max " +
                                     std::string{fields[kMaxColumn]});
    }
    return values;
}

core::Result<HallmarkWindow> parse_row(std::size_t line, std::string_view text) {
    const auto fields = split(text, '\t');
    if (fields.size() != kColumns) {
        return parse_error(line, "expected " + std::to_string(kColumns) +
                                     " tab-separated fields, got " + std::to_string(fields.size()));
    }
    const auto group = group_named(fields[0]);
    if (!group) {
        return parse_error(line, "unknown group '" + std::string{fields[0]} + "'");
    }
    const auto types = parse_types(fields[2]);
    if (!types || types->empty()) {
        return parse_error(line,
                           "DEPT types must be like {C, CH}, got '" + std::string{fields[2]} + "'");
    }
    const auto n = parse_positive_count(line, fields, kNColumn, "n");
    if (!n) {
        return n.error();
    }
    const auto molecules = parse_positive_count(line, fields, kMoleculesColumn, "molecules");
    if (!molecules) {
        return molecules.error();
    }
    const auto shifts = parse_shifts(line, fields);
    if (!shifts) {
        return shifts.error();
    }
    const bool small = molecules.value() < kMinMoleculesForPercentiles;
    const auto at = [&](std::size_t column) {
        return core::Ppm{shifts.value()[column - kMinColumn]};
    };
    return HallmarkWindow{
        .group = *group,
        .hallmark = std::string{fields[1]},
        .types = *types,
        .low = at(small ? kMinColumn : kP1Column),
        .high = at(small ? kMaxColumn : kP99Column),
        .n = n.value(),
        .molecules = molecules.value(),
        .source = small ? WindowSource::MeasuredRange : WindowSource::Percentiles,
    };
}

std::string lowercase(std::string_view text) {
    std::string lower{text};
    for (auto& c : lower) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return lower;
}

std::optional<chem::Group> group_named_ignoring_case(std::string_view name) {
    for (const auto group : chem::kAllGroups) {
        if (lowercase(chem::to_string(group)) == lowercase(name)) {
            return group;
        }
    }
    return std::nullopt;
}

std::string group_names() {
    std::string names;
    for (const auto group : chem::kAllGroups) {
        names.append(names.empty() ? "" : ", ").append(chem::to_string(group));
    }
    return names;
}

core::Error lab_error(std::size_t line, const std::string& message) {
    return {
        .code = core::ErrorCode::ParseError,
        .message = "lab windows, line " + std::to_string(line) + ": " + message,
    };
}

std::optional<nmr1d::CarbonType> type_named_ignoring_case(std::string_view name) {
    for (const auto type : kAllTypes) {
        if (lowercase(nmr1d::to_string(type)) == lowercase(name)) {
            return type;
        }
    }
    return std::nullopt;
}

// "C CH" (any case, spaces between); blank means every type.
core::Result<nmr1d::CarbonTypeSet> parse_dept_list(std::size_t line, std::string_view text) {
    nmr1d::CarbonTypeSet types;
    for (const auto token : split(text, ' ')) {
        if (token.empty()) {
            continue;
        }
        const auto type = type_named_ignoring_case(token);
        if (!type) {
            return lab_error(line, "unknown DEPT type '" + std::string{token} +
                                       "' (use C, CH, CH2 or CH3, separated by spaces)");
        }
        types = types | nmr1d::CarbonTypeSet{*type};
    }
    return types.empty() ? nmr1d::CarbonTypeSet::all() : types;
}

std::string default_label(chem::Group group) {
    std::string label;
    for (const auto name : chem::hallmarks(group)) {
        label.append(label.empty() ? "" : " or ").append(name);
    }
    return label.empty() ? "13C window" : label;
}

core::Result<HallmarkWindow> lab_window(const io::WindowRow& row) {
    const auto group = group_named_ignoring_case(row.group);
    if (!group) {
        return lab_error(row.line,
                         "unknown group '" + row.group + "'; the groups are: " + group_names());
    }
    auto types = parse_dept_list(row.line, row.dept);
    if (!types) {
        return std::move(types).error();
    }
    return HallmarkWindow{
        .group = *group,
        .hallmark = row.hallmark.empty() ? default_label(*group) : row.hallmark,
        .types = types.value(),
        .low = row.low,
        .high = row.high,
        .n = 0,
        .molecules = 0,
        .source = WindowSource::Lab,
    };
}

// One decimal, rounded half away from zero (no printf, so no locale).
std::string fixed1(double value) {
    const long long tenths = std::llround(value * 10.0);
    const long long magnitude = tenths < 0 ? -tenths : tenths;
    return (tenths < 0 ? "-" : "") + std::to_string(magnitude / 10) + "." +
           std::to_string(magnitude % 10);
}

// "C CH" for the lab table's dept column.
std::string dept_list(nmr1d::CarbonTypeSet types) {
    std::string list;
    for (const auto type : kAllTypes) {
        if (types.contains(type)) {
            list.append(list.empty() ? "" : " ").append(nmr1d::to_string(type));
        }
    }
    return list;
}

bool inside(const nmr1d::CarbonAssignment& peak, const HallmarkWindow& window) {
    return peak.shift >= window.low && peak.shift <= window.high &&
           !(peak.possible & window.types).empty();
}

}  // namespace

core::Result<std::vector<HallmarkWindow>> parse_hallmark_windows(std::string_view text) {
    std::vector<HallmarkWindow> windows;
    std::set<std::pair<chem::Group, std::string>> seen;
    std::size_t line = 0;
    for (auto content : split(text, '\n')) {
        ++line;
        if (content.ends_with('\r')) {
            content.remove_suffix(1);
        }
        if (content.empty() || content.front() == '#') {
            continue;
        }
        auto window = parse_row(line, content);
        if (!window) {
            return std::move(window).error();
        }
        if (!seen.insert({window.value().group, window.value().hallmark}).second) {
            return parse_error(line, std::string{chem::to_string(window.value().group)} + " / " +
                                         window.value().hallmark + " is listed twice");
        }
        windows.push_back(std::move(window).value());
    }
    if (windows.empty()) {
        return core::Error{
            .code = core::ErrorCode::ParseError,
            .message = "hallmark windows: the table has no rows",
        };
    }
    return windows;
}

core::Result<std::vector<HallmarkWindow>> measured_hallmark_windows() {
    return parse_hallmark_windows(detail::kMeasuredGroupShifts);
}

core::Result<std::vector<HallmarkWindow>> parse_lab_windows(std::string_view text) {
    auto rows = io::parse_window_csv(text);
    if (!rows) {
        return std::move(rows).error();
    }
    std::vector<HallmarkWindow> windows;
    std::set<std::pair<chem::Group, std::string>> seen;
    for (const auto& row : rows.value()) {
        auto window = lab_window(row);
        if (!window) {
            return std::move(window).error();
        }
        if (!seen.insert({window.value().group, window.value().hallmark}).second) {
            return lab_error(row.line, std::string{chem::to_string(window.value().group)} + " / " +
                                           window.value().hallmark + " is listed twice");
        }
        windows.push_back(std::move(window).value());
    }
    return windows;
}

std::vector<HallmarkWindow> merge_windows(std::span<const HallmarkWindow> built_in,
                                          std::span<const HallmarkWindow> lab) {
    std::vector<HallmarkWindow> merged;
    for (const auto group : chem::kAllGroups) {
        const bool lab_has_group = std::any_of(
            lab.begin(), lab.end(), [&](const HallmarkWindow& w) { return w.group == group; });
        for (const auto& window : lab_has_group ? lab : built_in) {
            if (window.group == group) {
                merged.push_back(window);
            }
        }
    }
    return merged;
}

std::string lab_table_template() {
    std::string text =
        "# ARISE 13C suggestion windows: a peak between low_ppm and high_ppm, from a carbon of\n"
        "# one of the dept types (blank = any), suggests the group. Edit the numbers, add rows,\n"
        "# or delete rows. A group you list replaces all of ARISE's windows for that group;\n"
        "# groups you leave out keep ARISE's. dept: C, CH, CH2, CH3 separated by spaces.\n"
        "# Built-in windows measured from nmrshiftdb2 (www.nmrshiftdb.org), made available\n"
        "# under the nmrshiftdb2 Database License.\n"
        "group,hallmark,dept,low_ppm,high_ppm\n";
    const auto windows = measured_hallmark_windows();
    for (const auto& w : windows.value_or(std::vector<HallmarkWindow>{})) {
        text.append(chem::to_string(w.group)).append(",").append(w.hallmark).append(",");
        text.append(dept_list(w.types)).append(",").append(fixed1(w.low.value())).append(",");
        text.append(fixed1(w.high.value())).append("\n");
    }
    return text;
}

std::vector<GroupSuggestion> suggest_groups(std::span<const nmr1d::CarbonAssignment> peaks,
                                            std::span<const HallmarkWindow> windows) {
    std::vector<GroupSuggestion> suggestions;
    for (const auto group : chem::kAllGroups) {
        GroupSuggestion suggestion{.group = group, .evidence = {}};
        for (std::size_t p = 0; p < peaks.size(); ++p) {
            for (const auto& window : windows) {
                if (window.group == group && inside(peaks[p], window)) {
                    suggestion.evidence.push_back(
                        {.peak = p, .shift = peaks[p].shift, .window = window});
                }
            }
        }
        if (!suggestion.evidence.empty()) {
            suggestions.push_back(std::move(suggestion));
        }
    }
    return suggestions;
}

std::string_view to_string(Tier tier) noexcept {
    switch (tier) {
        case Tier::Active:
            return "active";
        case Tier::BackBurner:
            return "back-burner";
        case Tier::Eliminated:
            return "eliminated";
    }
    return "unknown";
}

core::Result<Verdict> check_confirmed_groups(const chem::Molecule& molecule,
                                             const ConfirmedGroups& confirmed) {
    for (const auto group : confirmed.present) {
        if (std::find(confirmed.absent.begin(), confirmed.absent.end(), group) !=
            confirmed.absent.end()) {
            return core::Error{
                .code = core::ErrorCode::InvalidArgument,
                .message =
                    std::string{chem::to_string(group)} + " is confirmed both present and absent",
            };
        }
    }
    Verdict verdict;
    for (const auto group : confirmed.present) {
        if (chem::count_group(molecule, group) == 0) {
            verdict.reasons.push_back("confirmed present: " + std::string{chem::to_string(group)} +
                                      "; the candidate has none");
        }
    }
    for (const auto group : confirmed.absent) {
        if (const int count = chem::count_group(molecule, group); count > 0) {
            verdict.reasons.push_back("confirmed absent: " + std::string{chem::to_string(group)} +
                                      "; the candidate has " + std::to_string(count));
        }
    }
    if (!verdict.reasons.empty()) {
        verdict.tier = Tier::Eliminated;
    }
    return verdict;
}

}  // namespace arise::filter

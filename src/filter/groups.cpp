#include "arise/filter/groups.hpp"

#include "arise/core/error.hpp"
#include "arise/io/number.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <set>
#include <string>
#include <utility>

namespace arise::filter {

namespace {

constexpr std::size_t kColumns = 11;
constexpr std::size_t kNColumn = 3;
constexpr std::size_t kLowColumn = 5;   // 1st percentile
constexpr std::size_t kHighColumn = 9;  // 99th percentile

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
        for (const auto type : {
                 nmr1d::CarbonType::C,
                 nmr1d::CarbonType::CH,
                 nmr1d::CarbonType::CH2,
                 nmr1d::CarbonType::CH3,
             }) {
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
    const auto n = parse_count(fields[kNColumn]);
    if (!n || *n == 0) {
        return parse_error(
            line, "n must be a positive whole number, got '" + std::string{fields[kNColumn]} + "'");
    }
    std::vector<double> values;
    for (std::size_t column = kNColumn + 1; column < kColumns; ++column) {
        const auto value = io::parse_decimal(fields[column]);
        if (!value || !std::isfinite(value.value())) {
            return parse_error(line, "not a finite number: '" + std::string{fields[column]} + "'");
        }
        values.push_back(value.value());
    }
    const double low = values[kLowColumn - kNColumn - 1];
    const double high = values[kHighColumn - kNColumn - 1];
    if (low > high) {
        return parse_error(line, "1st percentile " + std::string{fields[kLowColumn]} +
                                     " is above the 99th " + std::string{fields[kHighColumn]});
    }
    return HallmarkWindow{
        .group = *group,
        .hallmark = std::string{fields[1]},
        .types = *types,
        .low = core::Ppm{low},
        .high = core::Ppm{high},
        .n = *n,
    };
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

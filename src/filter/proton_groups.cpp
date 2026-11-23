#include "arise/filter/proton_groups.hpp"

#include "arise/core/error.hpp"
#include "arise/io/number.hpp"
#include "group_proton_shifts_data.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <tuple>
#include <utility>

namespace arise::filter {
namespace {

constexpr std::size_t kFields = 21;

// Columns of data/group_proton_shifts.tsv.
constexpr std::size_t kGroup = 0;
constexpr std::size_t kHallmark = 1;
constexpr std::size_t kType = 2;
constexpr std::size_t kN = 3;
constexpr std::size_t kMolecules = 4;
constexpr std::size_t kLow = 12;
constexpr std::size_t kHigh = 13;
constexpr std::size_t kPrecision = 16;
constexpr std::size_t kBaseRate = 19;
constexpr std::size_t kRound = 20;

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
        .message = "group 1H windows, line " + std::to_string(line) + ": " + message,
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

// A protonated carbon type: CH, CH2 or CH3.
std::optional<nmr1d::CarbonType> protonated_type_named(std::string_view name) {
    for (const auto type :
         {nmr1d::CarbonType::CH, nmr1d::CarbonType::CH2, nmr1d::CarbonType::CH3}) {
        if (nmr1d::to_string(type) == name) {
            return type;
        }
    }
    return std::nullopt;
}

int hydrogens_of(nmr1d::CarbonType type) noexcept {
    return static_cast<int>(type);
}

std::optional<std::size_t> parse_count(std::string_view text) {
    if (text.empty() || text.size() > 9 ||
        !std::all_of(text.begin(), text.end(), [](char c) { return c >= '0' && c <= '9'; })) {
        return std::nullopt;
    }
    std::size_t value = 0;
    for (const char c : text) {
        value = (value * 10) + static_cast<std::size_t>(c - '0');
    }
    return value;
}

// The numbers of one row: every numeric column, in order; error naming the bad field.
core::Result<std::array<double, kFields>> parse_numbers(std::size_t line,
                                                        const std::vector<std::string_view>& f) {
    std::array<double, kFields> numbers{};
    for (std::size_t column = kN; column < kRound; ++column) {
        const auto value = io::parse_decimal(f[column]);
        if (!value.ok() || !std::isfinite(value.value())) {
            return parse_error(line, "column " + std::to_string(column + 1) + " '" +
                                         std::string{f[column]} + "' is not a number");
        }
        numbers[column] = value.value();
    }
    return numbers;
}

core::Result<ProtonWindow> parse_row(std::size_t line, std::string_view text) {
    const auto f = split(text, '\t');
    if (f.size() != kFields) {
        return parse_error(line, "expected 21 fields, found " + std::to_string(f.size()));
    }
    const auto group = group_named(f[kGroup]);
    if (!group.has_value()) {
        return parse_error(line, "unknown group '" + std::string{f[kGroup]} + "'");
    }
    const auto type = protonated_type_named(f[kType]);
    if (!type.has_value()) {
        return parse_error(line,
                           "H count must be CH, CH2 or CH3, not '" + std::string{f[kType]} + "'");
    }
    const auto n = parse_count(f[kN]);
    const auto molecules = parse_count(f[kMolecules]);
    if (!n.has_value() || !molecules.has_value()) {
        return parse_error(line, "n and molecules must be whole numbers");
    }
    const auto numbers = parse_numbers(line, f);
    if (!numbers.ok()) {
        return numbers.error();
    }
    const auto& v = numbers.value();
    if (v[kLow] > v[kHigh]) {
        return parse_error(
            line, "window low " + std::string{f[kLow]} + " is above high " + std::string{f[kHigh]});
    }
    if (f[kRound] != "1" && f[kRound] != "-") {
        return parse_error(line, "round must be '1' or '-', not '" + std::string{f[kRound]} + "'");
    }
    return ProtonWindow{
        .group = *group,
        .hallmark = std::string{f[kHallmark]},
        .type = *type,
        .low = core::Ppm{v[kLow]},
        .high = core::Ppm{v[kHigh]},
        .n = *n,
        .molecules = *molecules,
        .precision = v[kPrecision],
        .base_rate = v[kBaseRate],
        .round1 = f[kRound] == "1",
    };
}

bool matches(const io::ObservedPeak& peak, const ProtonWindow& window, ProtonOptions options) {
    if (!window.round1 || peak.shift < window.low || peak.shift > window.high) {
        return false;
    }
    const auto protons = proton_count(peak);
    return !options.use_counts || !protons.has_value() || *protons >= hydrogens_of(window.type);
}

}  // namespace

core::Result<std::vector<ProtonWindow>> parse_proton_windows(std::string_view text) {
    std::vector<ProtonWindow> windows;
    std::set<std::tuple<chem::Group, std::string, nmr1d::CarbonType>> seen;
    std::size_t number = 0;
    for (const auto line : split(text, '\n')) {
        ++number;
        if (line.empty() || line.front() == '#') {
            continue;
        }
        auto window = parse_row(number, line);
        if (!window.ok()) {
            return std::move(window).error();
        }
        if (!seen.emplace(window.value().group, window.value().hallmark, window.value().type)
                 .second) {
            return parse_error(number, "the window " +
                                           std::string{chem::to_string(window.value().group)} +
                                           " / " + window.value().hallmark + " is listed twice");
        }
        windows.push_back(std::move(window).value());
    }
    return windows;
}

core::Result<std::vector<ProtonWindow>> measured_proton_windows() {
    return parse_proton_windows(detail::kMeasuredGroupProtonShifts);
}

std::optional<int> proton_count(const io::ObservedPeak& peak) {
    if (!peak.integral.has_value()) {
        return std::nullopt;
    }
    // Rounded to the nearest whole number, halves up, at least 1: as hsqc counts the protons
    // of an exchangeable peak.
    const auto rounded = static_cast<int>(std::floor(*peak.integral + 0.5));
    return std::max(1, rounded);
}

std::vector<ProtonSuggestion> suggest_groups_from_protons(std::span<const io::ObservedPeak> peaks,
                                                          std::span<const ProtonWindow> windows,
                                                          ProtonOptions options) {
    std::vector<ProtonSuggestion> suggestions;
    for (const auto group : chem::kAllGroups) {
        ProtonSuggestion suggestion{.group = group, .evidence = {}};
        for (std::size_t p = 0; p < peaks.size(); ++p) {
            for (const auto& window : windows) {
                if (window.group == group && matches(peaks[p], window, options)) {
                    suggestion.evidence.push_back({
                        .peak = p,
                        .shift = peaks[p].shift,
                        .protons = proton_count(peaks[p]),
                        .window = window,
                    });
                }
            }
        }
        if (!suggestion.evidence.empty()) {
            suggestions.push_back(std::move(suggestion));
        }
    }
    return suggestions;
}

}  // namespace arise::filter

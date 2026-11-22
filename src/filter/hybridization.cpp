#include "arise/filter/hybridization.hpp"

#include "arise/core/error.hpp"
#include "arise/io/number.hpp"
#include "hybridization_shifts_data.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace arise::filter {

namespace {

// Columns of data/hybridization_shifts.tsv.
constexpr std::size_t kColumns = 10;
constexpr std::size_t kCarbonsColumn = 1;
constexpr std::size_t kMoleculesColumn = 2;
constexpr std::size_t kFirstShiftColumn = 3;  // min; then 0.1st, 1st, median, 99th, 99.9th, max
constexpr std::size_t kLowColumn = 4;         // 0.1st percentile
constexpr std::size_t kHighColumn = 8;        // 99.9th percentile

constexpr std::array<std::string_view, 3> kNames{"sp", "sp2", "sp3"};

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
        .message = "hybridization windows, line " + std::to_string(line) + ": " + message,
    };
}

bool positive_whole_number(std::string_view text) {
    return !text.empty() &&
           std::all_of(text.begin(), text.end(), [](char c) { return c >= '0' && c <= '9'; }) &&
           text.find_first_not_of('0') != std::string_view::npos;
}

// The range of one row, and its index into kNames.
struct Row {
    std::size_t index = 0;
    ShiftRange range;
};

core::Result<Row> parse_row(std::size_t line, std::string_view text) {
    const auto fields = split(text, '\t');
    if (fields.size() != kColumns) {
        return parse_error(line, "expected " + std::to_string(kColumns) + " fields, got " +
                                     std::to_string(fields.size()));
    }
    const auto* name = std::find(kNames.begin(), kNames.end(), fields[0]);
    if (name == kNames.end()) {
        return parse_error(line, "unknown hybridization '" + std::string{fields[0]} + "'");
    }
    for (const auto column : {kCarbonsColumn, kMoleculesColumn}) {
        if (!positive_whole_number(fields[column])) {
            return parse_error(line, "count must be a positive whole number, got '" +
                                         std::string{fields[column]} + "'");
        }
    }
    std::array<double, kColumns> values{};
    for (std::size_t column = kFirstShiftColumn; column < kColumns; ++column) {
        const auto value = io::parse_decimal(fields[column]);
        if (!value || !std::isfinite(value.value())) {
            return parse_error(line, "not a finite number: '" + std::string{fields[column]} + "'");
        }
        values[column] = value.value();
    }
    if (values[kLowColumn] > values[kHighColumn]) {
        return parse_error(line, "0.1st percentile " + std::string{fields[kLowColumn]} +
                                     " is above the 99.9th " + std::string{fields[kHighColumn]});
    }
    return Row{
        .index = static_cast<std::size_t>(name - kNames.begin()),
        .range = {.low = core::Ppm{values[kLowColumn]}, .high = core::Ppm{values[kHighColumn]}},
    };
}

// Shortest fixed-point form: 29.000000 -> "29", 45.500000 -> "45.5".
std::string format_ppm(core::Ppm shift) {
    std::string text = std::to_string(shift.value());
    text.erase(text.find_last_not_of('0') + 1);
    if (text.back() == '.') {
        text.pop_back();
    }
    return text;
}

bool inside(core::Ppm shift, const ShiftRange& range) noexcept {
    return range.low <= shift && shift <= range.high;
}

std::string plural(int count, std::string_view one, std::string_view many) {
    return std::to_string(count) + " " + std::string{count == 1 ? one : many};
}

}  // namespace

core::Result<HybridizationWindows> parse_hybridization_windows(std::string_view text) {
    std::array<ShiftRange, kNames.size()> ranges{};
    std::array<bool, kNames.size()> seen{};
    std::size_t line = 0;
    for (const auto raw : split(text, '\n')) {
        ++line;
        if (raw.empty() || raw.front() == '#') {
            continue;
        }
        const auto row = parse_row(line, raw);
        if (!row) {
            return row.error();
        }
        const auto index = row.value().index;
        if (seen[index]) {
            return parse_error(line, std::string{kNames[index]} + " listed twice");
        }
        seen[index] = true;
        ranges[index] = row.value().range;
    }
    for (std::size_t i = 0; i < kNames.size(); ++i) {
        if (!seen[i]) {
            return core::Error{
                .code = core::ErrorCode::ParseError,
                .message = "hybridization windows: no row for " + std::string{kNames[i]},
            };
        }
    }
    return HybridizationWindows{.sp = ranges[0], .sp2 = ranges[1], .sp3 = ranges[2]};
}

core::Result<HybridizationWindows> measured_hybridization_windows() {
    return parse_hybridization_windows(detail::kMeasuredHybridizationShifts);
}

core::Result<HybridizationSuggestion> suggest_hybridization(
    std::span<const nmr1d::CarbonAssignment> peaks, const HybridizationWindows& windows,
    std::optional<int> typed_sp_peaks) {
    if (typed_sp_peaks && *typed_sp_peaks < 0) {
        return core::Error{
            .code = core::ErrorCode::InvalidArgument,
            .message = "suggest_hybridization: the typed sp count is negative",
        };
    }
    HybridizationSuggestion suggestion;
    int only_sp2 = 0;
    int only_sp3 = 0;
    int outside = 0;
    for (const auto& peak : peaks) {
        const PeakHybridization options{
            .shift = peak.shift,
            .could_be_sp2 = inside(peak.shift, windows.sp2),
            .could_be_sp3 = inside(peak.shift, windows.sp3),
        };
        suggestion.sp2.max += options.could_be_sp2 ? 1 : 0;
        suggestion.sp3.max += options.could_be_sp3 ? 1 : 0;
        only_sp2 += options.could_be_sp2 && !options.could_be_sp3 ? 1 : 0;
        only_sp3 += options.could_be_sp3 && !options.could_be_sp2 ? 1 : 0;
        outside += !options.could_be_sp2 && !options.could_be_sp3 ? 1 : 0;
        suggestion.peaks.push_back(options);
    }
    const int sp = typed_sp_peaks.value_or(0);
    suggestion.sp2.min = std::max(0, only_sp2 - sp);
    suggestion.sp3.min = std::max(0, only_sp3 - sp);
    suggestion.notes.push_back(
        "sp carbons appear from " + format_ppm(windows.sp.low) + " to " +
        format_ppm(windows.sp.high) +
        " ppm in nmrshiftdb2, overlapping both ranges, so ARISE does not suggest their count; "
        "type it if known");
    if (sp > 0) {
        suggestion.notes.push_back("minimums lowered by the " +
                                   plural(sp, "typed sp peak", "typed sp peaks"));
    }
    if (outside > 0) {
        suggestion.notes.push_back(plural(outside, "peak lies", "peaks lie") +
                                   " outside both the sp2 and the sp3 range");
    }
    return suggestion;
}

}  // namespace arise::filter

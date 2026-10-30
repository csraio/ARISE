#include "arise/io/peak_csv.hpp"

#include "arise/core/error.hpp"
#include "arise/io/number.hpp"

#include <cstddef>
#include <span>
#include <string>
#include <utility>

namespace arise::io {
namespace {

constexpr std::string_view kShiftColumn = "shift_ppm";
constexpr std::string_view kIntegralColumn = "integral";
constexpr std::string_view kByteOrderMark = "\xEF\xBB\xBF";

// A line that holds data or the header, with its 1-based physical line number.
struct Line {
    std::size_t number;
    std::string_view text;
};

// Column positions resolved from the header.
struct Columns {
    std::size_t count;
    std::size_t shift;
    std::optional<std::size_t> integral;
};

core::Error error_at(core::ErrorCode code, std::size_t line, std::string_view what) {
    return core::Error{
        .code = code,
        .message = "line " + std::to_string(line) + ": " + std::string{what},
    };
}

std::string quoted(std::string_view text) {
    return "'" + std::string{text} + "'";
}

std::string_view trim(std::string_view text) noexcept {
    constexpr std::string_view kBlank = " \t";
    const auto first = text.find_first_not_of(kBlank);
    if (first == std::string_view::npos) {
        return {};
    }
    text.remove_suffix(text.size() - 1 - text.find_last_not_of(kBlank));
    text.remove_prefix(first);
    return text;
}

// Splits on LF (dropping one trailing CR) and keeps only header and data lines:
// blank lines and `#` comments are skipped but still counted.
std::vector<Line> content_lines(std::string_view text) {
    std::vector<Line> lines;
    std::size_t number = 0;
    while (!text.empty()) {
        ++number;
        const auto newline = text.find('\n');
        std::string_view raw = text.substr(0, newline);
        text = newline == std::string_view::npos ? std::string_view{} : text.substr(newline + 1);
        if (!raw.empty() && raw.back() == '\r') {
            raw.remove_suffix(1);
        }
        const auto content = trim(raw);
        if (content.empty() || content.front() == '#') {
            continue;
        }
        lines.push_back(Line{.number = number, .text = raw});
    }
    return lines;
}

std::vector<std::string_view> split_fields(std::string_view line) {
    std::vector<std::string_view> fields;
    std::size_t start = 0;
    while (true) {
        const auto comma = line.find(',', start);
        fields.push_back(trim(line.substr(start, comma - start)));
        if (comma == std::string_view::npos) {
            return fields;
        }
        start = comma + 1;
    }
}

core::Result<Columns> parse_header(const Line& header) {
    const auto names = split_fields(header.text);
    std::optional<std::size_t> shift;
    std::optional<std::size_t> integral;
    for (std::size_t i = 0; i < names.size(); ++i) {
        const auto name = names[i];
        std::optional<std::size_t>* slot = nullptr;
        if (name == kShiftColumn) {
            slot = &shift;
        } else if (name == kIntegralColumn) {
            slot = &integral;
        } else if (name.empty()) {
            return error_at(core::ErrorCode::ParseError, header.number,
                            "empty column name in column " + std::to_string(i + 1));
        } else {
            return error_at(core::ErrorCode::ParseError, header.number,
                            "unknown column " + quoted(name) +
                                " (allowed: " + std::string{kShiftColumn} + ", " +
                                std::string{kIntegralColumn} + ")");
        }
        if (slot->has_value()) {
            return error_at(core::ErrorCode::ParseError, header.number,
                            "duplicate column " + quoted(name));
        }
        *slot = i;
    }
    if (!shift.has_value()) {
        return error_at(core::ErrorCode::ParseError, header.number,
                        "missing required column " + quoted(kShiftColumn));
    }
    return Columns{.count = names.size(), .shift = *shift, .integral = integral};
}

core::Result<double> parse_field(const Line& line, std::string_view column,
                                 std::string_view field) {
    if (field.empty()) {
        return error_at(core::ErrorCode::ParseError, line.number,
                        "empty " + std::string{column} + " field");
    }
    if (field.find('"') != std::string_view::npos) {
        return error_at(
            core::ErrorCode::ParseError, line.number,
            std::string{column} + " field " + quoted(field) + ": quoted fields are not supported");
    }
    auto value = parse_decimal(field);
    if (!value.ok()) {
        return error_at(core::ErrorCode::ParseError, line.number,
                        std::string{column} + ": " + value.error().message);
    }
    return value;
}

core::Result<ObservedPeak> parse_row(const Line& line, const Columns& columns) {
    const auto fields = split_fields(line.text);
    if (fields.size() != columns.count) {
        return error_at(core::ErrorCode::ParseError, line.number,
                        "expected " + std::to_string(columns.count) +
                            (columns.count == 1 ? " field" : " fields") + ", got " +
                            std::to_string(fields.size()));
    }

    const auto shift_text = fields[columns.shift];
    auto shift_value = parse_field(line, kShiftColumn, shift_text);
    if (!shift_value.ok()) {
        return std::move(shift_value).error();
    }
    const core::Ppm shift{shift_value.value()};
    if (shift < kMinProtonShift || shift > kMaxProtonShift) {
        return error_at(core::ErrorCode::OutOfRange, line.number,
                        std::string{kShiftColumn} + " " + std::string{shift_text} +
                            " is outside the plausible 1H window [" +
                            std::to_string(kMinProtonShift.value()) + ", " +
                            std::to_string(kMaxProtonShift.value()) +
                            "] ppm; is it a 13C shift or a value in Hz?");
    }

    ObservedPeak peak{.shift = shift, .integral = std::nullopt};
    if (columns.integral.has_value()) {
        const auto integral_text = fields[*columns.integral];
        auto integral = parse_field(line, kIntegralColumn, integral_text);
        if (!integral.ok()) {
            return std::move(integral).error();
        }
        if (integral.value() <= 0.0) {
            return error_at(
                core::ErrorCode::OutOfRange, line.number,
                std::string{kIntegralColumn} + " must be > 0, got " + std::string{integral_text});
        }
        peak.integral = integral.value();
    }
    return peak;
}

}  // namespace

core::Result<std::vector<ObservedPeak>> parse_peak_csv(std::string_view text) {
    if (text.starts_with(kByteOrderMark)) {
        text.remove_prefix(kByteOrderMark.size());
    }
    const auto lines = content_lines(text);
    if (lines.empty()) {
        return core::Error{
            .code = core::ErrorCode::ParseError,
            .message = "missing header line (expected a " + quoted(kShiftColumn) + " column)",
        };
    }

    auto columns = parse_header(lines.front());
    if (!columns.ok()) {
        return std::move(columns).error();
    }

    std::vector<ObservedPeak> peaks;
    peaks.reserve(lines.size() - 1);
    for (const Line& line : std::span{lines}.subspan(1)) {
        auto peak = parse_row(line, columns.value());
        if (!peak.ok()) {
            return std::move(peak).error();
        }
        peaks.push_back(peak.value());
    }
    if (peaks.empty()) {
        return error_at(core::ErrorCode::ParseError, lines.front().number,
                        "no peaks after the header");
    }
    return peaks;
}

}  // namespace arise::io

#include "arise/io/pasted_peaks.hpp"

#include "arise/io/number.hpp"
#include "carbon_rows.hpp"
#include "csv_table.hpp"
#include "text_lines.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace arise::io {
namespace {

// Lines with content, numbered from 1. Unlike detail::content_lines, a lone CR also ends a
// line: pasted text can come from software that writes old Mac line ends.
std::vector<detail::Line> pasted_lines(std::string_view text) {
    constexpr std::string_view kByteOrderMark = "\xEF\xBB\xBF";
    if (text.starts_with(kByteOrderMark)) {
        text.remove_prefix(kByteOrderMark.size());
    }
    std::vector<detail::Line> lines;
    std::size_t number = 0;
    while (!text.empty()) {
        ++number;
        const auto end = text.find_first_of("\r\n");
        const auto raw = text.substr(0, end);
        if (end == std::string_view::npos) {
            text = {};
        } else {
            const bool crlf = text[end] == '\r' && end + 1 < text.size() && text[end + 1] == '\n';
            text.remove_prefix(end + (crlf ? 2 : 1));
        }
        const auto content = detail::trim(raw);
        if (!content.empty() && content.front() != '#') {
            lines.push_back({.number = number, .text = content});
        }
    }
    return lines;
}

// The fields of a line: runs of anything but tabs, semicolons and spaces.
std::vector<std::string_view> fields_of(std::string_view line) {
    constexpr std::string_view kSeparators = " \t;";
    std::vector<std::string_view> fields;
    std::size_t start = line.find_first_not_of(kSeparators);
    while (start != std::string_view::npos) {
        const auto end = line.find_first_of(kSeparators, start);
        fields.push_back(line.substr(start, end - start));
        start = end == std::string_view::npos ? end : line.find_first_not_of(kSeparators, end);
    }
    return fields;
}

// Whether `field` starts like a number: a digit, sign or decimal point ("29,5" does).
bool starts_like_a_number(std::string_view field) {
    const char c = field.front();  // fields are never empty
    return (c >= '0' && c <= '9') || c == '+' || c == '-' || c == '.';
}

bool is_letter(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

// Whether `line`, the first with content, is a header: it has a letter and no field starts
// like a number, so a mistyped first peak ("29,5") is an error, never skipped.
bool is_header(std::string_view line) {
    const auto fields = fields_of(line);
    return std::any_of(line.begin(), line.end(), is_letter) &&
           std::none_of(fields.begin(), fields.end(), starts_like_a_number);
}

// Whether a line must, or may, carry an intensity after its shift.
enum class Intensity : std::uint8_t { Optional, Required };

// The table row of one data line: a shift and, as `intensity` allows, an intensity.
core::Result<detail::TableRow> read_row(const detail::Line& line, Intensity intensity) {
    const auto fields = fields_of(line.text);
    const std::size_t most = 2;
    const std::size_t least = intensity == Intensity::Required ? 2 : 1;
    if (fields.size() < least || fields.size() > most) {
        const std::string expected = intensity == Intensity::Required
                                         ? "a shift and an intensity"
                                         : "a shift, optionally followed by an intensity";
        return detail::error_at(core::ErrorCode::ParseError, line.number,
                                "expected " + expected + ", found " +
                                    std::to_string(fields.size()) + " fields in '" +
                                    std::string{line.text} + "'");
    }
    detail::TableRow row{
        .line = line.number,
        .values = {std::nullopt, std::nullopt},
        .texts = {{}, {}},
    };
    constexpr std::array<std::string_view, 2> kNames{"shift", "intensity"};
    for (std::size_t i = 0; i < fields.size(); ++i) {
        auto value = parse_decimal(fields[i]);
        if (!value.ok()) {
            return detail::error_at(core::ErrorCode::ParseError, line.number,
                                    std::string{kNames[i]} + ": " + value.error().message);
        }
        row.values[i] = value.value();
        row.texts[i] = fields[i];
    }
    return row;
}

// Reads every data line of `text` into a peak with `convert`.
template <typename Peak, typename Convert>
core::Result<std::vector<Peak>> read_pasted(std::string_view text, Intensity intensity,
                                            Convert convert) {
    auto lines = pasted_lines(text);
    if (!lines.empty() && is_header(lines.front().text)) {
        lines.erase(lines.begin());
    }
    std::vector<detail::TableRow> rows;
    rows.reserve(lines.size());
    for (const auto& line : lines) {
        auto row = read_row(line, intensity);
        if (!row.ok()) {
            return std::move(row).error();
        }
        rows.push_back(std::move(row).value());
    }
    return detail::convert_rows<Peak>(std::move(rows), convert);
}

}  // namespace

core::Result<std::vector<CarbonPeak>> parse_pasted_carbon(std::string_view text) {
    auto peaks = read_pasted<CarbonPeak>(text, Intensity::Optional, detail::to_carbon_peak);
    if (peaks.ok() && peaks.value().empty()) {
        return core::Error{
            .code = core::ErrorCode::ParseError,
            .message = "no 13C peaks: paste one shift per line",
        };
    }
    return peaks;
}

core::Result<std::vector<DeptPeak>> parse_pasted_dept(std::string_view text) {
    return read_pasted<DeptPeak>(text, Intensity::Required, detail::to_dept_peak);
}

}  // namespace arise::io

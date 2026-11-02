#pragma once
// Internal to arise_io: the numeric CSV table format shared by every peak-list reader.

#include "arise/core/error.hpp"
#include "arise/core/result.hpp"
#include "arise/core/units.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace arise::io::detail {

/// What a column holds: a number (parsed with `parse_decimal`) or free text.
enum class ColumnKind : std::uint8_t { Number, Text };

/// One column a reader understands.
struct ColumnSpec {
    std::string_view name;
    bool required;
    ColumnKind kind = ColumnKind::Number;
};

/// Whether a table with a header but no data rows is valid.
enum class EmptyTable : std::uint8_t { Reject, Allow };

/// A data row. `values[i]` and `texts[i]` belong to `ColumnSpec` i and are empty when
/// that optional column is absent from the file (a text column has no value; an optional
/// text field may also be blank). Texts view the caller's input.
struct TableRow {
    std::size_t line;
    std::vector<std::optional<double>> values;
    std::vector<std::string_view> texts;
};

/// Reads a numeric table in the shared peak-list format:
/// - UTF-8, optional byte-order mark, LF or CRLF line endings.
/// - Blank lines and `#` comments are skipped but counted for line numbers.
/// - The first remaining line is the header; its names must be known columns, without
///   duplicates, and must include every required column.
/// - Each data row has one field per header column; fields are trimmed of spaces and
///   tabs and unquoted. Number fields are parsed with `parse_decimal`; text fields are
///   kept as written. Only optional text fields may be blank.
/// - At least one data row, unless `empty` is EmptyTable::Allow; the error calls the
///   rows `rows_name` ("no peaks after the header").
/// Errors are ParseError and name the 1-based physical line.
[[nodiscard]] core::Result<std::vector<TableRow>> read_numeric_table(
    std::string_view text, std::span<const ColumnSpec> columns,
    EmptyTable empty = EmptyTable::Reject, std::string_view rows_name = "peaks");

/// A plausibility window for one nucleus, with a hint for values outside it.
struct ShiftWindow {
    std::string_view nucleus;
    core::Ppm min;
    core::Ppm max;
    std::string_view hint;
};

/// The value of required column `column` (named `name`) as a shift inside `window`;
/// otherwise OutOfRange naming the line, the column, the original text and the hint.
[[nodiscard]] core::Result<core::Ppm> shift_in_window(const TableRow& row, std::size_t column,
                                                      std::string_view name,
                                                      const ShiftWindow& window);

// Converts every row with `convert`, stopping at the first error.
template <typename Peak, typename Convert>
core::Result<std::vector<Peak>> convert_rows(core::Result<std::vector<TableRow>> rows,
                                             Convert convert) {
    if (!rows.ok()) {
        return std::move(rows).error();
    }
    std::vector<Peak> peaks;
    peaks.reserve(rows.value().size());
    for (const auto& row : rows.value()) {
        auto peak = convert(row);
        if (!peak.ok()) {
            return std::move(peak).error();
        }
        peaks.push_back(peak.value());
    }
    return peaks;
}

/// "line N: what", with the given code.
[[nodiscard]] core::Error error_at(core::ErrorCode code, std::size_t line, std::string_view what);

}  // namespace arise::io::detail

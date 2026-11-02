#pragma once

#include "arise/core/result.hpp"
#include "arise/core/units.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace arise::io {

/// One row of a lab's 13C suggestion-window table, as written. Names are checked later,
/// by `filter`, which knows the groups; this reader checks the file's shape and numbers.
struct WindowRow {
    /// 1-based physical line, for error messages.
    std::size_t line = 0;
    std::string group;
    /// Empty when the column is absent or the field is blank.
    std::string hallmark;
    /// DEPT types as written (e.g. "C CH"); empty when absent or blank (any type).
    std::string dept;
    core::Ppm low;
    core::Ppm high;

    friend bool operator==(const WindowRow&, const WindowRow&) = default;
};

/// Parses a lab's suggestion-window table: the shared peak-list format (see
/// parse_peak_csv; comment lines, byte-order mark, CRLF, header names in any order) with
/// columns `group`, `low_ppm` and `high_ppm` (required) and `hallmark` and `dept`
/// (optional, may be blank). `group` must not be blank. Shifts must lie in the 13C window
/// [kMinCarbonShift, kMaxCarbonShift], and `low_ppm` must not exceed `high_ppm`
/// (OutOfRange otherwise). Rows are returned in file order; at least one is required.
[[nodiscard]] core::Result<std::vector<WindowRow>> parse_window_csv(std::string_view text);

}  // namespace arise::io

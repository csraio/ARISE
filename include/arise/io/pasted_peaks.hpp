#pragma once

#include "arise/core/result.hpp"
#include "arise/io/carbon_csv.hpp"
#include "arise/io/peak_csv.hpp"

#include <string_view>
#include <vector>

namespace arise::io {

/// Peak lists pasted by the chemist, e.g. columns copied from a spreadsheet or from the
/// spectrometer software's peak table. The same peaks the CSV readers return, with the
/// same value checks; only the text layout differs.
///
/// - Lines end in LF, CRLF or CR. Blank lines and lines starting with `#` are ignored.
/// - A line holds a shift, then optionally an intensity, separated by tabs, semicolons
///   or spaces. A comma is not a separator, so a comma decimal ("29,5") is an error
///   instead of being read as two numbers.
/// - The first non-blank line may be a header: it is skipped when it has a letter and none
///   of its fields starts like a number (a digit, sign or decimal point), e.g.
///   "ppm<TAB>Intensity". So a mistyped first peak ("29,5") is an error, never skipped.
///   Any later line that does not read is an error.
/// - Numbers are read with `parse_decimal`.
///
/// Errors name the 1-based line: ParseError for text that does not read (a third column,
/// a missing intensity), OutOfRange for values outside the CSV readers' limits.

/// 13C peaks: the intensity is optional on each line. At least one peak is required.
[[nodiscard]] core::Result<std::vector<CarbonPeak>> parse_pasted_carbon(std::string_view text);

/// DEPT peaks (DEPT-90 or DEPT-135): every line needs an intensity, whose sign is the
/// phase, so it must be nonzero. Text without peaks is a valid, empty spectrum.
[[nodiscard]] core::Result<std::vector<DeptPeak>> parse_pasted_dept(std::string_view text);

/// 1H peaks: a shift, then optionally the peak's integral as a number of protons (> 0), with
/// parse_peak_csv's limits. At least one peak is required.
[[nodiscard]] core::Result<std::vector<ObservedPeak>> parse_pasted_protons(std::string_view text);

}  // namespace arise::io

#pragma once

#include "arise/core/result.hpp"
#include "arise/core/units.hpp"

#include <optional>
#include <string_view>
#include <vector>

namespace arise::io {

/// One observed 1H resonance.
struct ObservedPeak {
    core::Ppm shift;
    /// Relative integral (proton count); empty when the file has no `integral` column.
    std::optional<double> integral;

    friend bool operator==(const ObservedPeak&, const ObservedPeak&) = default;
};

/// Plausibility window for 1H shifts, deliberately wide. Strongly shielded protons do
/// reach negative values (the internal methyls of trans-15,16-dimethyldihydropyrene sit
/// near -4.25 ppm; Boekelheide), and acidic or chelated OH can exceed 15 ppm. Values
/// outside this window almost always mean the wrong column: 13C shifts or Hz offsets.
inline constexpr core::Ppm kMinProtonShift{-10.0};
inline constexpr core::Ppm kMaxProtonShift{25.0};

/// Parses an observed 1H peak list (format version 1).
///
/// - UTF-8 text; a leading byte-order mark is skipped. Lines end in LF or CRLF.
/// - Blank lines and lines whose first non-blank character is `#` are ignored.
/// - The first remaining line is the header: comma-separated column names.
///   `shift_ppm` is required, `integral` is optional; any order, no duplicates,
///   no other names. Every data row has exactly one field per column.
/// - Fields are trimmed of spaces and tabs and parsed with `parse_decimal`; quoting is
///   not supported.
/// - Shifts must lie in [kMinProtonShift, kMaxProtonShift]; integrals must be > 0.
/// - Peaks are returned in file order. At least one peak is required.
///
/// Errors name the 1-based physical line: ParseError for malformed text, OutOfRange for
/// well-formed values outside their allowed range.
[[nodiscard]] core::Result<std::vector<ObservedPeak>> parse_peak_csv(std::string_view text);

}  // namespace arise::io

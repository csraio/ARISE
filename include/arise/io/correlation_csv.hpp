#pragma once

#include "arise/core/result.hpp"
#include "arise/core/units.hpp"

#include <optional>
#include <string_view>
#include <vector>

namespace arise::io {

/// One cross peak of a 1H-13C correlation spectrum (HSQC or HMBC).
struct CorrelationPeak {
    core::Ppm proton;
    core::Ppm carbon;
    /// Peak volume or height as exported; empty when the file has no `intensity` column.
    std::optional<double> intensity;

    friend bool operator==(const CorrelationPeak&, const CorrelationPeak&) = default;
};

/// Parses an HSQC peak list (one-bond C-H correlations): the shared peak-list format
/// with columns `h_ppm` and `c_ppm` (required) and `intensity` (optional).
/// - `h_ppm` must lie in the 1H window and `c_ppm` in the 13C window (OutOfRange
///   otherwise, with a hint that the two columns may be swapped).
/// - In a multiplicity-edited HSQC the sign of `intensity` separates CH2 from CH/CH3,
///   so a given intensity must be nonzero.
/// - At least one peak; leave the file out entirely when no HSQC was acquired.
[[nodiscard]] core::Result<std::vector<CorrelationPeak>> parse_hsqc_csv(std::string_view text);

/// Parses an HMBC peak list (two- and three-bond C-H correlations): same columns and
/// shift windows as parse_hsqc_csv. `intensity` may be any finite value (ARISE does not
/// use its sign). At least one peak.
[[nodiscard]] core::Result<std::vector<CorrelationPeak>> parse_hmbc_csv(std::string_view text);

}  // namespace arise::io

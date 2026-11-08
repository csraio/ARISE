#pragma once

#include "arise/core/result.hpp"
#include "arise/core/units.hpp"

#include <optional>
#include <string_view>
#include <vector>

namespace arise::io {

/// One peak of a broadband-decoupled 13C spectrum.
struct CarbonPeak {
    core::Ppm shift;
    /// Peak height or area as exported; empty when the file has no `intensity` column.
    std::optional<double> intensity;

    friend bool operator==(const CarbonPeak&, const CarbonPeak&) = default;
};

/// One peak of a DEPT spectrum (DEPT-90 or DEPT-135). The sign of `intensity` is the
/// phase: in DEPT-135, CH and CH3 point up (positive) and CH2 points down (negative).
struct DeptPeak {
    core::Ppm shift;
    double intensity{};

    friend bool operator==(const DeptPeak&, const DeptPeak&) = default;
};

/// Plausibility window for 13C shifts: a convention, deliberately wide. It must contain
/// every 13C shift nmrshiftdb2 lists for an in-scope molecule (measured by
/// scripts/build_reference_sets.py into tests/data/reference/measurements.tsv and checked
/// by a test). A value outside it usually means the wrong column, such as Hz offsets.
inline constexpr core::Ppm kMinCarbonShift{-50.0};
inline constexpr core::Ppm kMaxCarbonShift{300.0};

/// Parses a 13C peak list: the shared peak-list format (see parse_peak_csv) with
/// columns `shift_ppm` (required) and `intensity` (optional, any finite value).
/// Shifts must lie in [kMinCarbonShift, kMaxCarbonShift] (OutOfRange otherwise).
[[nodiscard]] core::Result<std::vector<CarbonPeak>> parse_carbon_csv(std::string_view text);

/// Parses a DEPT peak list: the shared peak-list format with columns `shift_ppm` and
/// `intensity`, both required. The intensity carries the phase, so it must be nonzero;
/// shifts must lie in the 13C window. Out-of-range values are OutOfRange. A header with
/// no rows is a valid, empty spectrum (a DEPT-90 of a molecule without CH carbons).
[[nodiscard]] core::Result<std::vector<DeptPeak>> parse_dept_csv(std::string_view text);

}  // namespace arise::io

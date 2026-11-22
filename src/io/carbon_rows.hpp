#pragma once
// Internal to arise_io: turning a numeric table row into a 13C or DEPT peak, shared by
// the CSV readers (carbon_csv.cpp) and the paste readers (pasted_peaks.cpp), so both
// apply the same value checks and messages.

#include "arise/io/carbon_csv.hpp"
#include "csv_table.hpp"

#include <cstddef>

namespace arise::io::detail {

// Indices into a peak row's values and texts.
inline constexpr std::size_t kShiftColumn = 0;
inline constexpr std::size_t kIntensityColumn = 1;

/// The 13C peak of `row` (shift required, intensity optional). OutOfRange for a shift
/// outside [kMinCarbonShift, kMaxCarbonShift].
[[nodiscard]] core::Result<CarbonPeak> to_carbon_peak(const TableRow& row);

/// The DEPT peak of `row` (shift and intensity required). OutOfRange for a shift outside
/// the 13C window or a zero intensity (no phase).
[[nodiscard]] core::Result<DeptPeak> to_dept_peak(const TableRow& row);

}  // namespace arise::io::detail

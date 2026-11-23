#pragma once
// Internal to arise_io: turning a numeric table row into a 1H peak, shared by the CSV reader
// (peak_csv.cpp) and the paste reader (pasted_peaks.cpp), so both apply the same value checks
// and messages. Row layout: shift (column 0), optional integral (column 1).

#include "arise/io/peak_csv.hpp"
#include "csv_table.hpp"

namespace arise::io::detail {

/// The 1H peak of `row`. OutOfRange for a shift outside [kMinProtonShift, kMaxProtonShift] or
/// an integral that is not > 0.
[[nodiscard]] core::Result<ObservedPeak> to_proton_peak(const TableRow& row);

}  // namespace arise::io::detail

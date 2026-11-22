#pragma once

#include "arise/core/result.hpp"
#include "arise/core/units.hpp"
#include "arise/nmr1d/dept.hpp"

#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

// sp2/sp3 peak counts for spectra without DEPT (plan 3c-b): ARISE suggests, the chemist
// confirms, and only confirmed counts change a candidate's tier (filter::check_carbon_fit).
namespace arise::filter {

/// An inclusive 13C shift range.
struct ShiftRange {
    core::Ppm low;
    core::Ppm high;

    friend bool operator==(const ShiftRange&, const ShiftRange&) = default;
};

/// The measured shift range of each carbon hybridization, from the 0.1st to the 99.9th
/// percentile (the chemist's decision, 2026-11-22). sp2 and sp3 suggest counts; sp is only
/// reported, because sp carbons overlap both, so their count is never suggested (the
/// chemist's decision, 2026-11-22).
struct HybridizationWindows {
    ShiftRange sp;
    ShiftRange sp2;
    ShiftRange sp3;

    friend bool operator==(const HybridizationWindows&, const HybridizationWindows&) = default;
};

/// Reads a table in the format of data/hybridization_shifts.tsv (written by
/// tools/measure_hybridization_shifts): '#' comment lines, then tab-separated rows of
/// hybridization ("sp", "sp2", "sp3"), n, molecules, min, 0.1st, 1st percentile, median,
/// 99th, 99.9th percentile, max. Each range runs from the 0.1st to the 99.9th percentile.
///
/// Errors (ParseError, citing the line): a row without 10 fields, an unknown or repeated
/// hybridization, a count that is not a positive whole number, a number that is not a plain
/// decimal or not finite, a 0.1st percentile above the 99.9th; and (ParseError) a table
/// lacking any of sp, sp2 and sp3.
[[nodiscard]] core::Result<HybridizationWindows> parse_hybridization_windows(std::string_view text);

/// The ranges ARISE ships with: data/hybridization_shifts.tsv, compiled in at configure time
/// and read with parse_hybridization_windows. Contains information from nmrshiftdb2
/// (www.nmrshiftdb.org), made available under the nmrshiftdb2 Database License.
[[nodiscard]] core::Result<HybridizationWindows> measured_hybridization_windows();

/// An inclusive range of counts.
struct CountRange {
    int min = 0;
    int max = 0;

    friend bool operator==(const CountRange&, const CountRange&) = default;
};

/// What its shift alone says about one peak's carbon.
struct PeakHybridization {
    core::Ppm shift;
    bool could_be_sp2 = false;  // inside the sp2 range
    bool could_be_sp3 = false;  // inside the sp3 range

    friend bool operator==(const PeakHybridization&, const PeakHybridization&) = default;
};

/// Suggested sp2 and sp3 peak counts, for the chemist to confirm or change.
struct HybridizationSuggestion {
    std::vector<PeakHybridization> peaks;  // in the order given
    CountRange sp2;
    CountRange sp3;
    /// The sp range (always), how the typed sp count lowered the minimums, and peaks
    /// outside both ranges.
    std::vector<std::string> notes;

    friend bool operator==(const HybridizationSuggestion&,
                           const HybridizationSuggestion&) = default;
};

/// Suggests how many 13C peaks come from sp2 and from sp3 carbons. Per type, `max` counts
/// the peaks inside its range and `min` the peaks inside its range only, less
/// `typed_sp_peaks` (the chemist's own sp count, if typed: an sp peak may sit in either
/// range), never below 0. A peak outside both ranges counts toward neither. Suggestions
/// change nothing by themselves. InvalidArgument if `typed_sp_peaks` is negative.
[[nodiscard]] core::Result<HybridizationSuggestion> suggest_hybridization(
    std::span<const nmr1d::CarbonAssignment> peaks, const HybridizationWindows& windows,
    std::optional<int> typed_sp_peaks = std::nullopt);

/// Minimum peak counts the chemist has confirmed (std::nullopt: not confirmed). The sp2 and
/// sp3 counts may come from suggest_hybridization; the sp count is only ever typed; the
/// DEPT-type counts are typed when the chemist knows them without a DEPT list.
struct ConfirmedCarbonCounts {
    std::optional<int> sp;
    std::optional<int> sp2;
    std::optional<int> sp3;
    std::optional<int> c;
    std::optional<int> ch;
    std::optional<int> ch2;
    std::optional<int> ch3;

    friend bool operator==(const ConfirmedCarbonCounts&, const ConfirmedCarbonCounts&) = default;
};

}  // namespace arise::filter

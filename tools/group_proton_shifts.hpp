#pragma once

#include "arise/chem/perception.hpp"
#include "arise/nmr1d/dept.hpp"
#include "group_shifts.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

// Development-only (never part of the engine): measures each group's 1H windows over
// nmrshiftdb2's assigned 1H spectra (tests/data/reference/nmrshiftdb2_1h_all.tsv, the same
// layout as the 13C extract: per atom, the mean 1H shift of the H on it), with ARISE's own
// key atoms (key_carbons). The 1H plan (2026-11-22): windows split by the key carbon's number
// of H, and each window's precision measured, because one 1H peak in a crowded range says
// little.
namespace arise::tools {

/// The 1H window of the H on one hallmark carbon type of one group, with how often a peak in
/// it means the group, over the even-ID records.
struct ProtonShifts {
    chem::Group group = chem::Group::CarboxylicAcid;
    std::string hallmark;  // as key_carbons names it ("O-C carbon (sp3)", ...)
    nmr1d::CarbonType type = nmr1d::CarbonType::CH;  // CH, CH2 or CH3: the H count
    std::size_t n = 0;                               // carbons (each key atom of each record once)
    std::size_t molecules = 0;                       // distinct nmrshiftdb2 IDs among them
    double min = 0;
    double p1 = 0;  // percentiles as Python statistics.quantiles(n=100, method="inclusive")
    double p5 = 0;
    double median = 0;
    double p95 = 0;
    double p99 = 0;
    double max = 0;
    /// The suggestion window: p1-p99, or min-max from fewer than
    /// filter::kMinMoleculesForPercentiles molecules (the 13C rule, the chemist's choice).
    double low = 0;
    double high = 0;
    /// Even-ID molecules with an H on a carbon of `type` (any carbon) inside the window, and
    /// how many of them have the group: precision = true_hits / hits.
    std::size_t hits = 0;
    std::size_t true_hits = 0;
    /// Even-ID molecules with the group, of all usable even-ID molecules: the base rate.
    std::size_t with_group = 0;
    std::size_t total = 0;
    /// Whether the chemist enabled the window for round 1 (2026-11-22: methoxy O-CH3 and
    /// benzene ring CH only).
    bool round1 = false;
};

/// Measures every group's 1H windows over the even-ID records of an extract read by
/// parse_shift_extract. Rows in kAllGroups order, then hallmark name, then H count.
[[nodiscard]] std::vector<ProtonShifts> measure_group_proton_shifts(
    const std::vector<ShiftRecord>& records);

/// The comment lines heading data/group_proton_shifts.tsv.
[[nodiscard]] std::string_view group_proton_shifts_header();

/// Data lines: group, hallmark, H count (CH/CH2/CH3), n, molecules, then min, 1st, 5th,
/// median, 95th, 99th, max, window low and high (ppm, two decimals), hits, true hits,
/// precision (three decimals), molecules with the group, molecules, base rate (three
/// decimals), round ("1" or "-").
[[nodiscard]] std::string to_tsv(const std::vector<ProtonShifts>& measurements);

/// `value` with exactly `decimals` decimals (rounded half away from zero), e.g. 3.81.
[[nodiscard]] std::string fixed(double value, int decimals);

}  // namespace arise::tools

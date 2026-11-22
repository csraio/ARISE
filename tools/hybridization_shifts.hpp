#pragma once

#include "arise/chem/perception.hpp"
#include "group_shifts.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

// Development-only (never part of the engine): measures the 13C shifts of sp, sp2 and sp3
// carbons over nmrshiftdb2's assigned spectra (tests/data/reference/nmrshiftdb2_13c_all.tsv),
// using ARISE's own hybridization (chem::carbon_hybridization).
namespace arise::tools {

/// Shift statistics of one carbon hybridization, in ppm.
struct HybridizationShifts {
    chem::Hybridization hybridization = chem::Hybridization::Sp3;
    std::size_t n = 0;          // carbons (each assigned carbon of each record once)
    std::size_t molecules = 0;  // distinct nmrshiftdb2 IDs among them
    double min = 0;
    double p0_1 = 0;  // percentiles as Python statistics.quantiles(n=1000, method="inclusive")
    double p1 = 0;
    double median = 0;
    double p99 = 0;
    double p99_9 = 0;
    double max = 0;
};

/// sp, sp2 and sp3 (those with any assigned carbon) measured over the records with even IDs
/// (the odd half is the hold-out for checking suggestions). Carbons ARISE leaves unassigned
/// (charged carbons) and records ARISE cannot read are left out.
[[nodiscard]] std::vector<HybridizationShifts> measure_hybridization_shifts(
    const std::vector<ShiftRecord>& records);

/// The provenance header of hybridization_shifts.tsv (comment lines).
[[nodiscard]] std::string_view hybridization_shifts_header();

/// TSV rows: hybridization, n, molecules, min, p0.1, p1, median, p99, p99.9, max (one decimal).
[[nodiscard]] std::string to_tsv(const std::vector<HybridizationShifts>& measurements);

}  // namespace arise::tools

#pragma once

#include "arise/chem/perception.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

// Development-only (never part of the engine): measures each group's hallmark 13C shifts
// over nmrshiftdb2's assigned spectra (tests/data/reference/nmrshiftdb2_13c_all.tsv),
// using ARISE's own key atoms (chem::find_group).
namespace arise::tools {

/// One record of the extract: per atom, its assigned 13C shift (empty for non-carbons and
/// carbons the database does not list).
struct ShiftRecord {
    long long id = 0;
    std::string_view smiles;
    std::vector<std::optional<double>> shifts;
};

/// The data lines of the extract (columns: id, SMILES, per-atom shift, '.' or '-' if none).
[[nodiscard]] std::vector<ShiftRecord> parse_shift_extract(std::string_view text);

/// Shift statistics of one hallmark of one group, in ppm. The O-C and N-C hallmarks
/// (ester, ether, amine) are measured separately for sp3 and for unsaturated (aromatic,
/// vinyl) carbons, as "O-C carbon (sp3)" and "O-C carbon (unsaturated)" (user decision
/// 2026-11-01): one range mixing both spans about 50-160 ppm.
struct HallmarkShifts {
    chem::Group group = chem::Group::CarboxylicAcid;
    std::string hallmark;
    std::size_t n = 0;  // carbons (each key atom of each record once)
    double min = 0;
    double p1 = 0;  // percentiles as Python statistics.quantiles(n=100, method="inclusive")
    double p5 = 0;
    double median = 0;
    double p95 = 0;
    double p99 = 0;
    double max = 0;
};

/// Every group's hallmarks measured over the records with even IDs (the odd half is the
/// hold-out for checking suggestions). Groups without carbon hallmarks are left out, and so
/// are hallmarks with no assigned carbon.
[[nodiscard]] std::vector<HallmarkShifts> measure_group_shifts(
    const std::vector<ShiftRecord>& records);

/// The provenance header of group_shifts.tsv (comment lines).
[[nodiscard]] std::string_view group_shifts_header();

/// TSV rows: group, hallmark, n, min, p1, p5, median, p95, p99, max (one decimal).
[[nodiscard]] std::string to_tsv(const std::vector<HallmarkShifts>& measurements);

}  // namespace arise::tools

#include "group_proton_shifts.hpp"

#include "arise/chem/molecule.hpp"
#include "arise/filter/groups.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <tuple>
#include <utility>

namespace arise::tools {

namespace {

// The windows the chemist enabled for round 1 (2026-11-22): the only two whose measured
// precision was above one half in the preview.
struct Enabled {
    chem::Group group;
    std::string_view hallmark;
    nmr1d::CarbonType type;
};
constexpr std::array kRound1{
    Enabled{
        .group = chem::Group::Methoxy,
        .hallmark = "O-CH3 carbon",
        .type = nmr1d::CarbonType::CH3,
    },
    Enabled{
        .group = chem::Group::BenzeneRing,
        .hallmark = "ring carbons",
        .type = nmr1d::CarbonType::CH,
    },
};

bool in_round1(chem::Group group, std::string_view hallmark, nmr1d::CarbonType type) {
    return std::any_of(kRound1.begin(), kRound1.end(), [&](const Enabled& e) {
        return e.group == group && e.hallmark == hallmark && e.type == type;
    });
}

// (group index, hallmark, H count)
using WindowKey = std::tuple<std::size_t, std::string, nmr1d::CarbonType>;

struct Values {
    std::vector<double> shifts;
    std::set<long long> molecules;
};

// One usable even-ID record: its carbons' 1H shifts with their H count, and its groups.
struct Molecule {
    std::vector<std::pair<nmr1d::CarbonType, double>> protons;
    std::vector<bool> has_group;  // per kAllGroups index
};

double percentile(const std::vector<double>& sorted, std::size_t i) {
    constexpr std::size_t kParts = 100;
    return quantile(sorted, i, kParts);
}

// The carbon's H count as a CH/CH2/CH3 type; nothing for a carbon without H.
std::optional<nmr1d::CarbonType> protonated_type(const chem::Molecule& molecule, std::size_t atom) {
    const auto type = chem::carbon_type(molecule, atom);
    if (!type.has_value() || *type == nmr1d::CarbonType::C) {
        return std::nullopt;
    }
    return type;
}

// Adds one even-ID record to `values` (its key carbons' 1H shifts) and returns its protons
// and groups, for the precision counts.
Molecule collect(const ShiftRecord& record, const chem::Molecule& molecule,
                 std::map<WindowKey, Values>& values) {
    Molecule usable;
    for (std::size_t atom = 0; atom < record.shifts.size(); ++atom) {
        const auto type = protonated_type(molecule, atom);
        const auto& shift = record.shifts[atom];
        if (type.has_value() && shift.has_value()) {
            usable.protons.emplace_back(*type, *shift);
        }
    }
    for (const auto group : chem::kAllGroups) {
        usable.has_group.push_back(chem::count_group(molecule, group) > 0);
    }
    for (const auto& key : key_carbons(molecule)) {
        const auto type = protonated_type(molecule, key.atom);
        const auto& shift = record.shifts[key.atom];
        if (type.has_value() && shift.has_value()) {
            auto& entry = values[{key.group_index, key.hallmark, *type}];
            entry.shifts.push_back(*shift);
            entry.molecules.insert(record.id);
        }
    }
    return usable;
}

// The statistics and window of one (group, hallmark, H count); `shifts` sorted.
ProtonShifts window_of(const WindowKey& key, const Values& entry) {
    const auto& shifts = entry.shifts;
    ProtonShifts row{
        .group = chem::kAllGroups[std::get<0>(key)],
        .hallmark = std::get<1>(key),
        .type = std::get<2>(key),
        .n = shifts.size(),
        .molecules = entry.molecules.size(),
        .min = shifts.front(),
        .p1 = percentile(shifts, 1),
        .p5 = percentile(shifts, 5),
        .median = median(shifts),
        .p95 = percentile(shifts, 95),
        .p99 = percentile(shifts, 99),
        .max = shifts.back(),
    };
    const bool percentiles = row.molecules >= filter::kMinMoleculesForPercentiles;
    row.low = percentiles ? row.p1 : row.min;
    row.high = percentiles ? row.p99 : row.max;
    row.round1 = in_round1(row.group, row.hallmark, row.type);
    return row;
}

// Fills `row`'s hits, true hits, molecules with the group and total over `molecules`.
void count_precision(ProtonShifts& row, std::size_t group_index,
                     const std::vector<Molecule>& molecules) {
    for (const auto& molecule : molecules) {
        const bool hit = std::any_of(
            molecule.protons.begin(), molecule.protons.end(), [&row](const auto& proton) {
                return proton.first == row.type && row.low <= proton.second &&
                       proton.second <= row.high;
            });
        const bool has = molecule.has_group[group_index];
        row.hits += hit ? 1U : 0U;
        row.true_hits += hit && has ? 1U : 0U;
        row.with_group += has ? 1U : 0U;
    }
    row.total = molecules.size();
}

}  // namespace

std::vector<ProtonShifts> measure_group_proton_shifts(const std::vector<ShiftRecord>& records) {
    std::map<WindowKey, Values> values;
    std::vector<Molecule> molecules;
    for (const auto& record : records) {
        if (record.id % 2 != 0) {
            continue;  // hold-out half
        }
        const auto parsed = chem::parse_smiles(record.smiles);
        if (parsed.ok() && parsed.value().atoms().size() == record.shifts.size()) {
            molecules.push_back(collect(record, parsed.value(), values));
        }
    }
    std::vector<ProtonShifts> result;
    for (auto& [key, entry] : values) {
        std::sort(entry.shifts.begin(), entry.shifts.end());
        auto row = window_of(key, entry);
        count_precision(row, std::get<0>(key), molecules);
        result.push_back(std::move(row));
    }
    return result;
}

std::string_view group_proton_shifts_header() {
    return "# Generated by tools/measure_group_proton_shifts.cpp from nmrshiftdb2_1h_all.tsv; do\n"
           "# not edit. Contains information from nmrshiftdb2 (www.nmrshiftdb.org), made\n"
           "# available under the nmrshiftdb2 Database License (NMRSHIFTDB2_LICENSE.txt).\n"
           "# Records with even IDs; each group's key carbons from ARISE's chem::find_group, and\n"
           "# the mean 1H shift of the H on each, split by the carbon's H count (CH/CH2/CH3).\n"
           "# Columns: group, hallmark, H count, n carbons, molecules, min, 1st, 5th percentile,\n"
           "# median, 95th, 99th percentile, max, window low and high (ppm; percentiles as\n"
           "# Python statistics.quantiles(n=100), inclusive; window 1st-99th, or min-max from\n"
           "# fewer than 100 molecules), hits (molecules with an H of that count on any carbon\n"
           "# in the window), true hits (of those, molecules with the group), precision, "
           "molecules\n"
           "# with the group, molecules, base rate, round (1 = suggested from 1H in round 1, the\n"
           "# chemist's decision 2026-11-22; '-' = measured only).\n";
}

std::string fixed(double value, int decimals) {
    const double scale = std::pow(10.0, decimals);
    const long long units = std::llround(value * scale);
    const long long magnitude = units < 0 ? -units : units;
    const auto whole = static_cast<long long>(scale);
    std::string fraction = std::to_string(magnitude % whole);
    fraction.insert(0, static_cast<std::size_t>(decimals) - fraction.size(), '0');
    return (units < 0 ? "-" : "") + std::to_string(magnitude / whole) + "." + fraction;
}

std::string to_tsv(const std::vector<ProtonShifts>& measurements) {
    constexpr int kShiftDecimals = 2;
    constexpr int kRateDecimals = 3;
    const auto rate = [](std::size_t part, std::size_t whole) {
        return whole == 0 ? 0.0 : static_cast<double>(part) / static_cast<double>(whole);
    };
    std::string text;
    for (const auto& m : measurements) {
        text.append(chem::to_string(m.group)).append("\t").append(m.hallmark).append("\t");
        text.append(nmr1d::to_string(m.type)).append("\t");
        text.append(std::to_string(m.n)).append("\t").append(std::to_string(m.molecules));
        for (const double value :
             {m.min, m.p1, m.p5, m.median, m.p95, m.p99, m.max, m.low, m.high}) {
            text.append("\t").append(fixed(value, kShiftDecimals));
        }
        text.append("\t").append(std::to_string(m.hits));
        text.append("\t").append(std::to_string(m.true_hits));
        text.append("\t").append(fixed(rate(m.true_hits, m.hits), kRateDecimals));
        text.append("\t").append(std::to_string(m.with_group));
        text.append("\t").append(std::to_string(m.total));
        text.append("\t").append(fixed(rate(m.with_group, m.total), kRateDecimals));
        text.append("\t").append(m.round1 ? "1" : "-").append("\n");
    }
    return text;
}

}  // namespace arise::tools

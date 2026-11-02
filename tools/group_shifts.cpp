#include "group_shifts.hpp"

#include "arise/chem/molecule.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>

namespace arise::tools {

namespace {

std::vector<std::string_view> split(std::string_view text, char separator) {
    std::vector<std::string_view> parts;
    while (true) {
        const auto end = text.find(separator);
        parts.push_back(text.substr(0, end));
        if (end == std::string_view::npos) {
            return parts;
        }
        text = text.substr(end + 1);
    }
}

std::optional<double> parse_shift(std::string_view text) {
    if (text == "." || text == "-" || text.empty()) {
        return std::nullopt;
    }
    return std::stod(std::string{text});
}

// Python's statistics.quantiles(data, n=100, method="inclusive"), percentile i (1..99).
double percentile(const std::vector<double>& sorted, std::size_t i) {
    constexpr std::size_t kParts = 100;
    const std::size_t m = sorted.size() - 1;
    const std::size_t j = i * m / kParts;
    const std::size_t delta = (i * m) - (j * kParts);
    if (j + 1 >= sorted.size()) {
        return sorted[j];
    }
    return ((sorted[j] * static_cast<double>(kParts - delta)) +
            (sorted[j + 1] * static_cast<double>(delta))) /
           static_cast<double>(kParts);
}

// Whether a hallmark is measured separately for sp3 and unsaturated carbons.
bool split_by_carbon_type(std::string_view hallmark) {
    return hallmark == "O-C carbon" || hallmark == "N-C carbon";
}

bool is_sp3(const chem::Molecule& molecule, std::size_t atom) {
    const auto bonds = molecule.bonds_of(atom);
    return std::all_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        return molecule.bonds()[bond].order == chem::BondOrder::Single;
    });
}

double median(const std::vector<double>& sorted) {
    const std::size_t n = sorted.size();
    return n % 2 == 1 ? sorted[n / 2] : (sorted[(n / 2) - 1] + sorted[n / 2]) / 2.0;
}

// One decimal, rounded half away from zero (no printf, so no locale or varargs).
std::string fixed1(double value) {
    const long long tenths = std::llround(value * 10.0);
    const long long magnitude = tenths < 0 ? -tenths : tenths;
    return (tenths < 0 ? "-" : "") + std::to_string(magnitude / 10) + "." +
           std::to_string(magnitude % 10);
}

// The shifts of one hallmark and the DEPT types of the carbons carrying them.
struct Values {
    std::vector<double> shifts;
    nmr1d::CarbonTypeSet types;
};

// Keyed by (group index, hallmark name with its carbon type where split).
using ValuesByHallmark = std::map<std::pair<std::size_t, std::string>, Values>;

// Adds one record's key-atom shifts to `values`.
void collect(const ShiftRecord& record, const chem::Molecule& molecule, ValuesByHallmark& values) {
    for (const auto& key : key_carbons(molecule)) {
        const auto& shift = record.shifts[key.atom];
        const auto type = chem::carbon_type(molecule, key.atom);
        if (shift && type) {
            auto& entry = values[{key.group_index, key.hallmark}];
            entry.shifts.push_back(*shift);
            entry.types = entry.types | nmr1d::CarbonTypeSet{*type};
        }
    }
}

}  // namespace

std::vector<KeyCarbon> key_carbons(const chem::Molecule& molecule) {
    std::vector<KeyCarbon> keys;
    for (std::size_t g = 0; g < chem::kAllGroups.size(); ++g) {
        const auto names = chem::hallmarks(chem::kAllGroups[g]);
        std::vector<std::set<std::size_t>> atoms(names.size());
        for (const auto& occurrence : chem::find_group(molecule, chem::kAllGroups[g])) {
            for (std::size_t h = 0; h < names.size(); ++h) {
                atoms[h].insert(occurrence.key_atoms[h].begin(), occurrence.key_atoms[h].end());
            }
        }
        for (std::size_t h = 0; h < names.size(); ++h) {
            for (const auto atom : atoms[h]) {
                std::string name{names[h]};
                if (split_by_carbon_type(name)) {
                    name += is_sp3(molecule, atom) ? " (sp3)" : " (unsaturated)";
                }
                keys.push_back({.group_index = g, .hallmark = std::move(name), .atom = atom});
            }
        }
    }
    return keys;
}

std::vector<ShiftRecord> parse_shift_extract(std::string_view text) {
    std::vector<ShiftRecord> records;
    for (const auto line : split(text, '\n')) {
        if (line.empty() || line.front() == '#') {
            continue;
        }
        const auto fields = split(line, '\t');
        if (fields.size() != 3) {
            continue;
        }
        ShiftRecord record{
            .id = std::stoll(std::string{fields[0]}),
            .smiles = fields[1],
            .shifts = {},
        };
        for (const auto value : split(fields[2], ' ')) {
            record.shifts.push_back(parse_shift(value));
        }
        records.push_back(std::move(record));
    }
    return records;
}

std::vector<HallmarkShifts> measure_group_shifts(const std::vector<ShiftRecord>& records) {
    ValuesByHallmark values;
    for (const auto& record : records) {
        if (record.id % 2 != 0) {
            continue;  // hold-out half
        }
        const auto molecule = chem::parse_smiles(record.smiles);
        if (molecule.ok() && molecule.value().atoms().size() == record.shifts.size()) {
            collect(record, molecule.value(), values);
        }
    }
    std::vector<HallmarkShifts> result;
    for (auto& [key, entry] : values) {
        auto& shifts = entry.shifts;
        std::sort(shifts.begin(), shifts.end());
        const auto group = chem::kAllGroups[key.first];
        result.push_back(HallmarkShifts{
            .group = group,
            .hallmark = key.second,
            .types = entry.types,
            .n = shifts.size(),
            .min = shifts.front(),
            .p1 = percentile(shifts, 1),
            .p5 = percentile(shifts, 5),
            .median = median(shifts),
            .p95 = percentile(shifts, 95),
            .p99 = percentile(shifts, 99),
            .max = shifts.back(),
        });
    }
    return result;
}

std::string_view group_shifts_header() {
    return "# Generated by tools/measure_group_shifts.cpp from nmrshiftdb2_13c_all.tsv; do not\n"
           "# edit. Contains information from nmrshiftdb2 (www.nmrshiftdb.org), made available\n"
           "# under the nmrshiftdb2 Database License (NMRSHIFTDB2_LICENSE.txt).\n"
           "# Records with even IDs; each group's key atoms from ARISE's chem::find_group.\n"
           "# O-C and N-C hallmarks are split into sp3 and unsaturated carbons.\n"
           "# Columns: group, hallmark, DEPT types of its carbons, n carbons, min, 1st, 5th\n"
           "# percentile, median, 95th, 99th percentile, max (ppm; percentiles as Python\n"
           "# statistics.quantiles(n=100), inclusive). Suggestion windows use 1st-99th (user\n"
           "# decision 2026-11-01).\n";
}

std::string to_tsv(const std::vector<HallmarkShifts>& measurements) {
    std::string text;
    for (const auto& m : measurements) {
        text.append(chem::to_string(m.group)).append("\t").append(m.hallmark).append("\t");
        text.append(nmr1d::to_string(m.types)).append("\t");
        text.append(std::to_string(m.n));
        for (const double value : {m.min, m.p1, m.p5, m.median, m.p95, m.p99, m.max}) {
            text.append("\t").append(fixed1(value));
        }
        text.append("\n");
    }
    return text;
}

}  // namespace arise::tools

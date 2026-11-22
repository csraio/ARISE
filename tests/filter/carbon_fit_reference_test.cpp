// filter::check_carbon_fit and order_by_fit on nmrshiftdb2's assigned 13C spectra (plan 3c).
// Each record's peak list is its distinct assigned shifts; each peak's types are those of
// the carbons assigned to it, as a complete DEPT (DEPT-90 and DEPT-135) would give them.
//   * The true structure stays Active at the measured rate; each exception (a split
//     ARISE's rules do not expect, or a database quirk) is printed for the chemist.
//   * For records that are PubChem candidates (nmrshiftdb2_pubchem_matches.tsv), every
//     candidate of the same formula (pubchem_rdkit_oracle.tsv) is ranked, and the place of
//     the true structure is reported: what this step buys.
// Both read the full extract, so they are measured natively and skipped under wasm.

#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "arise/filter/carbon_fit.hpp"
#include "group_shifts.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace {

#ifndef __EMSCRIPTEN__  // the helpers below serve tests that do not run under wasm

// The share of fully assigned nmrshiftdb2 records whose true structure check_carbon_fit
// keeps Active, measured by the test below after the chemist's nitrosamine, ring gem-pair
// and allene rules (2026-11-09) and their extensions to charge-separated nitrosamines,
// chiral allenes and sulfilimine S (2026-11-11): 27,814 of 27,901 (0.996882), recorded
// rounded down. The
// pass mark is this rate less three binomial standard errors (records are the units).
// Re-measure on purpose when the rules or the extract change.
constexpr double kMeasuredActiveRate = 0.996;

using arise::filter::CarbonFit;
using arise::filter::Tier;
using arise::nmr1d::CarbonAssignment;
using arise::nmr1d::CarbonTypeSet;

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

// The record's peak list, or none when a carbon has no assigned shift.
std::optional<std::vector<CarbonAssignment>> observed_peaks(
    const arise::chem::Molecule& molecule, const arise::tools::ShiftRecord& record) {
    if (record.shifts.size() != molecule.atoms().size()) {
        return std::nullopt;
    }
    std::map<double, CarbonTypeSet> by_shift;
    for (std::size_t atom = 0; atom < record.shifts.size(); ++atom) {
        if (molecule.atoms()[atom].element != arise::chem::Element::C) {
            continue;
        }
        const auto shift = record.shifts[atom];
        if (!shift) {
            return std::nullopt;
        }
        const auto type = arise::chem::carbon_type(molecule, atom);
        auto& possible = by_shift[*shift];
        possible = possible | (type ? CarbonTypeSet{*type} : CarbonTypeSet::all());
    }
    std::vector<CarbonAssignment> peaks;
    peaks.reserve(by_shift.size());
    for (const auto& [shift, possible] : by_shift) {
        peaks.push_back({.shift = arise::core::Ppm{shift}, .possible = possible});
    }
    return peaks;
}

const std::vector<arise::tools::ShiftRecord>& extract_records() {
    static const std::string extract = arise::test::join_chunks(arise::test::kNmrshiftdbAllChunks);
    static const auto records = arise::tools::parse_shift_extract(extract);
    return records;
}

// Element and hydrogen counts: candidates of one PubChem formula share them.
std::map<std::string, int> composition(const arise::chem::Molecule& molecule) {
    std::map<std::string, int> counts;
    for (std::size_t atom = 0; atom < molecule.atoms().size(); ++atom) {
        ++counts[std::string{arise::chem::symbol(molecule.atoms()[atom].element)}];
        counts["H"] += molecule.hydrogen_count(atom);
    }
    return counts;
}

struct Candidate {
    std::string smiles;
    arise::chem::Molecule molecule;
    std::map<std::string, int> composition;
};

std::vector<Candidate> oracle_candidates() {
    std::vector<Candidate> candidates;
    const std::string text = arise::test::join_chunks(arise::test::kRdkitOracleChunks);
    for (const auto line : split(text, '\n')) {
        if (line.empty() || line.front() == '#') {
            continue;
        }
        const auto smiles = line.substr(0, line.find('\t'));
        auto parsed = arise::chem::parse_smiles(smiles);
        if (parsed) {
            auto counts = composition(parsed.value());
            candidates.push_back({
                .smiles = std::string{smiles},
                .molecule = std::move(parsed).value(),
                .composition = std::move(counts),
            });
        }
    }
    return candidates;
}

// Where the true candidate lands: how many candidates rank strictly ahead of it, and how
// many share its ranking key (tier, overlaps, reasons) but come earlier in input order.
struct Place {
    bool active = false;
    std::size_t ahead = 0;
    std::size_t level = 0;
    std::size_t candidates = 0;
};

Place place_of(const std::vector<CarbonAssignment>& peaks, const std::vector<Candidate>& all,
               const Candidate& truth) {
    std::vector<CarbonFit> fits;
    std::optional<std::size_t> truth_index;
    for (const auto& candidate : all) {
        if (candidate.composition != truth.composition) {
            continue;
        }
        if (candidate.smiles == truth.smiles) {
            truth_index = fits.size();
        }
        const auto fit = arise::filter::check_carbon_fit(candidate.molecule, peaks);
        fits.push_back(fit ? fit.value() : CarbonFit{});
    }
    Place place{.candidates = fits.size()};
    if (!truth_index) {
        return place;
    }
    const auto& mine = fits[*truth_index];
    place.active = mine.verdict.tier == Tier::Active;
    const auto key = [](const CarbonFit& fit) {
        return std::tuple{static_cast<int>(fit.verdict.tier), fit.assumed_overlaps,
                          fit.verdict.reasons.size()};
    };
    for (std::size_t i = 0; i < fits.size(); ++i) {
        if (i == *truth_index) {
            continue;
        }
        if (key(fits[i]) < key(mine)) {
            ++place.ahead;
        } else if (key(fits[i]) == key(mine)) {
            ++place.level;
        }
    }
    return place;
}

using RecordsById = std::map<long long, const arise::tools::ShiftRecord*>;

// Ranks one nmrshiftdb2_pubchem_matches.tsv row (record id, candidate SMILES) among its
// PubChem candidates, printing its place; the rank of the true structure (worst case among
// ties), or none when the record is out of scope.
std::optional<std::size_t> rank_match(std::string_view line, const RecordsById& by_id,
                                      const std::vector<Candidate>& candidates) {
    const auto fields = split(line, '\t');
    EXPECT_EQ(fields.size(), 2U) << line;
    if (fields.size() != 2) {
        return std::nullopt;
    }
    const auto found = by_id.find(std::stoll(std::string{fields[0]}));
    if (found == by_id.end()) {
        return std::nullopt;  // outside the extract's scope (e.g. too many heavy atoms)
    }
    const auto molecule = arise::chem::parse_smiles(found->second->smiles);
    EXPECT_TRUE(molecule) << fields[0];
    if (!molecule) {
        return std::nullopt;
    }
    const auto peaks = observed_peaks(molecule.value(), *found->second);
    if (!peaks) {
        return std::nullopt;  // a carbon without an assigned shift
    }
    const auto truth = std::find_if(candidates.begin(), candidates.end(),
                                    [&](const Candidate& c) { return c.smiles == fields[1]; });
    EXPECT_NE(truth, candidates.end()) << fields[1];
    if (truth == candidates.end()) {
        return std::nullopt;
    }
    const auto place = place_of(*peaks, candidates, *truth);
    EXPECT_TRUE(place.active) << fields[0] << " " << fields[1];
    std::cout << "  " << fields[0] << " " << fields[1] << ": rank " << place.ahead + 1 << " to "
              << place.ahead + place.level + 1 << " of " << place.candidates << "\n";
    return place.ahead + place.level + 1;
}

#endif

TEST(CarbonFitReference, TrueStructuresStayActiveAtTheMeasuredRate) {
#ifdef __EMSCRIPTEN__
    GTEST_SKIP() << "measured natively over the full extract; far slower under wasm checks";
#else
    std::size_t checked = 0;
    std::size_t active = 0;
    std::map<int, std::size_t> overlaps;
    for (const auto& record : extract_records()) {
        const auto molecule = arise::chem::parse_smiles(record.smiles);
        if (!molecule) {
            continue;
        }
        const auto peaks = observed_peaks(molecule.value(), record);
        if (!peaks || peaks->empty()) {
            continue;
        }
        const auto fit = arise::filter::check_carbon_fit(molecule.value(), *peaks);
        ASSERT_TRUE(fit) << record.id;
        ++checked;
        ++overlaps[fit.value().assumed_overlaps];
        if (fit.value().verdict.tier == Tier::Active) {
            ++active;
        } else {
            std::cout << "  back burner: " << record.id << " " << record.smiles << ": "
                      << fit.value().verdict.reasons.front() << "\n";
        }
    }
    ASSERT_GT(checked, 20000U);
    const auto n = static_cast<double>(checked);
    const double rate = static_cast<double>(active) / n;
    const double pass_mark =
        kMeasuredActiveRate -
        (3.0 * std::sqrt(kMeasuredActiveRate * (1.0 - kMeasuredActiveRate) / n));
    std::cout << "[ carbon fit ] true structures checked " << checked << ", active " << active
              << " (" << rate << "), pass mark " << pass_mark << "; assumed overlaps:";
    for (const auto& [count, records] : overlaps) {
        std::cout << " " << count << ":" << records;
    }
    std::cout << "\n";
    EXPECT_GE(rate, pass_mark);
#endif
}

TEST(CarbonFitReference, TrueStructuresStayActiveAmongTheirPubChemCandidates) {
#ifdef __EMSCRIPTEN__
    GTEST_SKIP() << "measured natively over the full extract; far slower under wasm checks";
#else
    RecordsById by_id;
    for (const auto& record : extract_records()) {
        by_id.emplace(record.id, &record);
    }
    const auto candidates = oracle_candidates();
    const std::string text = arise::test::join_chunks(arise::test::kNmrshiftdbPubchemMatchesChunks);
    std::vector<std::size_t> ranks;
    for (const auto line : split(text, '\n')) {
        if (line.empty() || line.front() == '#') {
            continue;
        }
        if (const auto rank = rank_match(line, by_id, candidates)) {
            ranks.push_back(*rank);
        }
    }
    ASSERT_GT(ranks.size(), 50U);
    std::sort(ranks.begin(), ranks.end());
    std::cout << "[ carbon fit ] " << ranks.size() << " records; worst-case rank of the true "
              << "structure: median " << ranks[ranks.size() / 2] << ", 90th percentile "
              << ranks[(ranks.size() * 9) / 10] << ", worst " << ranks.back() << "\n";
#endif
}

}  // namespace

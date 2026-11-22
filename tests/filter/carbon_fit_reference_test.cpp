// filter::check_carbon_fit and order_by_fit on nmrshiftdb2's assigned 13C spectra (plan 3c).
// Each record's peak list is its distinct assigned shifts; each peak's types are those of
// the carbons assigned to it, as a complete DEPT (DEPT-90 and DEPT-135) would give them.
//   * The true structure stays Active at the measured rate; each exception (a split
//     ARISE's rules do not expect, or a database quirk) is printed for the chemist.
//   * Without DEPT (plan 3c-b): on the hold-out half (odd IDs), the chemist accepts ARISE's
//     suggested sp2/sp3 minimums and types the true sp peak count; the true structure must
//     stay Active at the measured rate. The PubChem ranking without DEPT, with and without
//     the accepted counts, is an on-demand report.
//   * For records that are PubChem candidates (nmrshiftdb2_pubchem_matches.tsv), every
//     candidate of the same formula (pubchem_rdkit_oracle.tsv) is ranked, and the place of
//     the true structure is reported: what this step buys.
// Both read the full extract, so they are measured natively and skipped under wasm.

#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "arise/filter/carbon_fit.hpp"
#include "arise/filter/hybridization.hpp"
#include "group_shifts.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

// The share of hold-out records (odd IDs) whose true structure, Active without confirmed
// counts, stays Active once the chemist accepts suggest_hybridization's minimums and types
// the true sp peak count: measured by the test below, 2026-11-22: 13,789 of 13,827
// (0.997252), recorded rounded down. Most that drop are CF3/CF2 carbons above the sp3
// range. Re-measure on purpose when the windows, the rules or the extract change.
constexpr double kMeasuredAcceptedSuggestionRate = 0.997;

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

// The peaks without DEPT: every type possible.
std::vector<CarbonAssignment> without_dept(std::vector<CarbonAssignment> peaks) {
    for (auto& peak : peaks) {
        peak.possible = CarbonTypeSet::all();
    }
    return peaks;
}

// What a chemist would confirm without DEPT: the suggested sp2/sp3 minimums, and the sp
// peak count typed from the true structure (the distinct shifts of its sp carbons).
arise::filter::ConfirmedCarbonCounts accepted_counts(const arise::chem::Molecule& molecule,
                                                     const arise::tools::ShiftRecord& record,
                                                     const std::vector<CarbonAssignment>& peaks,
                                                     const arise::filter::HybridizationWindows& w) {
    std::vector<double> sp_shifts;
    for (std::size_t atom = 0; atom < record.shifts.size(); ++atom) {
        const auto shift = record.shifts[atom];
        if (!shift ||
            arise::chem::carbon_hybridization(molecule, atom) != arise::chem::Hybridization::Sp) {
            continue;
        }
        if (std::find(sp_shifts.begin(), sp_shifts.end(), *shift) == sp_shifts.end()) {
            sp_shifts.push_back(*shift);
        }
    }
    const int sp = static_cast<int>(sp_shifts.size());
    const auto suggestion = arise::filter::suggest_hybridization(peaks, w, sp);
    EXPECT_TRUE(suggestion) << record.id;
    if (!suggestion) {
        return {};
    }
    arise::filter::ConfirmedCarbonCounts confirmed;
    if (sp > 0) {
        confirmed.sp = sp;
    }
    confirmed.sp2 = suggestion.value().sp2.min;
    confirmed.sp3 = suggestion.value().sp3.min;
    return confirmed;
}

const arise::filter::HybridizationWindows& windows() {
    static const auto measured = arise::filter::measured_hybridization_windows();
    EXPECT_TRUE(measured);
    static const auto value = measured ? measured.value() : arise::filter::HybridizationWindows{};
    return value;
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
               const Candidate& truth,
               const std::optional<arise::filter::ConfirmedCarbonCounts>& confirmed = {}) {
    std::vector<CarbonFit> fits;
    std::optional<std::size_t> truth_index;
    for (const auto& candidate : all) {
        if (candidate.composition != truth.composition) {
            continue;
        }
        if (candidate.smiles == truth.smiles) {
            truth_index = fits.size();
        }
        const auto fit =
            confirmed ? arise::filter::check_carbon_fit(candidate.molecule, peaks, *confirmed)
                      : arise::filter::check_carbon_fit(candidate.molecule, peaks);
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
// How the chemist's 13C evidence is given when ranking.
enum class Evidence : std::uint8_t {
    Dept,              // complete DEPT types
    NoDept,            // peaks only
    NoDeptWithCounts,  // peaks plus accepted hybridization counts (accepted_counts)
};

// The worst-case rank of `truth` from peaks without DEPT, with the accepted counts when
// `evidence` asks for them; prints a true structure that leaves Active.
std::size_t rank_without_dept(const arise::chem::Molecule& molecule,
                              const arise::tools::ShiftRecord& record,
                              const std::vector<CarbonAssignment>& peaks,
                              const std::vector<Candidate>& candidates, const Candidate& truth,
                              Evidence evidence) {
    const auto plain = without_dept(peaks);
    const auto place = evidence == Evidence::NoDept
                           ? place_of(plain, candidates, truth)
                           : place_of(plain, candidates, truth,
                                      accepted_counts(molecule, record, plain, windows()));
    if (!place.active) {
        std::cout << "  not active without DEPT: " << record.id << " " << truth.smiles << "\n";
    }
    return place.ahead + place.level + 1;
}

std::optional<std::size_t> rank_match(std::string_view line, const RecordsById& by_id,
                                      const std::vector<Candidate>& candidates,
                                      Evidence evidence = Evidence::Dept) {
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
    if (evidence != Evidence::Dept) {
        return rank_without_dept(molecule.value(), *found->second, *peaks, candidates, *truth,
                                 evidence);
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

#ifndef __EMSCRIPTEN__
// One hold-out record without DEPT: whether it counts (readable, fully assigned, Active
// without confirmed counts), whether it stays Active with the accepted counts, whether both
// suggested ranges were single numbers, and the first reason if it left Active.
struct AcceptedOutcome {
    bool counted = false;
    bool active = false;
    bool exact = false;
    std::string reason;
};

AcceptedOutcome accepted_outcome(const arise::tools::ShiftRecord& record) {
    const auto molecule = arise::chem::parse_smiles(record.smiles);
    if (!molecule) {
        return {};
    }
    const auto observed = observed_peaks(molecule.value(), record);
    if (!observed || observed->empty()) {
        return {};
    }
    const auto peaks = without_dept(*observed);
    const auto plain = arise::filter::check_carbon_fit(molecule.value(), peaks);
    EXPECT_TRUE(plain) << record.id;
    if (!plain || plain.value().verdict.tier != Tier::Active) {
        return {};  // measured by TrueStructuresStayActiveAtTheMeasuredRate
    }
    const auto confirmed = accepted_counts(molecule.value(), record, peaks, windows());
    const auto suggestion = arise::filter::suggest_hybridization(peaks, windows(), confirmed.sp);
    const auto fit = arise::filter::check_carbon_fit(molecule.value(), peaks, confirmed);
    EXPECT_TRUE(suggestion && fit) << record.id;
    if (!suggestion || !fit) {
        return {};
    }
    const auto& s = suggestion.value();
    AcceptedOutcome outcome{
        .counted = true,
        .active = fit.value().verdict.tier == Tier::Active,
        .exact = s.sp2.min == s.sp2.max && s.sp3.min == s.sp3.max,
        .reason = {},
    };
    if (!outcome.active) {
        outcome.reason = fit.value().verdict.reasons.front();
    }
    return outcome;
}
#endif

TEST(CarbonFitReference, AcceptedHybridizationSuggestionsKeepTrueStructuresActive) {
#ifdef __EMSCRIPTEN__
    GTEST_SKIP() << "measured natively over the full extract; far slower under wasm checks";
#else
    std::size_t checked = 0;
    std::size_t active = 0;
    std::size_t exact = 0;  // both suggested ranges a single number
    for (const auto& record : extract_records()) {
        if (record.id % 2 == 0) {
            continue;  // the windows were measured on even IDs
        }
        const auto outcome = accepted_outcome(record);
        checked += outcome.counted ? 1U : 0U;
        active += outcome.active ? 1U : 0U;
        exact += outcome.exact ? 1U : 0U;
        if (outcome.counted && !outcome.active) {
            std::cout << "  back burner: " << record.id << " " << record.smiles << ": "
                      << outcome.reason << "\n";
        }
    }
    ASSERT_GT(checked, 10000U);
    const auto n = static_cast<double>(checked);
    const double rate = static_cast<double>(active) / n;
    const double pass_mark = kMeasuredAcceptedSuggestionRate -
                             (3.0 * std::sqrt(kMeasuredAcceptedSuggestionRate *
                                              (1.0 - kMeasuredAcceptedSuggestionRate) / n));
    std::cout << "[ hybridization ] hold-out true structures checked " << checked << ", active "
              << active << " (" << rate << "), pass mark " << pass_mark
              << "; both ranges exact for " << exact << "\n";
    EXPECT_GE(rate, pass_mark);
#endif
}

// A report, not a check, taking about 5 minutes, so it runs only on demand: run
// build/dev/tests/filter/arise_filter_tests with --gtest_also_run_disabled_tests and
// --gtest_filter='*PubChemRanksWithoutDept*'.
// Measured 2026-11-22 over 114 records: median worst-case rank 340 without DEPT, 120 with
// accepted sp2/sp3 counts (38 with a complete DEPT).
TEST(CarbonFitReference, DISABLED_PubChemRanksWithoutDeptWithAndWithoutAcceptedCounts) {
#ifdef __EMSCRIPTEN__
    GTEST_SKIP() << "measured natively over the full extract; far slower under wasm checks";
#else
    RecordsById by_id;
    for (const auto& record : extract_records()) {
        by_id.emplace(record.id, &record);
    }
    const auto candidates = oracle_candidates();
    const std::string text = arise::test::join_chunks(arise::test::kNmrshiftdbPubchemMatchesChunks);
    for (const auto evidence : {Evidence::NoDept, Evidence::NoDeptWithCounts}) {
        std::vector<std::size_t> ranks;
        for (const auto line : split(text, '\n')) {
            if (line.empty() || line.front() == '#') {
                continue;
            }
            if (const auto rank = rank_match(line, by_id, candidates, evidence)) {
                ranks.push_back(*rank);
            }
        }
        ASSERT_GT(ranks.size(), 50U);
        std::sort(ranks.begin(), ranks.end());
        std::cout << "[ carbon fit ] without DEPT"
                  << (evidence == Evidence::NoDept ? "" : ", accepted sp2 and sp3 counts") << ": "
                  << ranks.size() << " records; worst-case rank: median " << ranks[ranks.size() / 2]
                  << ", 90th percentile " << ranks[(ranks.size() * 9) / 10] << ", worst "
                  << ranks.back() << "\n";
    }
#endif
}

}  // namespace

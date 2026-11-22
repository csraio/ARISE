// bench/prune_pipeline (plan 3d): engineering tests of the benchmark's chain and timing
// helpers. The chain must rank exactly as calling the filter functions directly; the
// benchmark case is read from tests/data/reference/benchmark_spectrum.tsv.

#include "prune_pipeline.hpp"

#include "arise/chem/perception.hpp"
#include "pubchem_fixture.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::bench::BenchmarkCase;
using arise::bench::evidence_for;
using arise::bench::parse_benchmark_case;
using arise::bench::Scenario;
using arise::core::ErrorCode;
using arise::filter::Tier;

std::string case_text() {
    return arise::test::join_chunks(arise::test::kBenchmarkSpectrumChunks);
}

BenchmarkCase benchmark_case() {
    const auto parsed = parse_benchmark_case(case_text());
    EXPECT_TRUE(parsed) << (parsed ? "" : parsed.error().message);
    return parsed ? parsed.value() : BenchmarkCase{};
}

// The value of the first line of `text` starting with "<kind>\t".
std::string line_value(const std::string& text, std::string_view kind) {
    const auto start = text.find("\n" + std::string{kind} + "\t");
    if (start == std::string::npos) {
        return {};
    }
    const auto from = start + kind.size() + 2;
    return text.substr(from, text.find('\n', from) - from);
}

// ---- The benchmark case ------------------------------------------------------------

TEST(BenchmarkCase, ReadsTheBuildersCase) {
    const auto text = case_text();
    const auto c = benchmark_case();
    EXPECT_EQ(c.formula, line_value(text, "formula"));
    EXPECT_EQ(std::to_string(c.record), line_value(text, "record"));
    EXPECT_EQ(c.structure, line_value(text, "structure"));
    std::size_t peak_lines = 0;
    for (std::size_t at = text.find("\npeak\t"); at != std::string::npos;
         at = text.find("\npeak\t", at + 1)) {
        ++peak_lines;
    }
    EXPECT_EQ(c.peaks.size(), peak_lines);
    EXPECT_GT(c.peaks.size(), 1U);
}

TEST(BenchmarkCase, TheTrueStructureHasTheCasesFormula) {
    const auto c = benchmark_case();
    const auto molecule = arise::chem::parse_smiles(c.structure);
    ASSERT_TRUE(molecule) << c.structure;
    int carbons = 0;
    for (const auto& atom : molecule.value().atoms()) {
        carbons += atom.element == arise::chem::Element::C ? 1 : 0;
    }
    EXPECT_TRUE(c.formula.starts_with("C" + std::to_string(carbons) + "H")) << c.formula;
}

void expect_case_error(const std::string& text, std::string_view needle) {
    const auto parsed = parse_benchmark_case(text);
    ASSERT_FALSE(parsed) << text;
    EXPECT_EQ(parsed.error().code, ErrorCode::ParseError);
    EXPECT_NE(parsed.error().message.find(needle), std::string::npos) << parsed.error().message;
}

TEST(BenchmarkCase, RejectsMalformedCases) {
    const auto text = case_text();
    const auto record = "\nrecord\t" + line_value(text, "record");
    const auto without = [&](const std::string& part) {
        auto copy = text;
        copy.erase(copy.find(part), part.size());
        return copy;
    };
    expect_case_error(without(record), "record");                // missing
    expect_case_error(text + record.substr(1) + "\n", "twice");  // repeated
    expect_case_error(text + "colour\tblue\n", "unknown");
    expect_case_error(text + "peak\tabc\tCH\n", "line");
    expect_case_error(text + "peak\t12.5\tCH5\n", "line");
    expect_case_error(text + "peak\t12.5\n", "line");
    expect_case_error(text + "record\t\n", "line");
    expect_case_error("# comments only\n", "formula");
}

// ---- Evidence ----------------------------------------------------------------------

TEST(BenchmarkEvidence, WorstConfirmsNothingButTheDeptTypes) {
    const auto c = benchmark_case();
    const auto evidence = evidence_for(Scenario::Worst, c);
    ASSERT_TRUE(evidence);
    EXPECT_EQ(evidence.value().peaks, c.peaks);
    EXPECT_TRUE(evidence.value().groups.present.empty());
    EXPECT_TRUE(evidence.value().groups.absent.empty());
    EXPECT_FALSE(evidence.value().exchangeable);
    EXPECT_FALSE(evidence.value().counts);
}

// Groups whose listing in `present` disagrees with the true structure having them.
std::vector<std::string> group_mismatches(const arise::chem::Molecule& truth,
                                          const std::vector<arise::chem::Group>& present) {
    std::vector<std::string> mismatches;
    for (const auto group : arise::chem::kAllGroups) {
        const bool listed = std::find(present.begin(), present.end(), group) != present.end();
        if (listed != (arise::chem::count_group(truth, group) > 0)) {
            mismatches.emplace_back(arise::chem::to_string(group));
        }
    }
    return mismatches;
}

TEST(BenchmarkEvidence, TypicalConfirmsTheTrueStructuresGroupsAndExchangeableCount) {
    const auto c = benchmark_case();
    const auto truth = arise::chem::parse_smiles(c.structure);
    ASSERT_TRUE(truth);
    const auto evidence = evidence_for(Scenario::Typical, c);
    ASSERT_TRUE(evidence);
    EXPECT_EQ(group_mismatches(truth.value(), evidence.value().groups.present),
              std::vector<std::string>{});
    const int oh_nh = arise::chem::exchangeable_hydrogens(truth.value());
    EXPECT_EQ(evidence.value().exchangeable,
              (arise::filter::ConfirmedExchangeable{.min = oh_nh, .max = oh_nh}));
}

// `peaks` with every DEPT type possible.
std::vector<arise::nmr1d::CarbonAssignment> untyped(
    std::vector<arise::nmr1d::CarbonAssignment> peaks) {
    for (auto& peak : peaks) {
        peak.possible = arise::nmr1d::CarbonTypeSet::all();
    }
    return peaks;
}

TEST(BenchmarkEvidence, NoDeptDropsTypesAndAcceptsTheSuggestedMinimums) {
    const auto c = benchmark_case();
    const auto evidence = evidence_for(Scenario::NoDept, c);
    ASSERT_TRUE(evidence);
    EXPECT_EQ(evidence.value().peaks, untyped(c.peaks));
    const auto windows = arise::filter::measured_hybridization_windows();
    ASSERT_TRUE(windows);
    const auto counts = evidence.value().counts.value_or(arise::filter::ConfirmedCarbonCounts{});
    EXPECT_TRUE(evidence.value().counts);
    const auto suggestion = arise::filter::suggest_hybridization(
        evidence.value().peaks, windows.value(), counts.sp.value_or(0));
    ASSERT_TRUE(suggestion);
    EXPECT_EQ(counts.sp2, suggestion.value().sp2.min);
    EXPECT_EQ(counts.sp3, suggestion.value().sp3.min);
}

TEST(BenchmarkEvidence, RejectsAnUnreadableStructure) {
    auto c = benchmark_case();
    c.structure = "";
    EXPECT_FALSE(evidence_for(Scenario::Typical, c));
}

// ---- The chain -----------------------------------------------------------------------

// The C4H8O2 search's candidates, read by the chain's first two stages.
arise::bench::Parsed small_search() {
    const auto imported = arise::bench::clean(arise::test::kPubChemC4H8O2);
    EXPECT_TRUE(imported);
    if (!imported) {
        return {};
    }
    auto parsed = arise::bench::parse(imported.value());
    EXPECT_EQ(parsed.molecules.size() + parsed.unreadable, imported.value().candidates.size());
    EXPECT_FALSE(parsed.molecules.empty());
    return parsed;
}

// check_carbon_fit called on each molecule directly.
std::vector<arise::filter::CarbonFit> direct_fits(
    const std::vector<arise::chem::Molecule>& molecules,
    const std::vector<arise::nmr1d::CarbonAssignment>& peaks) {
    std::vector<arise::filter::CarbonFit> fits;
    for (const auto& molecule : molecules) {
        const auto fit = arise::filter::check_carbon_fit(molecule, peaks);
        EXPECT_TRUE(fit);
        fits.push_back(fit ? fit.value() : arise::filter::CarbonFit{});
    }
    return fits;
}

TEST(BenchmarkChain, RanksExactlyAsTheFilterFunctionsDo) {
    // A small real search (C4H8O2) with the case's DEPT peaks: an engineering check that the
    // chain adds nothing of its own.
    const auto parsed = small_search();
    const auto evidence = evidence_for(Scenario::Worst, benchmark_case());
    ASSERT_TRUE(evidence);
    const auto fits = arise::bench::prune(parsed.molecules, evidence.value());
    ASSERT_TRUE(fits);
    const auto direct = direct_fits(parsed.molecules, evidence.value().peaks);
    EXPECT_EQ(fits.value(), direct);
    EXPECT_EQ(arise::bench::rank(fits.value()), arise::filter::order_by_fit(direct));
}

// Indices whose Eliminated verdict disagrees with the molecule having OH/NH protons.
std::vector<std::size_t> elimination_mismatches(const std::vector<arise::chem::Molecule>& molecules,
                                                const std::vector<arise::filter::CarbonFit>& fits) {
    std::vector<std::size_t> mismatches;
    for (std::size_t i = 0; i < molecules.size() && i < fits.size(); ++i) {
        const bool has_oh_nh = arise::chem::exchangeable_hydrogens(molecules[i]) > 0;
        if ((fits[i].verdict.tier == Tier::Eliminated) != has_oh_nh) {
            mismatches.push_back(i);
        }
    }
    return mismatches;
}

TEST(BenchmarkChain, HardChecksEliminateBeforeTheCarbonFit) {
    const auto parsed = small_search();
    arise::bench::Evidence evidence;
    evidence.peaks = benchmark_case().peaks;
    evidence.exchangeable = arise::filter::ConfirmedExchangeable{.min = 0, .max = 0};
    const auto fits = arise::bench::prune(parsed.molecules, evidence);
    ASSERT_TRUE(fits);
    ASSERT_EQ(fits.value().size(), parsed.molecules.size());
    EXPECT_EQ(elimination_mismatches(parsed.molecules, fits.value()), std::vector<std::size_t>{});
    const auto eliminated =
        std::count_if(fits.value().begin(), fits.value().end(),
                      [](const auto& f) { return f.verdict.tier == Tier::Eliminated; });
    EXPECT_GT(eliminated, 0);
    EXPECT_LT(static_cast<std::size_t>(eliminated), parsed.molecules.size());
}

TEST(BenchmarkChain, RejectsContradictoryGroups) {
    arise::bench::Evidence evidence;
    evidence.peaks = benchmark_case().peaks;
    evidence.groups.present = {arise::chem::Group::Ketone};
    evidence.groups.absent = {arise::chem::Group::Ketone};
    const auto parsed = small_search();
    EXPECT_FALSE(arise::bench::prune(parsed.molecules, evidence));
}

#ifndef __EMSCRIPTEN__
// The true structure's tier in each scenario over the full case, as "scenario: tier".
std::vector<std::string> true_structure_tiers(const BenchmarkCase& c,
                                              const arise::bench::Parsed& parsed) {
    std::vector<std::string> tiers;
    const auto truth = std::find(parsed.smiles.begin(), parsed.smiles.end(), c.structure);
    if (truth == parsed.smiles.end()) {
        return {"the true structure is not among the candidates"};
    }
    const auto index = static_cast<std::size_t>(truth - parsed.smiles.begin());
    for (const auto scenario : {Scenario::Worst, Scenario::Typical, Scenario::NoDept}) {
        const std::string name{arise::bench::to_string(scenario)};
        const auto evidence = evidence_for(scenario, c);
        if (!evidence) {
            tiers.push_back(name + ": " + evidence.error().message);
            continue;
        }
        const auto fits = arise::bench::prune(parsed.molecules, evidence.value());
        tiers.push_back(
            name + ": " +
            (fits ? std::string{arise::filter::to_string(fits.value()[index].verdict.tier)}
                  : fits.error().message));
    }
    return tiers;
}
#endif

TEST(BenchmarkReference, TheTrueStructureStaysActiveInEveryScenario) {
    // The full case: every candidate of the 10,000+-row search, as arise_bench runs it.
#ifdef __EMSCRIPTEN__
    GTEST_SKIP() << "the full case is timed by arise_bench; far slower under wasm checks";
#else
    const auto c = benchmark_case();
    const auto csv = arise::test::join_chunks(arise::test::kBenchmarkPubchemChunks);
    const auto imported = arise::bench::clean(csv);
    ASSERT_TRUE(imported);
    EXPECT_GE(imported.value().rows, 10000U);
    const auto parsed = arise::bench::parse(imported.value());
    EXPECT_EQ(true_structure_tiers(c, parsed),
              (std::vector<std::string>{"worst: active", "typical: active", "no-dept: active"}));
#endif
}

// ---- Timing helpers --------------------------------------------------------------------

TEST(BenchmarkSummary, GivesMinMedianMax) {
    const auto odd = arise::bench::summarize({5.0, 1.0, 3.0});
    ASSERT_TRUE(odd);
    EXPECT_EQ(odd.value().min, 1.0);
    EXPECT_EQ(odd.value().median, 3.0);
    EXPECT_EQ(odd.value().max, 5.0);
    const auto even = arise::bench::summarize({4.0, 1.0, 2.0, 8.0});
    ASSERT_TRUE(even);
    EXPECT_EQ(even.value().median, 3.0);
    const auto one = arise::bench::summarize({2.5});
    ASSERT_TRUE(one);
    EXPECT_EQ(one.value().min, 2.5);
    EXPECT_EQ(one.value().median, 2.5);
    EXPECT_EQ(one.value().max, 2.5);
}

TEST(BenchmarkSummary, RejectsNoSamplesAndBadSamples) {
    EXPECT_FALSE(arise::bench::summarize({}));
    EXPECT_FALSE(arise::bench::summarize({1.0, -0.5}));
    EXPECT_FALSE(arise::bench::summarize({1.0, std::numeric_limits<double>::quiet_NaN()}));
    EXPECT_FALSE(arise::bench::summarize({std::numeric_limits<double>::infinity()}));
}

arise::core::Result<arise::bench::Options> options(std::vector<std::string_view> args) {
    return arise::bench::parse_options(args);
}

TEST(BenchmarkOptions, DefaultsAndValues) {
    const auto none = options({});
    ASSERT_TRUE(none);
    EXPECT_EQ(none.value().runs, 15U);
    EXPECT_FALSE(none.value().budget_ms);
    const auto set = options({"--runs", "3", "--budget-ms", "100"});
    ASSERT_TRUE(set);
    EXPECT_EQ(set.value().runs, 3U);
    EXPECT_EQ(set.value().budget_ms, 100.0);
    const auto decimal = options({"--budget-ms", "2.5"});
    ASSERT_TRUE(decimal);
    EXPECT_EQ(decimal.value().budget_ms, 2.5);
}

TEST(BenchmarkOptions, RejectsBadArguments) {
    for (const auto& args : std::vector<std::vector<std::string_view>>{
             {"--runs"},
             {"--runs", "0"},
             {"--runs", "-2"},
             {"--runs", "two"},
             {"--runs", "1.5"},
             {"--budget-ms", "0"},
             {"--budget-ms", "-1"},
             {"--budget-ms", "fast"},
             {"--budget-ms"},
             {"--verbose"},
         }) {
        const auto result = options(args);
        ASSERT_FALSE(result) << args.front();
        EXPECT_EQ(result.error().code, ErrorCode::InvalidArgument);
    }
}

TEST(BenchmarkExitCode, FailsOnlyOverTheBudget) {
    EXPECT_EQ(arise::bench::exit_code(250.0, std::nullopt), 0);
    EXPECT_EQ(arise::bench::exit_code(99.9, 100.0), 0);
    EXPECT_EQ(arise::bench::exit_code(100.0, 100.0), 0);
    EXPECT_EQ(arise::bench::exit_code(100.1, 100.0), 1);
}

TEST(BenchmarkScenario, NamesEachScenario) {
    EXPECT_EQ(arise::bench::to_string(Scenario::Worst), "worst");
    EXPECT_EQ(arise::bench::to_string(Scenario::Typical), "typical");
    EXPECT_EQ(arise::bench::to_string(Scenario::NoDept), "no-dept");
}

}  // namespace

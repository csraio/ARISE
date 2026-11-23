// pipeline (plan "web client"): the stateless calls the page makes. Each must give exactly
// what calling the underlying modules directly gives; the real-data case is the committed
// C15H20O3 PubChem search with the benchmark's nmrshiftdb2 spectrum
// (tests/data/reference/benchmark_spectrum.tsv).

#include "arise/pipeline/pipeline.hpp"

#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "arise/core/error.hpp"
#include "arise/core/formula.hpp"
#include "arise/filter/carbon_fit.hpp"
#include "arise/filter/groups.hpp"
#include "arise/io/carbon_csv.hpp"
#include "arise/io/pasted_peaks.hpp"
#include "arise/io/pubchem.hpp"
#include "arise/nmr1d/dept.hpp"
#include "pubchem_fixture.hpp"
#include "reference_fixtures.hpp"
#include "test_molecules.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::filter::ConfirmedGroups;
using arise::filter::SuggestOptions;
using arise::filter::Tier;
using arise::nmr1d::CarbonType;
using arise::pipeline::group_named;
using arise::pipeline::rank;
using arise::pipeline::read_spectra;
using arise::pipeline::search_url;
using arise::pipeline::SpectraText;
using arise::pipeline::SpectrumText;
using arise::pipeline::suggest;
using arise::pipeline::TextFormat;
using arise::test::compound_cid;
using arise::test::compound_formula;
using arise::test::compound_smiles;

// ---- the benchmark case as the chemist would enter it ----------------------------

struct CasePeak {
    std::string shift;  // as written in benchmark_spectrum.tsv
    CarbonType type;
};

struct Case {
    std::string formula;
    std::vector<CasePeak> peaks;
};

// The fields of a tab-separated line.
std::vector<std::string_view> fields(std::string_view line) {
    std::vector<std::string_view> parts;
    while (true) {
        const auto tab = line.find('\t');
        parts.push_back(line.substr(0, tab));
        if (tab == std::string_view::npos) {
            return parts;
        }
        line.remove_prefix(tab + 1);
    }
}

CarbonType type_named(std::string_view name) {
    for (const auto type : {CarbonType::C, CarbonType::CH, CarbonType::CH2, CarbonType::CH3}) {
        if (arise::nmr1d::to_string(type) == name) {
            return type;
        }
    }
    ADD_FAILURE() << "unknown carbon type '" << name << "'";
    return CarbonType::C;
}

// benchmark_spectrum.tsv: "formula", "record", "structure" and "peak" lines; each peak has
// one type (RDKit's hydrogen count on the record's structure, written by the builder).
Case benchmark_case() {
    const std::string text = arise::test::join_chunks(arise::test::kBenchmarkSpectrumChunks);
    Case c;
    std::string_view rest = text;
    while (!rest.empty()) {
        const auto end = rest.find('\n');
        const auto line = rest.substr(0, end);
        rest = end == std::string_view::npos ? std::string_view{} : rest.substr(end + 1);
        if (line.empty() || line.front() == '#') {
            continue;
        }
        const auto parts = fields(line);
        if (parts[0] == "formula") {
            c.formula = std::string{parts[1]};
        } else if (parts[0] == "peak") {
            EXPECT_EQ(parts.size(), 3U) << line;
            c.peaks.push_back({.shift = std::string{parts[1]}, .type = type_named(parts[2])});
        }
    }
    EXPECT_FALSE(c.formula.empty());
    EXPECT_FALSE(c.peaks.empty());
    return c;
}

// The spectra the case's types imply, pasted as shift<TAB>intensity: DEPT-135 up for CH
// and CH3, down for CH2, absent for C; DEPT-90 up for CH only.
SpectraText pasted_spectra(const Case& c) {
    std::string carbon;
    std::string dept135;
    std::string dept90;
    for (const auto& peak : c.peaks) {
        carbon += peak.shift + "\n";
        if (peak.type == CarbonType::CH || peak.type == CarbonType::CH3) {
            dept135 += peak.shift + "\t100\n";
        } else if (peak.type == CarbonType::CH2) {
            dept135 += peak.shift + "\t-100\n";
        }
        if (peak.type == CarbonType::CH) {
            dept90 += peak.shift + "\t100\n";
        }
    }
    return {
        .carbon = {.format = TextFormat::Pasted, .text = carbon},
        .dept135 = SpectrumText{.format = TextFormat::Pasted, .text = dept135},
        .dept90 = SpectrumText{.format = TextFormat::Pasted, .text = dept90},
    };
}

// The text of a DEPT spectrum, or none when it is absent.
std::string text_of(const std::optional<SpectrumText>& spectrum) {
    return spectrum.has_value() ? spectrum->text : std::string{};
}

// The same spectra as CSV text in ARISE's peak-list format.
SpectraText csv_spectra(const Case& c) {
    const auto pasted = pasted_spectra(c);
    const auto to_csv = [](std::string text, std::string_view header) {
        for (auto& ch : text) {
            ch = ch == '\t' ? ',' : ch;
        }
        return SpectrumText{.format = TextFormat::Csv, .text = std::string{header} + text};
    };
    return {
        .carbon = to_csv(pasted.carbon.text, "shift_ppm\n"),
        .dept135 = to_csv(text_of(pasted.dept135), "shift_ppm,intensity\n"),
        .dept90 = to_csv(text_of(pasted.dept90), "shift_ppm,intensity\n"),
    };
}

// ---- search_url -------------------------------------------------------------------

TEST(PipelineSearchUrl, IsTheIoAddressForTheFormula) {
    const auto formula = compound_formula("Acetone");
    const auto parsed = arise::core::Formula::parse(formula);
    ASSERT_TRUE(parsed.ok()) << parsed.error().message;
    const auto url = search_url(formula);
    ASSERT_TRUE(url.ok()) << url.error().message;
    EXPECT_EQ(url.value(), arise::io::pubchem_formula_search_url(parsed.value()));
}

TEST(PipelineSearchUrl, AFormulaThatDoesNotReadIsAParseError) {
    const auto url = search_url("not a formula");
    ASSERT_FALSE(url.ok());
    EXPECT_EQ(url.error().code, ErrorCode::ParseError);
}

// ---- read_spectra -----------------------------------------------------------------

TEST(PipelineReadSpectra, ClassifiesTheBenchmarkPeaksAsTheirTypes) {
    const auto c = benchmark_case();
    const auto peaks = read_spectra(pasted_spectra(c));
    ASSERT_TRUE(peaks.ok()) << peaks.error().message;
    ASSERT_EQ(peaks.value().size(), c.peaks.size());
    for (std::size_t i = 0; i < c.peaks.size(); ++i) {
        EXPECT_EQ(peaks.value()[i].possible.resolved(), std::optional{c.peaks[i].type})
            << "peak " << c.peaks[i].shift;
    }
}

TEST(PipelineReadSpectra, CsvAndPastedTextGiveTheSamePeaks) {
    const auto c = benchmark_case();
    const auto pasted = read_spectra(pasted_spectra(c));
    const auto csv = read_spectra(csv_spectra(c));
    ASSERT_TRUE(pasted.ok()) << pasted.error().message;
    ASSERT_TRUE(csv.ok()) << csv.error().message;
    EXPECT_EQ(pasted.value(), csv.value());
}

TEST(PipelineReadSpectra, IsClassifyCarbonsOnTheParsedLists) {
    const auto c = benchmark_case();
    const auto spectra = pasted_spectra(c);
    const auto carbon = arise::io::parse_pasted_carbon(spectra.carbon.text);
    const auto dept135 = arise::io::parse_pasted_dept(text_of(spectra.dept135));
    const auto dept90 = arise::io::parse_pasted_dept(text_of(spectra.dept90));
    ASSERT_TRUE(carbon.ok() && dept135.ok() && dept90.ok());
    const auto direct = arise::nmr1d::classify_carbons(
        {.carbon = carbon.value(), .dept135 = dept135.value(), .dept90 = dept90.value()});
    ASSERT_TRUE(direct.ok()) << direct.error().message;
    const auto peaks = read_spectra(spectra);
    ASSERT_TRUE(peaks.ok()) << peaks.error().message;
    EXPECT_EQ(peaks.value(), direct.value());
}

TEST(PipelineReadSpectra, WithoutDeptEveryTypeIsPossible) {
    auto spectra = pasted_spectra(benchmark_case());
    spectra.dept135.reset();
    spectra.dept90.reset();
    const auto peaks = read_spectra(spectra);
    ASSERT_TRUE(peaks.ok()) << peaks.error().message;
    for (const auto& peak : peaks.value()) {
        EXPECT_EQ(peak.possible, arise::nmr1d::CarbonTypeSet::all());
    }
}

TEST(PipelineReadSpectra, ErrorsNameTheSpectrumAndKeepTheirCode) {
    auto spectra = pasted_spectra(benchmark_case());
    const auto expect_error = [](const SpectraText& s, ErrorCode code, std::string_view prefix) {
        const auto result = read_spectra(s);
        ASSERT_FALSE(result.ok());
        EXPECT_EQ(result.error().code, code) << result.error().message;
        EXPECT_TRUE(result.error().message.starts_with(prefix)) << result.error().message;
    };
    auto bad_carbon = spectra;
    bad_carbon.carbon.text = "29,5\n";
    expect_error(bad_carbon, ErrorCode::ParseError, "13C spectrum: line 1");
    auto far_carbon = spectra;
    far_carbon.carbon.text += "400\n";
    expect_error(far_carbon, ErrorCode::OutOfRange, "13C spectrum: line");
    auto bad_135 = spectra;
    bad_135.dept135 =
        SpectrumText{.format = TextFormat::Pasted, .text = text_of(spectra.dept135) + "29\n"};
    expect_error(bad_135, ErrorCode::ParseError, "DEPT-135 spectrum: line");
    auto bad_90 = spectra;
    bad_90.dept90 = SpectrumText{.format = TextFormat::Csv, .text = "shift_ppm\n29\n"};
    expect_error(bad_90, ErrorCode::ParseError, "DEPT-90 spectrum: line");
}

// ---- suggest ----------------------------------------------------------------------

TEST(PipelineSuggest, IsSuggestGroupsWithTheMeasuredWindows) {
    const auto spectra = pasted_spectra(benchmark_case());
    const auto peaks = read_spectra(spectra);
    const auto windows = arise::filter::measured_hallmark_windows();
    ASSERT_TRUE(peaks.ok() && windows.ok());
    for (const bool use_dept : {true, false}) {
        const SuggestOptions options{.use_dept = use_dept};
        const auto suggestions = suggest(spectra, options);
        ASSERT_TRUE(suggestions.ok()) << suggestions.error().message;
        EXPECT_EQ(suggestions.value(),
                  arise::filter::suggest_groups(peaks.value(), windows.value(), options));
    }
}

TEST(PipelineSuggest, PassesSpectrumErrorsOn) {
    auto spectra = pasted_spectra(benchmark_case());
    spectra.carbon.text.clear();
    const auto suggestions = suggest(spectra, {});
    ASSERT_FALSE(suggestions.ok());
    EXPECT_TRUE(suggestions.error().message.starts_with("13C spectrum"))
        << suggestions.error().message;
}

// ---- text formats -----------------------------------------------------------------

TEST(PipelineTextFormat, NamesReadBack) {
    for (const auto format : {TextFormat::Csv, TextFormat::Pasted}) {
        EXPECT_EQ(arise::pipeline::text_format_named(arise::pipeline::to_string(format)),
                  std::optional{format});
    }
    EXPECT_EQ(arise::pipeline::to_string(TextFormat::Csv), "csv");
    EXPECT_EQ(arise::pipeline::to_string(TextFormat::Pasted), "pasted");
    EXPECT_EQ(arise::pipeline::text_format_named("CSV"), std::nullopt);
    EXPECT_EQ(arise::pipeline::text_format_named(""), std::nullopt);
}

// ---- group_named ------------------------------------------------------------------

TEST(PipelineGroupNamed, ReadsEveryGroupByItsName) {
    for (const auto group : arise::chem::kAllGroups) {
        EXPECT_EQ(group_named(arise::chem::to_string(group)), std::optional{group});
    }
    EXPECT_EQ(group_named("no such group"), std::nullopt);
    EXPECT_EQ(group_named(""), std::nullopt);
}

// ---- rank -------------------------------------------------------------------------

// PubChem's CSV for the given compounds, as PubChem writes it.
std::string pubchem_csv(const std::vector<std::string_view>& names) {
    std::string csv = "\"CID\",\"ConnectivitySMILES\"\n";
    for (const auto name : names) {
        csv += std::to_string(compound_cid(name)) + ",\"" + compound_smiles(name) + "\"\n";
    }
    return csv;
}

// Two compounds with one formula (PubChem's, in test_compounds.tsv).
std::vector<std::string_view> same_formula() {
    return {"Acetone", "Propylene Oxide"};
}

std::string formula_of_same_formula() {
    return compound_formula("Acetone");
}

SpectraText some_spectrum() {
    // The benchmark's 13C list, as text: the carbon fit needs peaks, and these tests
    // compare with direct filter calls instead of expecting particular verdicts.
    auto spectra = pasted_spectra(benchmark_case());
    spectra.dept135.reset();
    spectra.dept90.reset();
    return spectra;
}

constexpr SuggestOptions kDefaults{};

ConfirmedGroups no_groups() {
    return {};
}

// check_confirmed_groups on the molecule `smiles`; fails the test if that cannot run.
arise::filter::Verdict group_verdict(const std::string& smiles, const ConfirmedGroups& confirmed) {
    const auto molecule = arise::chem::parse_smiles(smiles);
    const auto verdict = molecule.ok()
                             ? arise::filter::check_confirmed_groups(molecule.value(), confirmed)
                             : arise::core::Result<arise::filter::Verdict>{molecule.error()};
    if (!verdict.ok()) {
        ADD_FAILURE() << smiles << ": " << verdict.error().message;
        return {};
    }
    return verdict.value();
}

// Expects `candidate` to have the verdict check_confirmed_groups gives when that
// eliminates it (before any carbon fit), and no elimination otherwise.
void expect_group_verdict(const arise::pipeline::RankedCandidate& candidate,
                          const ConfirmedGroups& confirmed) {
    SCOPED_TRACE("CID " + std::to_string(candidate.cid));
    const auto expected = group_verdict(candidate.smiles, confirmed);
    if (expected.tier == Tier::Eliminated) {
        EXPECT_EQ(candidate.verdict, expected);
        EXPECT_EQ(candidate.environments_most, 0);
    } else {
        EXPECT_NE(candidate.verdict.tier, Tier::Eliminated);
    }
}

TEST(PipelineRank, AConfirmedGroupEliminatesTheCandidatesWithoutIt) {
    ASSERT_EQ(compound_formula("Acetone"), compound_formula("Propylene Oxide"));
    const ConfirmedGroups ketone{.present = {arise::chem::Group::Ketone}, .absent = {}};
    const auto ranking = rank(formula_of_same_formula(), pubchem_csv(same_formula()),
                              some_spectrum(), ketone, kDefaults);
    ASSERT_TRUE(ranking.ok()) << ranking.error().message;
    for (const auto& candidate : ranking.value().candidates) {
        expect_group_verdict(candidate, ketone);
    }
    // Exactly the candidate without a ketone is eliminated, and it ranks last.
    const auto& last = ranking.value().candidates.back();
    EXPECT_EQ(last.verdict.tier, Tier::Eliminated);
    EXPECT_EQ(last.cid, static_cast<std::uint64_t>(compound_cid("Propylene Oxide")));
    EXPECT_NE(ranking.value().candidates.front().verdict.tier, Tier::Eliminated);
}

TEST(PipelineRank, AGroupConfirmedBothWaysIsInvalid) {
    const ConfirmedGroups both{
        .present = {arise::chem::Group::Ketone},
        .absent = {arise::chem::Group::Ketone},
    };
    const auto ranking = rank(formula_of_same_formula(), pubchem_csv(same_formula()),
                              some_spectrum(), both, kDefaults);
    ASSERT_FALSE(ranking.ok());
    EXPECT_EQ(ranking.error().code, ErrorCode::InvalidArgument);
}

TEST(PipelineRank, AnUnreadableSmilesIsEliminatedWithItsReason) {
    const auto broken = arise::test::syntax_case("aromatic lowercase atom after a Kekule ring");
    ASSERT_FALSE(arise::chem::parse_smiles(broken).ok());
    const std::string csv = pubchem_csv(same_formula()) + "999999999,\"" + broken + "\"\n";
    const auto ranking =
        rank(formula_of_same_formula(), csv, some_spectrum(), no_groups(), kDefaults);
    ASSERT_TRUE(ranking.ok()) << ranking.error().message;
    const auto& candidates = ranking.value().candidates;
    ASSERT_EQ(candidates.size(), 3U);
    const auto& last = candidates.back();
    EXPECT_EQ(last.cid, 999999999U);
    EXPECT_EQ(last.verdict.tier, Tier::Eliminated);
    ASSERT_EQ(last.verdict.reasons.size(), 1U);
    EXPECT_NE(last.verdict.reasons[0].find("SMILES"), std::string::npos) << last.verdict.reasons[0];
}

TEST(PipelineRank, RecordsTheFormulaRowsAndRemovals) {
    // A mixture of the two compounds is removed by cleaning, and recorded.
    const std::string csv = pubchem_csv(same_formula()) + "5,\"" + compound_smiles("Acetone") +
                            "." + compound_smiles("Propylene Oxide") + "\"\n";
    const auto ranking =
        rank(formula_of_same_formula(), csv, some_spectrum(), no_groups(), kDefaults);
    ASSERT_TRUE(ranking.ok()) << ranking.error().message;
    const auto imported = arise::io::import_pubchem_csv(csv);
    ASSERT_TRUE(imported.ok());
    const auto hill = arise::core::Formula::parse(formula_of_same_formula());
    ASSERT_TRUE(hill.ok());
    EXPECT_EQ(ranking.value().formula, hill.value().hill());
    EXPECT_EQ(ranking.value().rows, imported.value().rows);
    EXPECT_EQ(ranking.value().removed, imported.value().removed);
    EXPECT_EQ(ranking.value().candidates.size(), imported.value().candidates.size());
}

TEST(PipelineRank, InputErrorsPassThrough) {
    const auto csv = pubchem_csv(same_formula());
    const auto bad_formula = rank("not a formula", csv, some_spectrum(), no_groups(), kDefaults);
    ASSERT_FALSE(bad_formula.ok());
    EXPECT_EQ(bad_formula.error().code, ErrorCode::ParseError);
    const auto bad_csv =
        rank(formula_of_same_formula(), "CID,SMILES\n", some_spectrum(), no_groups(), kDefaults);
    ASSERT_FALSE(bad_csv.ok());
    EXPECT_EQ(bad_csv.error().code, ErrorCode::ParseError);
    auto no_peaks = some_spectrum();
    no_peaks.carbon.text.clear();
    const auto bad_spectrum =
        rank(formula_of_same_formula(), csv, no_peaks, no_groups(), kDefaults);
    ASSERT_FALSE(bad_spectrum.ok());
    EXPECT_TRUE(bad_spectrum.error().message.starts_with("13C spectrum"))
        << bad_spectrum.error().message;
}

// What calling the filter functions directly gives for a search: each cleaned
// candidate's carbon fit, in import order (ascending CID), and their order_by_fit.
struct DirectRanking {
    std::vector<arise::io::Candidate> candidates;
    std::vector<arise::filter::CarbonFit> fits;
    std::vector<std::size_t> order;
};

// check_carbon_fit on the molecule `smiles`; fails the test if that cannot run.
arise::filter::CarbonFit fit_directly(const std::string& smiles,
                                      const std::vector<arise::nmr1d::CarbonAssignment>& peaks) {
    const auto molecule = arise::chem::parse_smiles(smiles);
    const auto fit = molecule.ok()
                         ? arise::filter::check_carbon_fit(molecule.value(), peaks, kDefaults)
                         : arise::core::Result<arise::filter::CarbonFit>{molecule.error()};
    if (!fit.ok()) {
        ADD_FAILURE() << smiles << ": " << fit.error().message;
        return {};
    }
    return fit.value();
}

DirectRanking direct_ranking(const std::string& csv, const SpectraText& spectra) {
    DirectRanking direct;
    const auto imported = arise::io::import_pubchem_csv(csv);
    const auto peaks = read_spectra(spectra);
    if (!imported.ok() || !peaks.ok()) {
        ADD_FAILURE() << "the search or the spectra do not read";
        return direct;
    }
    direct.candidates = imported.value().candidates;
    for (const auto& candidate : direct.candidates) {
        direct.fits.push_back(fit_directly(candidate.smiles, peaks.value()));
    }
    direct.order = arise::filter::order_by_fit(direct.fits);
    return direct;
}

// The ranked form of the direct result's candidate `index`.
arise::pipeline::RankedCandidate expected_candidate(const DirectRanking& direct,
                                                    std::size_t index) {
    const auto& candidate = direct.candidates[index];
    const auto& fit = direct.fits[index];
    constexpr std::string_view kPrefix = "CID ";
    EXPECT_TRUE(candidate.name.starts_with(kPrefix)) << candidate.name;
    return {
        .cid = std::stoull(candidate.name.substr(kPrefix.size())),
        .smiles = candidate.smiles,
        .verdict = fit.verdict,
        .assumed_overlaps = fit.assumed_overlaps,
        .environments_fewest = fit.environments_fewest,
        .environments_most = fit.environments_most,
        .notes = fit.notes,
    };
}

// What rank should give a candidate whose groups are only soft-confirmed: its carbon fit
// (always run), with the group's reasons first and the worse of the two tiers.
arise::pipeline::RankedCandidate expected_soft(
    const arise::pipeline::RankedCandidate& candidate,
    const std::vector<arise::nmr1d::CarbonAssignment>& peaks, const ConfirmedGroups& confirmed) {
    const auto fit = fit_directly(candidate.smiles, peaks);
    const auto groups = group_verdict(candidate.smiles, confirmed);
    auto reasons = groups.reasons;
    reasons.insert(reasons.end(), fit.verdict.reasons.begin(), fit.verdict.reasons.end());
    return {
        .cid = candidate.cid,
        .smiles = candidate.smiles,
        .verdict = {.tier = std::max(groups.tier, fit.verdict.tier), .reasons = reasons},
        .assumed_overlaps = fit.assumed_overlaps,
        .environments_fewest = fit.environments_fewest,
        .environments_most = fit.environments_most,
        .notes = fit.notes,
    };
}

// A group confirmed from 1H evidence is soft (the chemist's decision, 2026-11-22): a
// candidate without it is back-burnered and still gets its carbon fit; its tier is the
// worse of the two, with the group's reason first.
TEST(PipelineRank, ALikelyGroupBackBurnersAndKeepsTheCarbonFit) {
    ASSERT_EQ(compound_formula("Methyl formate"), compound_formula("Acetic Acid"));
    const ConfirmedGroups likely_methoxy{
        .present = {},
        .absent = {},
        .likely_present = {arise::chem::Group::Methoxy},
    };
    const auto spectra = some_spectrum();
    const auto ranking =
        rank(compound_formula("Acetic Acid"), pubchem_csv({"Methyl formate", "Acetic Acid"}),
             spectra, likely_methoxy, kDefaults);
    ASSERT_TRUE(ranking.ok()) << ranking.error().message;
    const auto peaks = read_spectra(spectra);
    ASSERT_TRUE(peaks.ok());
    for (const auto& candidate : ranking.value().candidates) {
        EXPECT_EQ(candidate, expected_soft(candidate, peaks.value(), likely_methoxy))
            << "CID " << candidate.cid;
    }
}

TEST(PipelineRank, ACandidateWithoutTheLikelyGroupIsBackBurneredNotEliminated) {
    const ConfirmedGroups likely_methoxy{
        .present = {},
        .absent = {},
        .likely_present = {arise::chem::Group::Methoxy},
    };
    // Two illustrative 13C peaks: acetic acid's own carbon fit is Active, so only the likely
    // group can move it to the back burner.
    const SpectraText two_peaks{
        .carbon = {.format = TextFormat::Pasted, .text = "21\n177\n"},
        .dept135 = std::nullopt,
        .dept90 = std::nullopt,
    };
    const auto ranking =
        rank(compound_formula("Acetic Acid"), pubchem_csv({"Methyl formate", "Acetic Acid"}),
             two_peaks, likely_methoxy, kDefaults);
    ASSERT_TRUE(ranking.ok()) << ranking.error().message;
    const auto acid = std::find_if(
        ranking.value().candidates.begin(), ranking.value().candidates.end(), [](const auto& c) {
            return c.cid == static_cast<std::uint64_t>(compound_cid("Acetic Acid"));
        });
    ASSERT_NE(acid, ranking.value().candidates.end());
    EXPECT_EQ(acid->verdict.tier, Tier::BackBurner);
    EXPECT_EQ(
        acid->verdict.reasons,
        std::vector<std::string>{"likely present (from 1H): methoxy; the candidate has none"});
    // Methyl formate has the group: its fit stands alone, and it ranks first.
    EXPECT_EQ(ranking.value().candidates.front().cid,
              static_cast<std::uint64_t>(compound_cid("Methyl formate")));
}

// The whole committed search (9,989 cleaned candidates): rank must give exactly the
// order and verdicts of calling the filter functions directly.
TEST(PipelineReference, RanksTheBenchmarkSearchAsTheFilterFunctionsDo) {
    const auto c = benchmark_case();
    const std::string csv = arise::test::join_chunks(arise::test::kBenchmarkPubchemChunks);
    const auto spectra = pasted_spectra(c);
    const auto ranking = rank(c.formula, csv, spectra, no_groups(), kDefaults);
    ASSERT_TRUE(ranking.ok()) << ranking.error().message;
    const auto direct = direct_ranking(csv, spectra);
    const auto& ranked = ranking.value().candidates;
    ASSERT_EQ(ranked.size(), direct.order.size());
    ASSERT_FALSE(ranked.empty());
    for (std::size_t k = 0; k < ranked.size(); ++k) {
        EXPECT_EQ(ranked[k], expected_candidate(direct, direct.order[k])) << "rank " << k;
    }
}

}  // namespace

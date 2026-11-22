// filter::parse_hybridization_windows, suggest_hybridization and check_carbon_fit with
// confirmed counts (plan 3c-b). The windows and shifts in the suggestion tests are
// illustrative inputs, not data; the built-in windows are held to
// data/hybridization_shifts.tsv. Expected hybridizations come from RDKit
// (test_compounds.tsv, compound_hybridizations).

#include "arise/filter/hybridization.hpp"

#include "arise/chem/molecule.hpp"
#include "arise/filter/carbon_fit.hpp"
#include "reference_fixtures.hpp"
#include "test_molecules.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using arise::test::compound_hybridizations;
using arise::test::compound_smiles;

using arise::core::ErrorCode;
using arise::core::Ppm;
using arise::filter::CarbonFit;
using arise::filter::check_carbon_fit;
using arise::filter::ConfirmedCarbonCounts;
using arise::filter::CountRange;
using arise::filter::HybridizationWindows;
using arise::filter::parse_hybridization_windows;
using arise::filter::ShiftRange;
using arise::filter::suggest_hybridization;
using arise::filter::Tier;
using arise::nmr1d::CarbonAssignment;
using arise::nmr1d::CarbonTypeSet;

std::string measured_text() {
    return arise::test::join_chunks(arise::test::kHybridizationShiftsChunks);
}

// The measured table's data lines.
std::vector<std::string> measured_rows() {
    std::vector<std::string> rows;
    const auto text = measured_text();
    std::size_t start = 0;
    while (start < text.size()) {
        const auto end = text.find('\n', start);
        auto line = text.substr(start, end == std::string::npos ? end : end - start);
        if (!line.empty() && line.front() != '#') {
            rows.push_back(std::move(line));
        }
        start = end == std::string::npos ? text.size() : end + 1;
    }
    return rows;
}

// Field `i` (0-based) of a tab-separated row.
std::string field(const std::string& row, std::size_t i) {
    std::size_t start = 0;
    for (std::size_t k = 0; k < i; ++k) {
        start = row.find('\t', start) + 1;
    }
    return row.substr(start, row.find('\t', start) - start);
}

// `row` with field `i` replaced by `value`.
std::string with_field(const std::string& row, std::size_t i, std::string_view value) {
    std::size_t start = 0;
    for (std::size_t k = 0; k < i; ++k) {
        start = row.find('\t', start) + 1;
    }
    const auto end = row.find('\t', start);
    return row.substr(0, start) + std::string{value} +
           (end == std::string::npos ? "" : row.substr(end));
}

std::string table(const std::vector<std::string>& rows) {
    std::string text = "# comment\n";
    for (const auto& row : rows) {
        text += row + "\n";
    }
    return text;
}

void expect_parse_error(const std::string& text, std::string_view needle) {
    const auto result = parse_hybridization_windows(text);
    ASSERT_FALSE(result) << text;
    EXPECT_EQ(result.error().code, ErrorCode::ParseError);
    EXPECT_NE(result.error().message.find(needle), std::string::npos) << result.error().message;
}

// ---- Windows ---------------------------------------------------------------------

TEST(HybridizationWindows, BuiltInWindowsAreTheMeasuredTable) {
    const auto built_in = arise::filter::measured_hybridization_windows();
    ASSERT_TRUE(built_in) << built_in.error().message;
    const auto parsed = parse_hybridization_windows(measured_text());
    ASSERT_TRUE(parsed) << parsed.error().message;
    EXPECT_EQ(built_in.value(), parsed.value());
}

const ShiftRange& range_named(const HybridizationWindows& windows, std::string_view name) {
    if (name == "sp") {
        return windows.sp;
    }
    return name == "sp2" ? windows.sp2 : windows.sp3;
}

TEST(HybridizationWindows, RangesRunFromThePointOneToTheNinetyNinePointNinePercentile) {
    const auto parsed = parse_hybridization_windows(measured_text());
    ASSERT_TRUE(parsed) << parsed.error().message;
    for (const auto& row : measured_rows()) {
        const auto name = field(row, 0);
        const auto& range = range_named(parsed.value(), name);
        EXPECT_EQ(range.low.value(), std::stod(field(row, 4))) << row;
        EXPECT_EQ(range.high.value(), std::stod(field(row, 8))) << row;
    }
    EXPECT_EQ(measured_rows().size(), 3U);
}

TEST(HybridizationWindows, RowOrderDoesNotMatter) {
    auto rows = measured_rows();
    std::reverse(rows.begin(), rows.end());
    const auto reversed = parse_hybridization_windows(table(rows));
    const auto parsed = parse_hybridization_windows(measured_text());
    ASSERT_TRUE(reversed && parsed);
    EXPECT_EQ(reversed.value(), parsed.value());
}

TEST(HybridizationWindows, RejectsMalformedTables) {
    const auto rows = measured_rows();
    ASSERT_EQ(rows.size(), 3U);
    expect_parse_error(table({rows[0], rows[1]}), "sp3");  // a hybridization missing
    expect_parse_error(table({rows[0], rows[1], rows[2], rows[2]}), "line 5");  // repeated
    expect_parse_error(table({with_field(rows[0], 0, "sp4"), rows[1], rows[2]}), "line 2");
    expect_parse_error(table({rows[0].substr(0, rows[0].rfind('\t')), rows[1], rows[2]}),
                       "line 2");  // 9 fields
    expect_parse_error(table({with_field(rows[0], 1, "0"), rows[1], rows[2]}), "line 2");
    expect_parse_error(table({with_field(rows[0], 2, "1.5"), rows[1], rows[2]}), "line 2");
    expect_parse_error(table({with_field(rows[0], 4, "low"), rows[1], rows[2]}), "line 2");
    expect_parse_error(table({with_field(rows[0], 4, "nan"), rows[1], rows[2]}), "line 2");
    // 0.1st percentile above the 99.9th.
    expect_parse_error(
        table({
            with_field(with_field(rows[1], 4, field(rows[1], 8)), 8, field(rows[1], 4)),
            rows[0],
            rows[2],
        }),
        "line 2");
    expect_parse_error("# only comments\n", "sp");
}

// ---- Suggestions -----------------------------------------------------------------

// Illustrative ranges: sp3 0-110, sp2 100-200 (overlapping 100-110), sp 50-150.
HybridizationWindows windows() {
    return {
        .sp = {.low = Ppm{50.0}, .high = Ppm{150.0}},
        .sp2 = {.low = Ppm{100.0}, .high = Ppm{200.0}},
        .sp3 = {.low = Ppm{0.0}, .high = Ppm{110.0}},
    };
}

CarbonAssignment peak(double shift) {
    return {.shift = Ppm{shift}, .possible = CarbonTypeSet::all()};
}

TEST(SuggestHybridization, CountsPeaksInsideEachRangeAndInsideOnlyOne) {
    const std::vector peaks{peak(20.0), peak(30.0), peak(105.0), peak(150.0)};
    const auto result = suggest_hybridization(peaks, windows());
    ASSERT_TRUE(result) << result.error().message;
    EXPECT_EQ(result.value().sp2, (CountRange{.min = 1, .max = 2}));
    EXPECT_EQ(result.value().sp3, (CountRange{.min = 2, .max = 3}));
    ASSERT_EQ(result.value().peaks.size(), 4U);
    EXPECT_FALSE(result.value().peaks[0].could_be_sp2);
    EXPECT_TRUE(result.value().peaks[0].could_be_sp3);
    EXPECT_TRUE(result.value().peaks[2].could_be_sp2);
    EXPECT_TRUE(result.value().peaks[2].could_be_sp3);
    EXPECT_EQ(result.value().peaks[3].shift, Ppm{150.0});
}

TEST(SuggestHybridization, RangeBoundsAreInclusive) {
    const std::vector peaks{peak(0.0), peak(200.0)};
    const auto result = suggest_hybridization(peaks, windows());
    ASSERT_TRUE(result);
    EXPECT_EQ(result.value().sp3, (CountRange{.min = 1, .max = 1}));
    EXPECT_EQ(result.value().sp2, (CountRange{.min = 1, .max = 1}));
}

TEST(SuggestHybridization, AlwaysNotesTheSpRangeAndThatSpIsNotSuggested) {
    const auto result = suggest_hybridization(std::vector{peak(20.0)}, windows());
    ASSERT_TRUE(result);
    ASSERT_FALSE(result.value().notes.empty());
    const auto& note = result.value().notes.front();
    EXPECT_NE(note.find("50"), std::string::npos) << note;
    EXPECT_NE(note.find("150"), std::string::npos) << note;
    EXPECT_NE(note.find("sp carbons"), std::string::npos) << note;
}

TEST(SuggestHybridization, TypedSpPeaksLowerTheMinimumsButNotBelowZero) {
    const std::vector peaks{peak(20.0), peak(30.0), peak(150.0)};
    const auto result = suggest_hybridization(peaks, windows(), 1);
    ASSERT_TRUE(result);
    EXPECT_EQ(result.value().sp3, (CountRange{.min = 1, .max = 2}));
    EXPECT_EQ(result.value().sp2, (CountRange{.min = 0, .max = 1}));
    const auto& notes = result.value().notes;
    EXPECT_TRUE(std::any_of(notes.begin(), notes.end(), [](const std::string& n) {
        return n.find("1 typed sp peak") != std::string::npos;
    }));
    const auto zero = suggest_hybridization(peaks, windows(), 0);
    ASSERT_TRUE(zero);
    EXPECT_EQ(zero.value().sp2.min, 1);
}

TEST(SuggestHybridization, PeaksOutsideBothRangesCountTowardNeitherAndAreNoted) {
    const std::vector peaks{peak(-20.0), peak(250.0), peak(20.0)};
    const auto result = suggest_hybridization(peaks, windows());
    ASSERT_TRUE(result);
    EXPECT_EQ(result.value().sp3, (CountRange{.min = 1, .max = 1}));
    EXPECT_EQ(result.value().sp2, (CountRange{.min = 0, .max = 0}));
    const auto& notes = result.value().notes;
    EXPECT_TRUE(std::any_of(notes.begin(), notes.end(), [](const std::string& n) {
        return n.find("2 peaks lie outside both") != std::string::npos;
    }));
}

TEST(SuggestHybridization, PeakOrderDoesNotChangeTheCounts) {
    std::vector peaks{peak(150.0), peak(105.0), peak(30.0), peak(20.0)};
    const auto forward = suggest_hybridization(peaks, windows());
    std::reverse(peaks.begin(), peaks.end());
    const auto backward = suggest_hybridization(peaks, windows());
    ASSERT_TRUE(forward && backward);
    EXPECT_EQ(forward.value().sp2, backward.value().sp2);
    EXPECT_EQ(forward.value().sp3, backward.value().sp3);
    EXPECT_EQ(forward.value().notes, backward.value().notes);
}

TEST(SuggestHybridization, RejectsANegativeTypedSpCount) {
    const auto result = suggest_hybridization(std::vector{peak(20.0)}, windows(), -1);
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, ErrorCode::InvalidArgument);
}

TEST(SuggestHybridization, AnEmptyListSuggestsNothing) {
    const auto result = suggest_hybridization(std::vector<CarbonAssignment>{}, windows());
    ASSERT_TRUE(result);
    EXPECT_EQ(result.value().sp2, CountRange{});
    EXPECT_EQ(result.value().sp3, CountRange{});
}

// ---- Confirmed counts in check_carbon_fit ----------------------------------------

using Field = std::optional<int> ConfirmedCarbonCounts::*;
constexpr Field kSp = &ConfirmedCarbonCounts::sp;
constexpr Field kSp2 = &ConfirmedCarbonCounts::sp2;
constexpr Field kSp3 = &ConfirmedCarbonCounts::sp3;
constexpr Field kC = &ConfirmedCarbonCounts::c;
constexpr Field kCh = &ConfirmedCarbonCounts::ch;
constexpr Field kCh2 = &ConfirmedCarbonCounts::ch2;
constexpr Field kCh3 = &ConfirmedCarbonCounts::ch3;

// Confirmed counts with only the listed fields set.
ConfirmedCarbonCounts counts(std::initializer_list<std::pair<Field, int>> values) {
    ConfirmedCarbonCounts confirmed;
    for (const auto& [field, value] : values) {
        confirmed.*field = value;
    }
    return confirmed;
}

arise::chem::Molecule molecule(std::string_view name) {
    auto parsed = arise::chem::parse_smiles(compound_smiles(name));
    EXPECT_TRUE(parsed) << name;
    return parsed ? std::move(parsed).value() : arise::chem::Molecule{};
}

// RDKit's number of carbons of hybridization `code` (1 sp, 2 sp2, 3 sp3) in `name`.
int rdkit_carbons(std::string_view name, int code) {
    const auto codes = compound_hybridizations(name);
    return static_cast<int>(std::count(codes.begin(), codes.end(), code));
}

CarbonFit fit(std::string_view name, const ConfirmedCarbonCounts& confirmed) {
    const auto result = check_carbon_fit(molecule(name), std::vector{peak(20.0)}, confirmed);
    EXPECT_TRUE(result) << name << ": " << (result ? "" : result.error().message);
    return result ? result.value() : CarbonFit{};
}

TEST(ConfirmedCarbonCounts, CountsTheCandidateHasKeepItActive) {
    // Ethyl acetate has no symmetry: each carbon is its own environment.
    const auto result = fit("Ethyl Acetate", counts({
                                                 {kSp2, rdkit_carbons("Ethyl Acetate", 2)},
                                                 {kSp3, rdkit_carbons("Ethyl Acetate", 3)},
                                             }));
    EXPECT_EQ(result.verdict.tier, Tier::Active);
    EXPECT_TRUE(result.verdict.reasons.empty());
}

TEST(ConfirmedCarbonCounts, UnconfirmedCountsChangeNothing) {
    const auto plain = check_carbon_fit(molecule("Ethyl Acetate"), std::vector{peak(20.0)});
    ASSERT_TRUE(plain);
    EXPECT_EQ(fit("Ethyl Acetate", {}), plain.value());
}

TEST(ConfirmedCarbonCounts, AShortfallGoesToTheBackBurnerWithItsReason) {
    const int sp2 = rdkit_carbons("Ethyl Acetate", 2);
    const auto result = fit("Ethyl Acetate", counts({{kSp2, sp2 + 1}}));
    EXPECT_EQ(result.verdict.tier, Tier::BackBurner);
    EXPECT_EQ(result.verdict.reasons,
              std::vector<std::string>{"confirmed at least " + std::to_string(sp2 + 1) +
                                       " sp2 peaks; the candidate has at most " +
                                       std::to_string(sp2) + " sp2 carbon environments"});
}

TEST(ConfirmedCarbonCounts, ATypedSpCountIsCheckedToo) {
    EXPECT_EQ(fit("Ethyl Acetate", counts({{kSp, 1}})).verdict.tier, Tier::BackBurner);
    ASSERT_EQ(rdkit_carbons("Acetonitrile", 1), 1);
    EXPECT_EQ(fit("Acetonitrile", counts({{kSp, 1}})).verdict.tier, Tier::Active);
}

TEST(ConfirmedCarbonCounts, EnvironmentsAreCountedNotCarbons) {
    // Benzene's six sp2 carbons are one environment.
    ASSERT_EQ(rdkit_carbons("Benzene", 2), 6);
    EXPECT_EQ(fit("Benzene", counts({{kSp2, 1}})).verdict.tier, Tier::Active);
    EXPECT_EQ(fit("Benzene", counts({{kSp2, 2}})).verdict.tier, Tier::BackBurner);
}

TEST(ConfirmedCarbonCounts, AnUnassignedChargedCarbonStandsInForAnyOneType) {
    // Methyl isocyanide: a charged carbon (unassigned) and an sp3 methyl.
    EXPECT_EQ(fit("Methyl isocyanide", counts({{kSp, 1}, {kSp3, 1}})).verdict.tier, Tier::Active);
    EXPECT_EQ(fit("Methyl isocyanide", counts({{kSp3, 2}})).verdict.tier, Tier::Active);
    EXPECT_EQ(fit("Methyl isocyanide", counts({{kSp2, 1}, {kSp3, 1}})).verdict.tier, Tier::Active);
    // ...but only once: three types cannot all be met by two environments.
    const auto result = fit("Methyl isocyanide", counts({{kSp, 1}, {kSp2, 1}, {kSp3, 1}}));
    EXPECT_EQ(result.verdict.tier, Tier::BackBurner);
    ASSERT_EQ(result.verdict.reasons.size(), 1U);
    EXPECT_NE(result.verdict.reasons.front().find("together"), std::string::npos)
        << result.verdict.reasons.front();
}

TEST(ConfirmedCarbonCounts, TypedDeptCountsAreCheckedAgainstCarbonTypes) {
    // Ethyl acetate: two CH3, one CH2, one C.
    EXPECT_EQ(fit("Ethyl Acetate", counts({{kC, 1}, {kCh2, 1}, {kCh3, 2}})).verdict.tier,
              Tier::Active);
    EXPECT_EQ(fit("Ethyl Acetate", counts({{kCh3, 3}})).verdict.tier, Tier::BackBurner);
    EXPECT_EQ(fit("Ethyl Acetate", counts({{kCh, 1}})).verdict.tier, Tier::BackBurner);
}

TEST(ConfirmedCarbonCounts, MethaneHasNoDeptTypeSoStandsInForAny) {
    EXPECT_EQ(fit("Methane", counts({{kCh3, 1}})).verdict.tier, Tier::Active);
    EXPECT_EQ(fit("Methane", counts({{kCh3, 1}, {kCh2, 1}})).verdict.tier, Tier::BackBurner);
}

TEST(ConfirmedCarbonCounts, KeepsThePeakCountChecksReasonsFirst) {
    // Ethyl acetate has four environments; five peaks are already a shortfall.
    const std::vector peaks{peak(10.0), peak(20.0), peak(30.0), peak(40.0), peak(50.0)};
    const auto result = check_carbon_fit(molecule("Ethyl Acetate"), peaks, counts({{kSp, 1}}));
    ASSERT_TRUE(result);
    EXPECT_EQ(result.value().verdict.tier, Tier::BackBurner);
    ASSERT_EQ(result.value().verdict.reasons.size(), 2U);
    EXPECT_NE(result.value().verdict.reasons[0].find("13C list has 5 peaks"), std::string::npos);
    EXPECT_NE(result.value().verdict.reasons[1].find("sp peaks"), std::string::npos);
}

TEST(ConfirmedCarbonCounts, RejectsNegativeCounts) {
    for (const auto& confirmed : {counts({{kSp, -1}}), counts({{kSp3, -1}}), counts({{kCh, -2}})}) {
        const auto result =
            check_carbon_fit(molecule("Ethyl Acetate"), std::vector{peak(20.0)}, confirmed);
        ASSERT_FALSE(result);
        EXPECT_EQ(result.error().code, ErrorCode::InvalidArgument);
    }
    const auto empty = check_carbon_fit(molecule("Ethyl Acetate"), std::vector<CarbonAssignment>{},
                                        ConfirmedCarbonCounts{});
    ASSERT_FALSE(empty);
    EXPECT_EQ(empty.error().code, ErrorCode::InvalidArgument);
}

}  // namespace

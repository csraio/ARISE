// Group suggestions and confirmed-group elimination (plan 3b, steps 3 and 4). These are
// engineering cases (bounds, ordering, errors) on small hand-made windows that make no
// chemical claim; the real windows are checked against nmrshiftdb2's hold-out half in
// holdout_test.cpp. Molecules here are re-checked by RDKit (spec_rdkit_oracle.tsv).
#include "arise/filter/groups.hpp"

#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "arise/core/error.hpp"
#include "arise/core/units.hpp"
#include "arise/nmr1d/dept.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::chem::Group;
using arise::core::ErrorCode;
using arise::core::Ppm;
using arise::filter::ConfirmedGroups;
using arise::filter::GroupSuggestion;
using arise::filter::HallmarkWindow;
using arise::filter::parse_hallmark_windows;
using arise::filter::suggest_groups;
using arise::filter::Tier;
using arise::nmr1d::CarbonAssignment;
using arise::nmr1d::CarbonType;
using arise::nmr1d::CarbonTypeSet;

HallmarkWindow window(Group group, std::string hallmark, CarbonTypeSet types, double low,
                      double high) {
    return {
        .group = group,
        .hallmark = std::move(hallmark),
        .types = types,
        .low = Ppm{low},
        .high = Ppm{high},
        .n = 1,
    };
}

CarbonAssignment peak(double shift, CarbonTypeSet possible = CarbonTypeSet::all()) {
    return {.shift = Ppm{shift}, .possible = possible};
}

std::vector<Group> suggested(const std::vector<GroupSuggestion>& suggestions) {
    std::vector<Group> groups;
    groups.reserve(suggestions.size());
    for (const auto& suggestion : suggestions) {
        groups.push_back(suggestion.group);
    }
    return groups;
}

// ---- parse_hallmark_windows ------------------------------------------------------------

std::vector<std::string_view> split(std::string_view text, char separator) {
    std::vector<std::string_view> parts;
    for (auto at = text.find(separator); at != std::string_view::npos; at = text.find(separator)) {
        parts.push_back(text.substr(0, at));
        text.remove_prefix(at + 1);
    }
    parts.push_back(text);
    return parts;
}

// The table's rows (not comments or blank lines), split into fields.
std::vector<std::vector<std::string_view>> table_rows(std::string_view text) {
    std::vector<std::vector<std::string_view>> rows;
    for (const auto line : split(text, '\n')) {
        if (!line.empty() && line.front() != '#') {
            rows.push_back(split(line, '\t'));
        }
    }
    return rows;
}

// A window's group, hallmark, types, n, low and high, as the table's columns 1-4, 6 and 10.
std::vector<std::string> window_fields(const HallmarkWindow& w) {
    return {
        std::string{arise::chem::to_string(w.group)},
        w.hallmark,
        arise::nmr1d::to_string(w.types),
        std::to_string(w.n),
        std::to_string(w.low.value()),
        std::to_string(w.high.value()),
    };
}

std::vector<std::string> row_fields(const std::vector<std::string_view>& row) {
    if (row.size() != 11) {
        return {"malformed row"};
    }
    return {
        std::string{row[0]},
        std::string{row[1]},
        std::string{row[2]},
        std::string{row[3]},
        std::to_string(std::stod(std::string{row[5]})),
        std::to_string(std::stod(std::string{row[9]})),
    };
}

TEST(HallmarkWindows, ReadsEveryRowOfTheMeasuredTable) {
    // The measured table (tools/measure_group_shifts): one window per row, from the 1st
    // (column 6) to the 99th percentile (column 10), with the row's group, hallmark,
    // types and n.
    const auto text = arise::test::join_chunks(arise::test::kGroupShiftsChunks);
    const auto windows = parse_hallmark_windows(text);
    ASSERT_TRUE(windows.ok()) << windows.error().message;
    const auto rows = table_rows(text);
    ASSERT_EQ(windows.value().size(), rows.size());
    EXPECT_GT(rows.size(), 20U);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        EXPECT_EQ(window_fields(windows.value()[i]), row_fields(rows[i]));
    }
}

constexpr std::string_view kRow = "ketone\tC=O carbon\t{C}\t10\t1.0\t2.0\t3.0\t4.0\t5.0\t6.0\t7.0";

TEST(HallmarkWindows, SkipsCommentsAndBlankLinesAndAcceptsCrlf) {
    const auto windows = parse_hallmark_windows("# header\n\n" + std::string{kRow} + "\r\n");
    ASSERT_TRUE(windows.ok()) << windows.error().message;
    ASSERT_EQ(windows.value().size(), 1U);
    auto expected = window(Group::Ketone, "C=O carbon", {CarbonType::C}, 2.0, 6.0);
    expected.n = 10;
    EXPECT_EQ(windows.value()[0], expected);
}

void expect_parse_error(const std::string& text, std::string_view fragment) {
    const auto windows = parse_hallmark_windows(text);
    ASSERT_FALSE(windows.ok()) << "accepted: " << text;
    EXPECT_EQ(windows.error().code, ErrorCode::ParseError);
    EXPECT_NE(windows.error().message.find(fragment), std::string::npos) << windows.error().message;
}

TEST(HallmarkWindows, RejectsMalformedTables) {
    const auto row = [](std::string_view group, std::string_view types, std::string_view n,
                        std::string_view p1, std::string_view p99) {
        return std::string{group} + "\tC=O carbon\t" + std::string{types} + "\t" + std::string{n} +
               "\t1.0\t" + std::string{p1} + "\t3.0\t4.0\t5.0\t" + std::string{p99} + "\t7.0\n";
    };
    expect_parse_error("# only comments\n", "no rows");
    expect_parse_error("", "no rows");
    expect_parse_error("ketone\tC=O carbon\t{C}\t10\n", "line 1");
    expect_parse_error("# x\n" + row("ketones", "{C}", "10", "2.0", "6.0"), "line 2");
    expect_parse_error(row("ketones", "{C}", "10", "2.0", "6.0"), "unknown group");
    expect_parse_error(row("ketone", "{}", "10", "2.0", "6.0"), "DEPT types");
    expect_parse_error(row("ketone", "{C, CH4}", "10", "2.0", "6.0"), "DEPT types");
    expect_parse_error(row("ketone", "C", "10", "2.0", "6.0"), "DEPT types");
    expect_parse_error(row("ketone", "{C}", "0", "2.0", "6.0"), "n must");
    expect_parse_error(row("ketone", "{C}", "ten", "2.0", "6.0"), "n must");
    expect_parse_error(row("ketone", "{C}", "10", "2,0", "6.0"), "number");
    expect_parse_error(row("ketone", "{C}", "10", "2.0", "1e999"), "number");
    expect_parse_error(row("ketone", "{C}", "10", "6.5", "6.0"), "above");
    expect_parse_error(
        row("ketone", "{C}", "10", "2.0", "6.0") + row("ketone", "{C}", "10", "2.0", "6.0"),
        "twice");
}

TEST(HallmarkWindows, AcceptsAOnePointWindow) {
    const auto windows =
        parse_hallmark_windows("ketone\tC=O carbon\t{C}\t1\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\n");
    ASSERT_TRUE(windows.ok()) << windows.error().message;
    EXPECT_EQ(windows.value()[0].low, windows.value()[0].high);
}

// ---- suggest_groups --------------------------------------------------------------------

TEST(SuggestGroups, PeakInsideAWindowSuggestsItsGroupWithEvidence) {
    const std::vector windows{window(Group::Ketone, "C=O carbon", {CarbonType::C}, 190, 220)};
    const std::vector peaks{peak(20.0), peak(205.0)};
    const auto suggestions = suggest_groups(peaks, windows);
    ASSERT_EQ(suggestions.size(), 1U);
    EXPECT_EQ(suggestions[0].group, Group::Ketone);
    ASSERT_EQ(suggestions[0].evidence.size(), 1U);
    EXPECT_EQ(suggestions[0].evidence[0].peak, 1U);
    EXPECT_EQ(suggestions[0].evidence[0].shift, Ppm{205.0});
    EXPECT_EQ(suggestions[0].evidence[0].window, windows[0]);
}

TEST(SuggestGroups, WindowBoundsAreInclusive) {
    const std::vector windows{window(Group::Nitrile, "C#N carbon", {CarbonType::C}, 110, 120)};
    EXPECT_EQ(suggest_groups(std::vector{peak(110.0)}, windows).size(), 1U);
    EXPECT_EQ(suggest_groups(std::vector{peak(120.0)}, windows).size(), 1U);
    EXPECT_TRUE(suggest_groups(std::vector{peak(109.99)}, windows).empty());
    EXPECT_TRUE(suggest_groups(std::vector{peak(120.01)}, windows).empty());
}

TEST(SuggestGroups, PeakMustShareADeptTypeWithTheWindow) {
    const std::vector windows{
        window(Group::Aldehyde, "C=O carbon", {CarbonType::CH}, 190, 205),
        window(Group::Ketone, "C=O carbon", {CarbonType::C}, 190, 220),
    };
    EXPECT_EQ(suggested(suggest_groups(std::vector{peak(200.0, {CarbonType::CH})}, windows)),
              std::vector{Group::Aldehyde});
    EXPECT_EQ(suggested(suggest_groups(std::vector{peak(200.0, {CarbonType::C})}, windows)),
              std::vector{Group::Ketone});
    // DEPT-135 alone: CH or CH3.
    EXPECT_EQ(suggested(suggest_groups(std::vector{peak(200.0, {CarbonType::CH, CarbonType::CH3})},
                                       windows)),
              std::vector{Group::Aldehyde});
    EXPECT_TRUE(suggest_groups(std::vector{peak(200.0, {CarbonType::CH2})}, windows).empty());
    EXPECT_TRUE(suggest_groups(std::vector{peak(200.0, CarbonTypeSet{})}, windows).empty());
}

TEST(SuggestGroups, UnknownDeptTypeMatchesOnShiftAlone) {
    const std::vector windows{
        window(Group::Aldehyde, "C=O carbon", {CarbonType::CH}, 190, 205),
        window(Group::Ketone, "C=O carbon", {CarbonType::C}, 190, 220),
    };
    EXPECT_EQ(suggested(suggest_groups(std::vector{peak(200.0)}, windows)),
              (std::vector{Group::Aldehyde, Group::Ketone}));
}

TEST(SuggestGroups, OverlappingWindowsSuggestEveryGroupInGroupOrder) {
    // Windows deliberately listed out of group order.
    const std::vector windows{
        window(Group::CarboxylicEster, "C=O carbon", {CarbonType::C}, 160, 178),
        window(Group::Carboxamide, "C=O carbon", {CarbonType::C}, 155, 180),
        window(Group::CarboxylicAcid, "C=O carbon", {CarbonType::C}, 165, 185),
    };
    EXPECT_EQ(suggested(suggest_groups(std::vector{peak(170.0)}, windows)),
              (std::vector{Group::CarboxylicAcid, Group::CarboxylicEster, Group::Carboxamide}));
}

TEST(SuggestGroups, EvidenceListsEveryPeakAndWindowByPeakThenWindow) {
    const std::vector windows{
        window(Group::Ether, "O-C carbon (sp3)", {CarbonType::CH2, CarbonType::CH3}, 50, 110),
        window(Group::Ether, "O-C carbon (unsaturated)", {CarbonType::C, CarbonType::CH}, 100, 180),
    };
    const std::vector peaks{peak(160.0), peak(105.0), peak(30.0), peak(60.0)};
    const auto suggestions = suggest_groups(peaks, windows);
    ASSERT_EQ(suggestions.size(), 1U);
    const auto& evidence = suggestions[0].evidence;
    ASSERT_EQ(evidence.size(), 4U);
    EXPECT_EQ(evidence[0].peak, 0U);
    EXPECT_EQ(evidence[0].window.hallmark, "O-C carbon (unsaturated)");
    EXPECT_EQ(evidence[1].peak, 1U);
    EXPECT_EQ(evidence[1].window.hallmark, "O-C carbon (sp3)");
    EXPECT_EQ(evidence[2].peak, 1U);
    EXPECT_EQ(evidence[2].window.hallmark, "O-C carbon (unsaturated)");
    EXPECT_EQ(evidence[3].peak, 3U);
    EXPECT_EQ(evidence[3].window.hallmark, "O-C carbon (sp3)");
}

TEST(SuggestGroups, NothingToMatchSuggestsNothing) {
    const std::vector windows{window(Group::Ketone, "C=O carbon", {CarbonType::C}, 190, 220)};
    EXPECT_TRUE(suggest_groups(std::vector<CarbonAssignment>{}, windows).empty());
    EXPECT_TRUE(suggest_groups(std::vector{peak(200.0)}, std::vector<HallmarkWindow>{}).empty());
}

// ---- check_confirmed_groups ------------------------------------------------------------

arise::filter::Verdict verdict(std::string_view smiles, const ConfirmedGroups& confirmed) {
    const auto molecule = arise::chem::parse_smiles(smiles);
    EXPECT_TRUE(molecule.ok()) << smiles;
    const auto result = arise::filter::check_confirmed_groups(molecule.value(), confirmed);
    EXPECT_TRUE(result.ok()) << result.error().message;
    return result.value();
}

TEST(ConfirmedGroups, ConfirmedEsterKeepsEstersAndEliminatesTheRest) {
    // The plan's end-to-end case (C4H8O2, "carboxylic ester" confirmed).
    const ConfirmedGroups ester{.present = {Group::CarboxylicEster}, .absent = {}};
    EXPECT_EQ(verdict("CCOC(C)=O", ester).tier, Tier::Active);  // ethyl acetate
    EXPECT_EQ(verdict("CCC(=O)OC", ester).tier, Tier::Active);  // methyl propanoate
    const auto methoxyacetone = verdict("COCC(C)=O", ester);
    EXPECT_EQ(methoxyacetone.tier, Tier::Eliminated);
    EXPECT_EQ(
        methoxyacetone.reasons,
        std::vector<std::string>{"confirmed present: carboxylic ester; the candidate has none"});
}

TEST(ConfirmedGroups, GroupRuledOutEliminatesCandidatesThatHaveIt) {
    const ConfirmedGroups no_ketone{.present = {}, .absent = {Group::Ketone}};
    EXPECT_EQ(verdict("CCOC(C)=O", no_ketone).tier, Tier::Active);
    EXPECT_EQ(verdict("COCC(C)=O", no_ketone).reasons,
              std::vector<std::string>{"confirmed absent: ketone; the candidate has 1"});
    EXPECT_EQ(verdict("CC(=O)CC(C)=O", no_ketone).reasons,
              std::vector<std::string>{"confirmed absent: ketone; the candidate has 2"});
}

TEST(ConfirmedGroups, EveryViolationIsAReasonInTheOrderGiven) {
    const ConfirmedGroups confirmed{
        .present = {Group::Alcohol, Group::CarboxylicEster},
        .absent = {Group::Ether, Group::Ketone},
    };
    const auto result = verdict("COCC(C)=O", confirmed);
    EXPECT_EQ(result.tier, Tier::Eliminated);
    EXPECT_EQ(result.reasons, (std::vector<std::string>{
                                  "confirmed present: alcohol; the candidate has none",
                                  "confirmed present: carboxylic ester; the candidate has none",
                                  "confirmed absent: ether; the candidate has 1",
                                  "confirmed absent: ketone; the candidate has 1",
                              }));
}

TEST(ConfirmedGroups, NothingConfirmedKeepsEveryCandidateActive) {
    EXPECT_EQ(verdict("COCC(C)=O", {}), arise::filter::Verdict{});
}

TEST(ConfirmedGroups, GroupConfirmedBothPresentAndAbsentIsAnError) {
    const auto molecule = arise::chem::parse_smiles("CCOC(C)=O");
    ASSERT_TRUE(molecule.ok());
    const auto result = arise::filter::check_confirmed_groups(
        molecule.value(),
        {.present = {Group::CarboxylicEster}, .absent = {Group::Ketone, Group::CarboxylicEster}});
    ASSERT_FALSE(result.ok());
    EXPECT_EQ(result.error().code, ErrorCode::InvalidArgument);
    EXPECT_NE(result.error().message.find("carboxylic ester"), std::string::npos)
        << result.error().message;
}

TEST(Tier, Names) {
    EXPECT_EQ(arise::filter::to_string(Tier::Active), "active");
    EXPECT_EQ(arise::filter::to_string(Tier::BackBurner), "back-burner");
    EXPECT_EQ(arise::filter::to_string(Tier::Eliminated), "eliminated");
}

}  // namespace

// filter/proton_groups (the 1H plan): reading the measured 1H windows and suggesting groups
// from 1H peaks. Windows and shifts below are illustrative inputs for the rules, not
// measured values; the measured table is checked against data/group_proton_shifts.tsv.

#include "arise/filter/proton_groups.hpp"

#include "arise/chem/perception.hpp"
#include "arise/core/error.hpp"
#include "arise/core/units.hpp"
#include "arise/io/peak_csv.hpp"
#include "arise/nmr1d/dept.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::chem::Group;
using arise::core::ErrorCode;
using arise::core::Ppm;
using arise::filter::parse_proton_windows;
using arise::filter::proton_count;
using arise::filter::ProtonEvidence;
using arise::filter::ProtonOptions;
using arise::filter::ProtonSuggestion;
using arise::filter::ProtonWindow;
using arise::filter::suggest_groups_from_protons;
using arise::io::ObservedPeak;
using arise::nmr1d::CarbonType;

// ---- reading the table ----------------------------------------------------------------

// Two rows in data/group_proton_shifts.tsv's layout (21 tab-separated columns).
constexpr std::string_view kTwoRows =
    "# comment\n"
    "methoxy\tO-CH3 carbon\tCH3\t1082\t807\t2.28\t3.04\t3.28\t3.77\t4.00\t4.31\t6.20\t3.04\t"
    "4.31\t1172\t843\t0.719\t843\t6020\t0.140\t1\n"
    "ether\tO-C carbon (sp3)\tCH2\t900\t534\t2.10\t2.70\t3.00\t3.60\t4.50\t5.99\t6.40\t2.70\t"
    "5.99\t3000\t864\t0.288\t1330\t6020\t0.221\t-\n";

TEST(ProtonWindows, ReadsTheTableLayout) {
    const auto windows = parse_proton_windows(kTwoRows);
    ASSERT_TRUE(windows.ok()) << windows.error().message;
    ASSERT_EQ(windows.value().size(), 2U);
    const auto& methoxy = windows.value()[0];
    EXPECT_EQ(methoxy.group, Group::Methoxy);
    EXPECT_EQ(methoxy.hallmark, "O-CH3 carbon");
    EXPECT_EQ(methoxy.type, CarbonType::CH3);
    EXPECT_EQ(methoxy.low, Ppm{3.04});
    EXPECT_EQ(methoxy.high, Ppm{4.31});
    EXPECT_EQ(methoxy.n, 1082U);
    EXPECT_EQ(methoxy.molecules, 807U);
    EXPECT_DOUBLE_EQ(methoxy.precision, 0.719);
    EXPECT_DOUBLE_EQ(methoxy.base_rate, 0.140);
    EXPECT_TRUE(methoxy.round1);
    EXPECT_EQ(windows.value()[1].type, CarbonType::CH2);
    EXPECT_FALSE(windows.value()[1].round1);
}

void expect_parse_error(const std::string& text, std::string_view fragment) {
    const auto result = parse_proton_windows(text);
    ASSERT_FALSE(result.ok()) << "accepted:\n" << text;
    EXPECT_EQ(result.error().code, ErrorCode::ParseError);
    EXPECT_NE(result.error().message.find(fragment), std::string::npos) << result.error().message;
}

// kTwoRows' first data row with column `column` replaced by `value`.
std::string with_field(std::size_t column, std::string_view value) {
    const std::string text{kTwoRows};
    const auto start = text.find('\n') + 1;
    const auto end = text.find('\n', start);
    std::string row = text.substr(start, end - start);
    std::vector<std::string> fields;
    std::size_t from = 0;
    while (true) {
        const auto tab = row.find('\t', from);
        fields.push_back(row.substr(from, tab - from));
        if (tab == std::string::npos) {
            break;
        }
        from = tab + 1;
    }
    fields[column] = std::string{value};
    std::string joined;
    for (const auto& field : fields) {
        joined += (joined.empty() ? "" : "\t") + field;
    }
    return joined + "\n";
}

TEST(ProtonWindows, MalformedRowsAreParseErrorsNamingTheLine) {
    expect_parse_error(with_field(0, "no such group"), "line 1");
    expect_parse_error(with_field(2, "CH4"), "H count");
    expect_parse_error(with_field(2, "C"), "H count");  // no H, no 1H window
    expect_parse_error(with_field(12, "abc"), "line 1");
    expect_parse_error(with_field(12, "4.50"), "low");  // above the 4.31 high
    expect_parse_error(with_field(20, "2"), "round");
    expect_parse_error("methoxy\tO-CH3 carbon\tCH3\n", "21 fields");
    const std::string row = with_field(20, "1");
    expect_parse_error(row + row, "line 2");  // the same window twice
}

TEST(ProtonWindows, MeasuredWindowsAreTheDataFile) {
    const auto measured = arise::filter::measured_proton_windows();
    ASSERT_TRUE(measured.ok()) << measured.error().message;
    const auto file =
        parse_proton_windows(arise::test::join_chunks(arise::test::kGroupProtonShiftsChunks));
    ASSERT_TRUE(file.ok()) << file.error().message;
    EXPECT_EQ(measured.value(), file.value());
    EXPECT_GT(measured.value().size(), 20U);
}

TEST(ProtonWindows, RoundOneIsMethoxyAndBenzeneOnly) {
    // The chemist's decision, 2026-11-22.
    const auto measured = arise::filter::measured_proton_windows();
    ASSERT_TRUE(measured.ok());
    std::vector<std::pair<Group, CarbonType>> enabled;
    for (const auto& window : measured.value()) {
        if (window.round1) {
            enabled.emplace_back(window.group, window.type);
        }
    }
    EXPECT_EQ(enabled, (std::vector<std::pair<Group, CarbonType>>{
                           {Group::BenzeneRing, CarbonType::CH},
                           {Group::Methoxy, CarbonType::CH3},
                       }));
}

// ---- proton counts ----------------------------------------------------------------------

TEST(ProtonCount, RoundsTheIntegralHalvesUpAtLeastOne) {
    const auto count = [](std::optional<double> integral) {
        return proton_count({.shift = Ppm{3.8}, .integral = integral});
    };
    EXPECT_EQ(count(3.0), 3);
    EXPECT_EQ(count(2.6), 3);
    EXPECT_EQ(count(2.5), 3);
    EXPECT_EQ(count(2.4), 2);
    EXPECT_EQ(count(0.3), 1);
    EXPECT_EQ(count(std::nullopt), std::nullopt);
}

// ---- suggestions -------------------------------------------------------------------------

ProtonWindow window(Group group, CarbonType type, double low, double high, bool round1) {
    return {
        .group = group,
        .hallmark = "key carbon",
        .type = type,
        .low = Ppm{low},
        .high = Ppm{high},
        .n = 0,
        .molecules = 0,
        .precision = 0.7,
        .base_rate = 0.1,
        .round1 = round1,
    };
}

const std::vector<ProtonWindow>& windows() {
    static const std::vector<ProtonWindow> all{
        window(Group::Methoxy, CarbonType::CH3, 3.0, 4.3, true),
        window(Group::BenzeneRing, CarbonType::CH, 6.3, 8.9, true),
        window(Group::Ether, CarbonType::CH2, 2.7, 6.0, false),  // measured, not enabled
    };
    return all;
}

ObservedPeak peak(double shift, std::optional<double> integral) {
    return {.shift = Ppm{shift}, .integral = integral};
}

std::vector<Group> suggested(const std::vector<ObservedPeak>& peaks, ProtonOptions options = {}) {
    std::vector<Group> groups;
    for (const auto& s : suggest_groups_from_protons(peaks, windows(), options)) {
        groups.push_back(s.group);
    }
    return groups;
}

TEST(ProtonSuggestions, AThreeProtonPeakInTheWindowSuggestsMethoxy) {
    EXPECT_EQ(suggested({peak(3.8, 3.0)}), std::vector<Group>{Group::Methoxy});
}

TEST(ProtonSuggestions, TooFewProtonsDoNotSuggestAMethylWindow) {
    EXPECT_EQ(suggested({peak(3.8, 2.0)}), std::vector<Group>{});
    EXPECT_EQ(suggested({peak(3.8, 2.0)}, {.use_counts = false}),
              std::vector<Group>{Group::Methoxy});
    // More protons than the window's H count still match (overlapping peaks).
    EXPECT_EQ(suggested({peak(3.8, 5.0)}), std::vector<Group>{Group::Methoxy});
}

TEST(ProtonSuggestions, APeakWithoutAnIntegralMatchesOnShiftAlone) {
    EXPECT_EQ(suggested({peak(3.8, std::nullopt)}), std::vector<Group>{Group::Methoxy});
}

TEST(ProtonSuggestions, OnlyEnabledWindowsSuggest) {
    // 5.5 ppm, 2H: inside the ether CH2 window, which is not enabled for round 1.
    EXPECT_EQ(suggested({peak(5.5, 2.0)}), std::vector<Group>{});
}

TEST(ProtonSuggestions, BoundsAreInclusive) {
    EXPECT_EQ(suggested({peak(3.0, 3.0)}), std::vector<Group>{Group::Methoxy});
    EXPECT_EQ(suggested({peak(4.3, 3.0)}), std::vector<Group>{Group::Methoxy});
    EXPECT_EQ(suggested({peak(2.99, 3.0)}), std::vector<Group>{});
    EXPECT_EQ(suggested({peak(8.91, 1.0)}), std::vector<Group>{});
}

TEST(ProtonSuggestions, GroupsComeInEngineOrderWithEveryPeakAsEvidence) {
    const std::vector<ObservedPeak> peaks{
        peak(7.3, 5.0),
        peak(3.9, 3.0),
        peak(7.1, 1.0),
        peak(3.8, 3.0),
    };
    const auto suggestions = suggest_groups_from_protons(peaks, windows());
    ASSERT_EQ(suggestions.size(), 2U);
    // kAllGroups order: benzene ring (row 23) before methoxy (row 28).
    EXPECT_EQ(suggestions[0].group, Group::BenzeneRing);
    EXPECT_EQ(suggestions[1].group, Group::Methoxy);
    const std::vector<ProtonEvidence> benzene{
        {.peak = 0, .shift = Ppm{7.3}, .protons = 5, .window = windows()[1]},
        {.peak = 2, .shift = Ppm{7.1}, .protons = 1, .window = windows()[1]},
    };
    EXPECT_EQ(suggestions[0].evidence, benzene);
    EXPECT_EQ(suggestions[1].evidence.size(), 2U);
    EXPECT_EQ(suggestions[1].evidence[0].peak, 1U);
}

TEST(ProtonSuggestions, NoPeaksNoSuggestions) {
    EXPECT_EQ(suggest_groups_from_protons({}, windows()), std::vector<ProtonSuggestion>{});
}

}  // namespace

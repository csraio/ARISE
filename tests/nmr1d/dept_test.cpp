#include "arise/nmr1d/dept.hpp"

#include "arise/core/error.hpp"
#include "arise/core/units.hpp"
#include "arise/io/carbon_csv.hpp"

#include <gtest/gtest.h>

#include <limits>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::core::Ppm;
using arise::io::CarbonPeak;
using arise::io::DeptPeak;
using arise::nmr1d::CarbonAssignment;
using arise::nmr1d::CarbonType;
using arise::nmr1d::CarbonTypeSet;
using arise::nmr1d::classify_carbons;
using arise::nmr1d::DeptSpectra;

std::vector<CarbonPeak> carbons(std::initializer_list<double> shifts) {
    std::vector<CarbonPeak> peaks;
    peaks.reserve(shifts.size());
    for (const double shift : shifts) {
        peaks.push_back({.shift = Ppm{shift}, .intensity = std::nullopt});
    }
    return peaks;
}

DeptPeak dept(double shift, double intensity) {
    return {.shift = Ppm{shift}, .intensity = intensity};
}

std::vector<CarbonTypeSet> possible_types(const std::vector<CarbonAssignment>& assignments) {
    std::vector<CarbonTypeSet> sets;
    sets.reserve(assignments.size());
    for (const auto& assignment : assignments) {
        sets.push_back(assignment.possible);
    }
    return sets;
}

void expect_invalid(const DeptSpectra& spectra, const std::vector<std::string>& fragments) {
    const auto result = classify_carbons(spectra);
    ASSERT_FALSE(result.ok()) << "accepted inconsistent spectra";
    EXPECT_EQ(result.error().code, ErrorCode::InvalidArgument) << result.error().message;
    for (const auto& fragment : fragments) {
        EXPECT_NE(result.error().message.find(fragment), std::string::npos)
            << "message '" << result.error().message << "' lacks '" << fragment << "'";
    }
}

// Illustrative shifts, not measured data, for ethylbenzene's carbon types (one C, three
// aromatic CH environments, CH2, CH3). The classification depends only on the pattern.
std::vector<CarbonPeak> ethylbenzene() {
    return carbons({144.2, 128.4, 127.9, 125.7, 29.0, 15.6});
}
std::vector<DeptPeak> ethylbenzene_dept135() {
    return {
        dept(128.4, 1.0), dept(127.9, 0.9), dept(125.7, 0.5), dept(29.0, -0.8), dept(15.6, 0.7),
    };
}
std::vector<DeptPeak> ethylbenzene_dept90() {
    return {
        dept(128.4, 1.0),
        dept(127.9, 0.9),
        dept(125.7, 0.5),
    };
}

constexpr CarbonTypeSet kC{CarbonType::C};
constexpr CarbonTypeSet kCh{CarbonType::CH};
constexpr CarbonTypeSet kCh2{CarbonType::CH2};
constexpr CarbonTypeSet kCh3{CarbonType::CH3};

// ---- CarbonType and CarbonTypeSet --------------------------------------------------

static_assert(static_cast<int>(CarbonType::C) == 0);
static_assert(static_cast<int>(CarbonType::CH3) == 3);  // enumerator value = H count
static_assert(sizeof(CarbonTypeSet) == 1);

TEST(CarbonType, NamesFollowAttachedHydrogens) {
    EXPECT_EQ(arise::nmr1d::to_string(CarbonType::C), "C");
    EXPECT_EQ(arise::nmr1d::to_string(CarbonType::CH), "CH");
    EXPECT_EQ(arise::nmr1d::to_string(CarbonType::CH2), "CH2");
    EXPECT_EQ(arise::nmr1d::to_string(CarbonType::CH3), "CH3");
}

TEST(CarbonTypeSet, BasicSetOperations) {
    constexpr CarbonTypeSet kEmpty;
    static_assert(kEmpty.empty());
    static_assert(CarbonTypeSet::all().size() == 4);
    constexpr CarbonTypeSet kChOrCh3{CarbonType::CH, CarbonType::CH3};
    static_assert(kChOrCh3.contains(CarbonType::CH3) && !kChOrCh3.contains(CarbonType::CH2));
    static_assert((kChOrCh3 & CarbonTypeSet{CarbonType::CH}) == CarbonTypeSet{CarbonType::CH});
    static_assert((kChOrCh3 & CarbonTypeSet{CarbonType::CH2}).empty());
    static_assert((CarbonTypeSet{CarbonType::CH} | CarbonTypeSet{CarbonType::CH3}) == kChOrCh3);
    static_assert((kChOrCh3 | kEmpty) == kChOrCh3);
    static_assert((kChOrCh3 | CarbonTypeSet{CarbonType::C, CarbonType::CH2}) ==
                  CarbonTypeSet::all());
    EXPECT_EQ(kChOrCh3.size(), 2);
}

TEST(CarbonTypeSet, ResolvedOnlyWhenExactlyOneTypeRemains) {
    EXPECT_EQ(CarbonTypeSet{CarbonType::CH2}.resolved(), CarbonType::CH2);
    EXPECT_EQ((CarbonTypeSet{CarbonType::CH, CarbonType::CH3}.resolved()), std::nullopt);
    EXPECT_EQ(CarbonTypeSet{}.resolved(), std::nullopt);
}

TEST(CarbonTypeSet, PrintsInHydrogenOrder) {
    EXPECT_EQ(arise::nmr1d::to_string(CarbonTypeSet{CarbonType::CH3, CarbonType::CH}), "{CH, CH3}");
    EXPECT_EQ(arise::nmr1d::to_string(CarbonTypeSet::all()), "{C, CH, CH2, CH3}");
    EXPECT_EQ(arise::nmr1d::to_string(CarbonTypeSet{}), "{}");
}

// ---- Classification -----------------------------------------------------------------

TEST(Dept, AllThreeSpectraResolveEveryCarbon) {
    const auto carbon = ethylbenzene();
    const auto dept90 = ethylbenzene_dept90();
    const auto dept135 = ethylbenzene_dept135();
    const auto result = classify_carbons({.carbon = carbon, .dept135 = dept135, .dept90 = dept90});
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<CarbonTypeSet> expected{kC, kCh, kCh, kCh, kCh2, kCh3};
    EXPECT_EQ(possible_types(result.value()), expected);
    EXPECT_EQ(result.value()[0].shift, Ppm{144.2});
}

TEST(Dept, CarbonOnlyLeavesEveryTypePossible) {
    const auto carbon = ethylbenzene();
    const auto result = classify_carbons({.carbon = carbon});
    ASSERT_TRUE(result.ok()) << result.error().message;
    for (const auto& assignment : result.value()) {
        EXPECT_EQ(assignment.possible, CarbonTypeSet::all());
    }
}

TEST(Dept, Dept135AloneCannotSeparateMethineFromMethyl) {
    const auto carbon = ethylbenzene();
    const auto dept135 = ethylbenzene_dept135();
    const auto result = classify_carbons({.carbon = carbon, .dept135 = dept135});
    ASSERT_TRUE(result.ok()) << result.error().message;
    const CarbonTypeSet ch_or_ch3{CarbonType::CH, CarbonType::CH3};
    const std::vector<CarbonTypeSet> expected{kC, ch_or_ch3, ch_or_ch3, ch_or_ch3, kCh2, ch_or_ch3};
    EXPECT_EQ(possible_types(result.value()), expected);
}

TEST(Dept, Dept90AloneFindsOnlyMethines) {
    const auto carbon = ethylbenzene();
    const auto dept90 = ethylbenzene_dept90();
    const auto result = classify_carbons({.carbon = carbon, .dept90 = dept90});
    ASSERT_TRUE(result.ok()) << result.error().message;
    const CarbonTypeSet not_ch{CarbonType::C, CarbonType::CH2, CarbonType::CH3};
    const std::vector<CarbonTypeSet> expected{not_ch, kCh, kCh, kCh, not_ch, not_ch};
    EXPECT_EQ(possible_types(result.value()), expected);
}

TEST(Dept, EmptyDept90MeansNoMethines) {
    // Illustrative shifts, not measured data, for 2-butanone's carbons (C=O, CH2, 2 CH3).
    const auto butanone = carbons({209.3, 36.9, 29.4, 7.9});
    const std::vector<DeptPeak> dept135{dept(36.9, -1.0), dept(29.4, 1.0), dept(7.9, 0.9)};
    const auto result = classify_carbons(
        {.carbon = butanone, .dept135 = dept135, .dept90 = std::span<const DeptPeak>{}});
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<CarbonTypeSet> expected{kC, kCh2, kCh3, kCh3};
    EXPECT_EQ(possible_types(result.value()), expected);
}

TEST(Dept, KeepsThe13COrderWhateverTheDeptOrder) {
    const auto carbon = ethylbenzene();
    const std::vector<DeptPeak> shuffled{
        dept(15.6, 0.7), dept(29.0, -0.8), dept(125.7, 0.5), dept(128.4, 1.0), dept(127.9, 0.9),
    };
    const auto result = classify_carbons({.carbon = carbon, .dept135 = shuffled});
    ASSERT_TRUE(result.ok()) << result.error().message;
    ASSERT_EQ(result.value().size(), carbon.size());
    for (std::size_t i = 0; i < carbon.size(); ++i) {
        EXPECT_EQ(result.value()[i].shift, carbon[i].shift);
    }
    EXPECT_EQ(result.value()[4].possible, kCh2);
}

TEST(Dept, MatchesWithinToleranceInclusive) {
    // 29.0625 - 29.0 = 0.0625 exactly (a power of two), so "inclusive" is testable.
    const auto carbon = carbons({29.0});
    const std::vector<DeptPeak> at_edge{dept(29.0625, -1.0)};
    const auto result =
        classify_carbons({.carbon = carbon, .dept135 = at_edge, .tolerance = Ppm{0.0625}});
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_EQ(result.value()[0].possible, kCh2);

    const std::vector<DeptPeak> beyond{dept(29.125, -1.0)};
    expect_invalid({.carbon = carbon, .dept135 = beyond, .tolerance = Ppm{0.0625}},
                   {"DEPT-135", "29.125", "no 13C peak"});
}

TEST(Dept, DeptPeakWithoutCarbonPartnerIsError) {
    const auto carbon = ethylbenzene();
    const std::vector<DeptPeak> extra{dept(128.4, 1.0), dept(50.1, 1.0)};
    expect_invalid({.carbon = carbon, .dept135 = extra}, {"DEPT-135", "50.1"});
    expect_invalid({.carbon = carbon, .dept90 = extra}, {"DEPT-90", "50.1"});
}

TEST(Dept, DeptPeakMatchingTwoCarbonPeaksIsAmbiguous) {
    const auto close = carbons({128.40, 128.43});
    const std::vector<DeptPeak> between{dept(128.42, 1.0)};
    expect_invalid({.carbon = close, .dept135 = between}, {"128.42", "ambiguous", "tolerance"});
}

TEST(Dept, TwoDeptPeaksOnOneCarbonPeakIsError) {
    const auto carbon = carbons({128.40});
    const std::vector<DeptPeak> two{dept(128.39, 1.0), dept(128.41, 1.0)};
    expect_invalid({.carbon = carbon, .dept135 = two}, {"DEPT-135", "128.4"});
}

TEST(Dept, NegativeDept90PeakIsError) {
    const auto carbon = ethylbenzene();
    // DEPT-90 shows only CH, and CH is positive: a negative peak means wrong phasing.
    const std::vector<DeptPeak> inverted{dept(128.4, -1.0)};
    expect_invalid({.carbon = carbon, .dept90 = inverted}, {"DEPT-90", "128.4", "negative"});
}

TEST(Dept, ContradictingSpectraAreErrors) {
    // DEPT-90 says CH, DEPT-135 says CH2 (negative): impossible together.
    const auto carbon = carbons({29.0});
    const std::vector<DeptPeak> ch2{dept(29.0, -1.0)};
    const std::vector<DeptPeak> ch{dept(29.0, 1.0)};
    expect_invalid({.carbon = carbon, .dept135 = ch2, .dept90 = ch}, {"29", "DEPT-90", "DEPT-135"});
    // DEPT-90 says CH, DEPT-135 shows nothing (C): also impossible.
    expect_invalid({.carbon = carbon, .dept135 = std::span<const DeptPeak>{}, .dept90 = ch},
                   {"29", "DEPT-90", "DEPT-135"});
}

TEST(Dept, EmptyCarbonListIsError) {
    expect_invalid({.carbon = std::span<const CarbonPeak>{}}, {"13C"});
}

TEST(Dept, InvalidToleranceIsError) {
    const auto carbon = ethylbenzene();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const double tolerance : {0.0, -0.05, nan}) {
        expect_invalid({.carbon = carbon, .tolerance = Ppm{tolerance}}, {"tolerance"});
    }
}

}  // namespace

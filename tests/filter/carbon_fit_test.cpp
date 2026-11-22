// filter::check_carbon_fit and order_by_fit (plan 3c): engineering tests of the rules in
// carbon_fit.hpp. Each molecule's environment counts come from chem::nmr_equivalence
// (tested in tests/chem against nmrshiftdb2); the shifts here are illustrative, not data.

#include "arise/filter/carbon_fit.hpp"

#include "arise/chem/molecule.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <initializer_list>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using arise::core::Ppm;
using arise::filter::CarbonFit;
using arise::filter::check_carbon_fit;
using arise::filter::order_by_fit;
using arise::filter::SuggestOptions;
using arise::filter::Tier;
using arise::nmr1d::CarbonAssignment;
using arise::nmr1d::CarbonType;
using arise::nmr1d::CarbonTypeSet;

arise::chem::Molecule molecule(std::string_view smiles) {
    auto parsed = arise::chem::parse_smiles(smiles);
    EXPECT_TRUE(parsed) << smiles;
    return parsed ? std::move(parsed).value() : arise::chem::Molecule{};
}

CarbonAssignment peak(double shift, CarbonTypeSet possible) {
    return {.shift = Ppm{shift}, .possible = possible};
}

CarbonFit fit(std::string_view smiles, const std::vector<CarbonAssignment>& peaks,
              SuggestOptions options = {}) {
    const auto result = check_carbon_fit(molecule(smiles), peaks, options);
    EXPECT_TRUE(result) << smiles << ": " << (result ? "" : result.error().message);
    return result ? result.value() : CarbonFit{};
}

constexpr CarbonTypeSet kC{CarbonType::C};
constexpr CarbonTypeSet kCH2{CarbonType::CH2};
constexpr CarbonTypeSet kCH3{CarbonType::CH3};
constexpr auto kAny = CarbonTypeSet::all();

// Ethyl acetate, CCOC(=O)C: four environments always (CH3, CH2, C=O, CH3).
constexpr std::string_view kEthylAcetate = "CCOC(=O)C";

std::vector<CarbonAssignment> ethyl_acetate_peaks() {
    return {peak(14.0, kCH3), peak(60.0, kCH2), peak(171.0, kC), peak(21.0, kCH3)};
}

TEST(CarbonFit, AnExactFitIsActiveWithNoOverlaps) {
    const auto result = fit(kEthylAcetate, ethyl_acetate_peaks());
    EXPECT_EQ(result.verdict.tier, Tier::Active);
    EXPECT_TRUE(result.verdict.reasons.empty());
    EXPECT_EQ(result.assumed_overlaps, 0);
    EXPECT_EQ(result.environments_fewest, 4);
    EXPECT_EQ(result.environments_most, 4);
    EXPECT_TRUE(result.notes.empty());
}

TEST(CarbonFit, FewerPeaksStayActiveAndCountAssumedOverlaps) {
    auto peaks = ethyl_acetate_peaks();
    peaks.pop_back();  // the two CH3 peaks taken as one
    const auto result = fit(kEthylAcetate, peaks);
    EXPECT_EQ(result.verdict.tier, Tier::Active);
    EXPECT_EQ(result.assumed_overlaps, 1);
    EXPECT_EQ(result.notes, std::vector<std::string>{"assumes 1 overlapping peak"});
}

TEST(CarbonFit, SeveralOverlapsAreCountedAndStillOnlyNoted) {
    const auto result = fit(kEthylAcetate, {peak(60.0, kAny), peak(171.0, kAny)});
    EXPECT_EQ(result.verdict.tier, Tier::Active);
    EXPECT_EQ(result.assumed_overlaps, 2);
    EXPECT_EQ(result.notes, std::vector<std::string>{"assumes 2 overlapping peaks"});
}

TEST(CarbonFit, MorePeaksThanEnvironmentsGoToTheBackBurner) {
    auto peaks = ethyl_acetate_peaks();
    peaks.push_back(peak(30.0, kAny));
    const auto result = fit(kEthylAcetate, peaks);
    EXPECT_EQ(result.verdict.tier, Tier::BackBurner);
    EXPECT_EQ(result.verdict.reasons,
              std::vector<std::string>{
                  "13C list has 5 peaks; the candidate has at most 4 carbon environments"});
}

TEST(CarbonFit, APeakWithNoEnvironmentOfItsTypeLeftGoesToTheBackBurner) {
    // Two CH2 peaks, but ethyl acetate has one CH2 environment.
    const auto result =
        fit(kEthylAcetate, {peak(14.0, kCH3), peak(60.0, kCH2), peak(45.5, kCH2), peak(171.0, kC)});
    EXPECT_EQ(result.verdict.tier, Tier::BackBurner);
    ASSERT_EQ(result.verdict.reasons.size(), 1U);
    const auto& reason = result.verdict.reasons.front();
    EXPECT_NE(reason.find("{CH2}"), std::string::npos) << reason;
    EXPECT_TRUE(reason.find("45.5 ppm") != std::string::npos ||
                reason.find("60 ppm") != std::string::npos)
        << reason;
}

TEST(CarbonFit, TheMatchingFindsAFitWhereAGreedyChoiceWouldNot) {
    // The first peak could be CH3 or CH2; taking the CH2 environment for it would leave the
    // second (CH2 only) without one. A maximum matching gives it the CH3.
    const CarbonTypeSet ch2_or_ch3{CarbonType::CH2, CarbonType::CH3};
    const auto result = fit(kEthylAcetate, {
                                               peak(20.0, ch2_or_ch3),
                                               peak(60.0, kCH2),
                                               peak(171.0, kC),
                                               peak(14.0, kCH3),
                                           });
    EXPECT_EQ(result.verdict.tier, Tier::Active);
}

TEST(CarbonFit, WithDeptOffOnlyTheCountIsCompared) {
    const SuggestOptions shift_only{.use_dept = false};
    const auto result =
        fit(kEthylAcetate, {peak(14.0, kCH3), peak(60.0, kCH2), peak(45.5, kCH2), peak(171.0, kC)},
            shift_only);
    EXPECT_EQ(result.verdict.tier, Tier::Active);
}

TEST(CarbonFit, PeaksOfUnknownTypeMatchAnyEnvironment) {
    const auto result = fit(
        kEthylAcetate, {peak(14.0, kAny), peak(60.0, kAny), peak(171.0, kAny), peak(21.0, kAny)});
    EXPECT_EQ(result.verdict.tier, Tier::Active);
}

TEST(CarbonFit, ARangeCandidateFitsAnyCountInsideItsRange) {
    // N-acetylpyrrole, CC(=O)N1C=CC=C1: 4 environments always, up to 6 (the ring halves
    // may differ: the chemist's rule, 2026-11-08).
    const auto four = fit("CC(=O)N1C=CC=C1", {
                                                 peak(22.0, kAny),
                                                 peak(168.0, kAny),
                                                 peak(119.0, kAny),
                                                 peak(113.0, kAny),
                                             });
    EXPECT_EQ(four.environments_fewest, 4);
    EXPECT_EQ(four.environments_most, 6);
    EXPECT_EQ(four.verdict.tier, Tier::Active);
    EXPECT_EQ(four.assumed_overlaps, 0);
    std::vector<CarbonAssignment> six;
    for (const double shift : {22.0, 168.0, 119.0, 118.0, 113.0, 112.0}) {
        six.push_back(peak(shift, kAny));
    }
    EXPECT_EQ(fit("CC(=O)N1C=CC=C1", six).verdict.tier, Tier::Active);
    six.push_back(peak(130.0, kAny));
    EXPECT_EQ(fit("CC(=O)N1C=CC=C1", six).verdict.tier, Tier::BackBurner);
}

TEST(CarbonFit, PeakOrderDoesNotChangeTheResult) {
    auto peaks = ethyl_acetate_peaks();
    const auto forward = fit(kEthylAcetate, peaks);
    std::reverse(peaks.begin(), peaks.end());
    EXPECT_EQ(fit(kEthylAcetate, peaks), forward);
}

TEST(CarbonFit, AnEmptyPeakListIsAnError) {
    const auto result = check_carbon_fit(molecule(kEthylAcetate), {});
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, arise::core::ErrorCode::InvalidArgument);
}

// ---- order_by_fit ------------------------------------------------------------------------

CarbonFit ranked(Tier tier, int overlaps, std::size_t reasons) {
    CarbonFit result;
    result.verdict.tier = tier;
    result.verdict.reasons.assign(reasons, "reason");
    result.assumed_overlaps = overlaps;
    return result;
}

TEST(CarbonFitOrder, ActiveFirstThenFewerOverlapsThenFewerReasonsThenInputOrder) {
    const std::vector<CarbonFit> fits{
        ranked(Tier::BackBurner, 0, 1),  // 0
        ranked(Tier::Active, 2, 0),      // 1
        ranked(Tier::Active, 0, 0),      // 2
        ranked(Tier::BackBurner, 0, 2),  // 3
        ranked(Tier::Active, 0, 0),      // 4: ties with 2, keeps input order
        ranked(Tier::Eliminated, 0, 1),  // 5
        ranked(Tier::BackBurner, 1, 1),  // 6
    };
    EXPECT_EQ(order_by_fit(fits), (std::vector<std::size_t>{2, 4, 1, 0, 3, 6, 5}));
}

TEST(CarbonFitOrder, FewerReasonsComeFirstEvenWhenListedLater) {
    const std::vector<CarbonFit> fits{
        ranked(Tier::BackBurner, 0, 2),
        ranked(Tier::BackBurner, 0, 1),
    };
    EXPECT_EQ(order_by_fit(fits), (std::vector<std::size_t>{1, 0}));
}

TEST(CarbonFitOrder, AnEmptyListGivesAnEmptyOrder) {
    EXPECT_TRUE(order_by_fit({}).empty());
}

}  // namespace

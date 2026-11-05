// Engineering tests of hsqc::classify_protons and suggest_exchangeable_count: the rules as
// documented in exchangeable.hpp, their order, the tolerance bound and every error. No
// chemistry claim rests on these values: solvent and water shifts are read from
// nmr1d::solvent_peaks, and the other shifts are placed well away from them. The
// chemistry check is bmrb_test.cpp.

#include "arise/hsqc/exchangeable.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::core::Ppm;
using arise::hsqc::classify_protons;
using arise::hsqc::ExchangeableCount;
using arise::hsqc::kDefaultProtonTolerance;
using arise::hsqc::ProtonClass;
using arise::hsqc::ProtonClassification;
using arise::hsqc::ProtonSpectra;
using arise::hsqc::ProtonVerdict;
using arise::hsqc::suggest_exchangeable_count;
using arise::io::CorrelationPeak;
using arise::io::ObservedPeak;
using arise::nmr1d::Solvent;

constexpr Solvent kSolvent = Solvent::CDCl3;

const arise::nmr1d::SolventPeaks& signals() {
    static const auto peaks = arise::nmr1d::solvent_peaks(kSolvent).value();
    return peaks;
}

Ppm residual() {
    return signals().residual.front();
}
Ppm water() {
    return signals().water.front();
}

// A shift at least 1 ppm from every solvent and water signal of kSolvent.
Ppm clear_shift(double near) {
    double shift = near;
    const auto close = [&shift](Ppm signal) { return std::abs(signal.value() - shift) < 1.0; };
    while (std::any_of(signals().residual.begin(), signals().residual.end(), close) ||
           std::any_of(signals().water.begin(), signals().water.end(), close)) {
        shift += 1.0;
    }
    return Ppm{shift};
}

ObservedPeak proton(Ppm shift, std::optional<double> integral = std::nullopt) {
    return {.shift = shift, .integral = integral};
}

CorrelationPeak cross(Ppm proton_shift) {
    return {.proton = proton_shift, .carbon = Ppm{30.0}, .intensity = std::nullopt};
}

ProtonClassification classify(const std::vector<ObservedPeak>& protons,
                              const std::vector<CorrelationPeak>& hsqc,
                              const std::optional<std::vector<ObservedPeak>>& after_d2o = {},
                              Ppm tolerance = kDefaultProtonTolerance) {
    const auto result = classify_protons({
        .proton = protons,
        .hsqc = hsqc,
        .solvent = kSolvent,
        .after_d2o =
            after_d2o ? std::optional{std::span<const ObservedPeak>{*after_d2o}} : std::nullopt,
        .tolerance = tolerance,
    });
    EXPECT_TRUE(result) << (result ? "" : result.error().message);
    return result ? result.value() : ProtonClassification{};
}

ProtonVerdict only_verdict(const ProtonClassification& classification) {
    EXPECT_EQ(classification.peaks.size(), 1U);
    return classification.peaks.empty() ? ProtonVerdict::NoHsqcPartner
                                        : classification.peaks.front().verdict;
}

// ---- verdicts ----------------------------------------------------------------------------

TEST(HsqcClassifyProtons, PeakWithAnHsqcPartnerIsOnCarbon) {
    const auto at = clear_shift(3.0);
    const auto result = classify({proton(at)}, {cross(at)});
    EXPECT_EQ(only_verdict(result), ProtonVerdict::OnCarbon);
    EXPECT_EQ(result.peaks.front().hsqc_partners, std::vector{cross(at)});
    EXPECT_TRUE(result.peaks.front().solvent_signals.empty());
}

TEST(HsqcClassifyProtons, PeakWithNeitherPartnerNorSolventIsNoHsqcPartner) {
    const auto at = clear_shift(3.0);
    const auto result = classify({proton(at)}, {cross(Ppm{at.value() + 1.0})});
    EXPECT_EQ(only_verdict(result), ProtonVerdict::NoHsqcPartner);
    EXPECT_TRUE(result.peaks.front().hsqc_partners.empty());
    EXPECT_FALSE(result.peaks.front().gone_after_d2o);
}

TEST(HsqcClassifyProtons, PeakAtTheResidualSolventSignalIsSolventOrWater) {
    const auto result = classify({proton(residual())}, {cross(clear_shift(3.0))});
    EXPECT_EQ(only_verdict(result), ProtonVerdict::SolventOrWater);
    const std::vector<arise::hsqc::SolventSignal> expected{{.shift = residual(), .water = false}};
    EXPECT_EQ(result.peaks.front().solvent_signals, expected);
}

TEST(HsqcClassifyProtons, PeakAtTheWaterSignalIsSolventOrWaterMarkedAsWater) {
    const auto result = classify({proton(water())}, {cross(clear_shift(3.0))});
    EXPECT_EQ(only_verdict(result), ProtonVerdict::SolventOrWater);
    const std::vector<arise::hsqc::SolventSignal> expected{{.shift = water(), .water = true}};
    EXPECT_EQ(result.peaks.front().solvent_signals, expected);
}

TEST(HsqcClassifyProtons, PeakAtTheSolventWithAnHsqcPartnerIsUndetermined) {
    const auto result = classify({proton(residual())}, {cross(residual())});
    EXPECT_EQ(only_verdict(result), ProtonVerdict::Undetermined);
    EXPECT_EQ(result.peaks.front().hsqc_partners.size(), 1U);
    EXPECT_EQ(result.peaks.front().solvent_signals.size(), 1U);
}

TEST(HsqcClassifyProtons, PeakGoneAfterTheD2oShakeIsExchangeableWhateverElseMatches) {
    // Rule 1 overrides an HSQC partner and a solvent signal alike.
    const auto at = clear_shift(3.0);
    const auto other = clear_shift(at.value() + 1.0);
    const auto result = classify({proton(at), proton(residual())}, {cross(at), cross(residual())},
                                 std::vector{proton(other)});
    ASSERT_EQ(result.peaks.size(), 2U);
    for (const auto& peak : result.peaks) {
        EXPECT_EQ(peak.verdict, ProtonVerdict::Exchangeable);
        EXPECT_TRUE(peak.gone_after_d2o);
    }
}

TEST(HsqcClassifyProtons, PeakStillThereAfterTheD2oShakeKeepsItsOtherVerdict) {
    const auto at = clear_shift(3.0);
    const auto lone = clear_shift(at.value() + 1.0);
    const auto result =
        classify({proton(at), proton(lone)}, {cross(at)}, std::vector{proton(at), proton(lone)});
    ASSERT_EQ(result.peaks.size(), 2U);
    EXPECT_EQ(result.peaks[0].verdict, ProtonVerdict::OnCarbon);
    EXPECT_EQ(result.peaks[1].verdict, ProtonVerdict::NoHsqcPartner);
    EXPECT_FALSE(result.peaks[0].gone_after_d2o);
    EXPECT_FALSE(result.peaks[1].gone_after_d2o);
}

TEST(HsqcClassifyProtons, NewPeaksAfterTheD2oShakeAreReportedNotUsed) {
    // The shake makes a new HDO peak; it matches no original peak.
    const auto at = clear_shift(3.0);
    const auto fresh = proton(clear_shift(at.value() + 2.0));
    const auto result = classify({proton(at)}, {cross(at)}, std::vector{proton(at), fresh});
    EXPECT_EQ(only_verdict(result), ProtonVerdict::OnCarbon);
    EXPECT_EQ(result.new_after_d2o, std::vector{fresh});
}

TEST(HsqcClassifyProtons, WithoutAD2oListNothingIsNewOrGone) {
    const auto at = clear_shift(3.0);
    const auto result = classify({proton(at)}, {cross(at)});
    EXPECT_TRUE(result.new_after_d2o.empty());
    EXPECT_FALSE(result.peaks.front().gone_after_d2o);
}

// ---- tolerance ---------------------------------------------------------------------------

TEST(HsqcClassifyProtons, DefaultToleranceIsTheChemistsChoice) {
    // 0.15 ppm, chosen by the chemist on 2026-11-05 from the BMRB measurement.
    EXPECT_EQ(kDefaultProtonTolerance, Ppm{0.15});
}

TEST(HsqcClassifyProtons, ToleranceIsInclusive) {
    // 0.25 and 0.5 are exact in binary, so the bound is tested without rounding.
    const auto at = clear_shift(3.0);
    const Ppm tolerance{0.25};
    const auto edge = classify({proton(at)}, {cross(Ppm{at.value() + 0.25})}, {}, tolerance);
    EXPECT_EQ(only_verdict(edge), ProtonVerdict::OnCarbon);
    const auto beyond = classify({proton(at)}, {cross(Ppm{at.value() + 0.5})}, {}, tolerance);
    EXPECT_EQ(only_verdict(beyond), ProtonVerdict::NoHsqcPartner);
}

TEST(HsqcClassifyProtons, ToleranceAppliesToSolventSignalsAndTheD2oList) {
    const Ppm tolerance{0.25};
    const auto near_residual = proton(Ppm{residual().value() + 0.25});
    const auto result = classify({near_residual}, {cross(clear_shift(3.0))}, {}, tolerance);
    EXPECT_EQ(only_verdict(result), ProtonVerdict::SolventOrWater);

    const auto at = clear_shift(3.0);
    const auto kept =
        classify({proton(at)}, {cross(at)}, std::vector{proton(Ppm{at.value() + 0.25})}, tolerance);
    EXPECT_EQ(only_verdict(kept), ProtonVerdict::OnCarbon);
    const auto gone =
        classify({proton(at)}, {cross(at)}, std::vector{proton(Ppm{at.value() + 0.5})}, tolerance);
    EXPECT_EQ(only_verdict(gone), ProtonVerdict::Exchangeable);
}

// ---- evidence and order ------------------------------------------------------------------

TEST(HsqcClassifyProtons, EveryPartnerWithinToleranceIsReportedInHsqcOrder) {
    const auto at = clear_shift(3.0);
    const auto above = cross(Ppm{at.value() + 0.1});
    const auto below = cross(Ppm{at.value() - 0.1});
    const auto far = cross(Ppm{at.value() + 1.0});
    const auto result = classify({proton(at)}, {above, far, below});
    EXPECT_EQ(result.peaks.front().hsqc_partners, (std::vector{above, below}));
}

TEST(HsqcClassifyProtons, VerdictsFollowTheProtonInputOrder) {
    const auto a = clear_shift(3.0);
    const auto b = clear_shift(a.value() + 1.0);
    const auto forward = classify({proton(a), proton(b), proton(residual())}, {cross(a)});
    const auto backward = classify({proton(residual()), proton(b), proton(a)}, {cross(a)});
    ASSERT_EQ(forward.peaks.size(), 3U);
    ASSERT_EQ(backward.peaks.size(), 3U);
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_EQ(forward.peaks[i], backward.peaks[2 - i]);
    }
    EXPECT_EQ(forward.peaks[0].peak, proton(a));
}

TEST(HsqcClassifyProtons, HsqcOrderChangesNoVerdict) {
    const auto a = clear_shift(3.0);
    const auto b = clear_shift(a.value() + 1.0);
    const auto one = classify({proton(a), proton(b)}, {cross(a), cross(Ppm{b.value() + 5.0})});
    const auto two = classify({proton(a), proton(b)}, {cross(Ppm{b.value() + 5.0}), cross(a)});
    ASSERT_EQ(one.peaks.size(), two.peaks.size());
    for (std::size_t i = 0; i < one.peaks.size(); ++i) {
        EXPECT_EQ(one.peaks[i].verdict, two.peaks[i].verdict);
    }
}

TEST(HsqcClassifyProtons, VerdictNamesAreDistinctAndNonEmpty) {
    const std::vector verdicts{
        ProtonVerdict::OnCarbon,      ProtonVerdict::SolventOrWater, ProtonVerdict::Exchangeable,
        ProtonVerdict::NoHsqcPartner, ProtonVerdict::Undetermined,
    };
    for (std::size_t i = 0; i < verdicts.size(); ++i) {
        EXPECT_FALSE(arise::hsqc::to_string(verdicts[i]).empty());
        for (std::size_t j = i + 1; j < verdicts.size(); ++j) {
            EXPECT_NE(arise::hsqc::to_string(verdicts[i]), arise::hsqc::to_string(verdicts[j]));
        }
    }
}

// ---- errors ------------------------------------------------------------------------------

ErrorCode error_of(const std::vector<ObservedPeak>& protons,
                   const std::vector<CorrelationPeak>& hsqc,
                   const std::optional<std::vector<ObservedPeak>>& after_d2o = {},
                   Ppm tolerance = kDefaultProtonTolerance) {
    const auto result = classify_protons({
        .proton = protons,
        .hsqc = hsqc,
        .solvent = kSolvent,
        .after_d2o =
            after_d2o ? std::optional{std::span<const ObservedPeak>{*after_d2o}} : std::nullopt,
        .tolerance = tolerance,
    });
    EXPECT_FALSE(result);
    return result ? ErrorCode::NotFound : result.error().code;
}

TEST(HsqcClassifyProtons, EmptyListsAreRejected) {
    const auto at = clear_shift(3.0);
    EXPECT_EQ(error_of({}, {cross(at)}), ErrorCode::InvalidArgument);
    EXPECT_EQ(error_of({proton(at)}, {}), ErrorCode::InvalidArgument);
}

TEST(HsqcClassifyProtons, NonFiniteShiftsAreRejected) {
    const auto at = clear_shift(3.0);
    const Ppm nan{std::numeric_limits<double>::quiet_NaN()};
    const Ppm inf{std::numeric_limits<double>::infinity()};
    EXPECT_EQ(error_of({proton(nan)}, {cross(at)}), ErrorCode::InvalidArgument);
    EXPECT_EQ(error_of({proton(at)}, {cross(inf)}), ErrorCode::InvalidArgument);
    EXPECT_EQ(error_of({proton(at)}, {cross(at)}, std::vector{proton(nan)}),
              ErrorCode::InvalidArgument);
}

TEST(HsqcClassifyProtons, NonPositiveOrNonFiniteToleranceIsRejected) {
    const auto at = clear_shift(3.0);
    for (const double tolerance : {
             0.0,
             -0.1,
             std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::infinity(),
         }) {
        EXPECT_EQ(error_of({proton(at)}, {cross(at)}, {}, Ppm{tolerance}),
                  ErrorCode::InvalidArgument)
            << tolerance;
    }
}

// ---- suggested count ---------------------------------------------------------------------

ProtonClass verdict(ProtonVerdict kind, std::optional<double> integral = std::nullopt) {
    return {
        .peak = proton(Ppm{3.0}, integral),
        .verdict = kind,
        .hsqc_partners = {},
        .solvent_signals = {},
        .gone_after_d2o = kind == ProtonVerdict::Exchangeable,
    };
}

ExchangeableCount count(const std::vector<ProtonClass>& classes, Solvent solvent = kSolvent) {
    const auto result = suggest_exchangeable_count(classes, solvent);
    EXPECT_TRUE(result) << (result ? "" : result.error().message);
    return result ? result.value() : ExchangeableCount{};
}

TEST(HsqcExchangeableCount, NoPeaksSuggestNone) {
    EXPECT_EQ(count({}), (ExchangeableCount{.min = 0, .max = 0}));
}

TEST(HsqcExchangeableCount, EachUnpartneredOrExchangeablePeakCountsOneWithoutAnIntegral) {
    const auto found = count({
        verdict(ProtonVerdict::NoHsqcPartner),
        verdict(ProtonVerdict::Exchangeable),
        verdict(ProtonVerdict::NoHsqcPartner),
    });
    EXPECT_EQ(found, (ExchangeableCount{.min = 3, .max = 3}));
}

TEST(HsqcExchangeableCount, OnCarbonAndSolventPeaksDoNotCount) {
    const auto found = count({
        verdict(ProtonVerdict::OnCarbon, 3.0),
        verdict(ProtonVerdict::SolventOrWater, 2.0),
        verdict(ProtonVerdict::NoHsqcPartner),
    });
    EXPECT_EQ(found, (ExchangeableCount{.min = 1, .max = 1}));
}

TEST(HsqcExchangeableCount, UndeterminedPeaksRaiseOnlyTheMaximum) {
    const auto found =
        count({verdict(ProtonVerdict::NoHsqcPartner), verdict(ProtonVerdict::Undetermined, 2.0)});
    EXPECT_EQ(found, (ExchangeableCount{.min = 1, .max = 3}));
}

TEST(HsqcExchangeableCount, IntegralsRoundToTheNearestWholeProtonAndAtLeastOne) {
    EXPECT_EQ(count({verdict(ProtonVerdict::NoHsqcPartner, 1.6)}).min, 2);
    EXPECT_EQ(count({verdict(ProtonVerdict::NoHsqcPartner, 2.4)}).min, 2);
    EXPECT_EQ(count({verdict(ProtonVerdict::NoHsqcPartner, 0.3)}).min, 1);
    EXPECT_EQ(count({verdict(ProtonVerdict::Exchangeable, 2.5)}).min, 3);
}

bool exchanges(Solvent solvent) {
    return arise::nmr1d::solvent_peaks(solvent).value().exchanges_labile_h;
}

TEST(HsqcExchangeableCount, ExchangingSolventsLeaveTheMaximumOpen) {
    // The solvents that swap labile H for D, as nmr1d's table marks them.
    int exchanging = 0;
    for (const auto solvent : arise::nmr1d::kAllSolvents) {
        const auto found = count({verdict(ProtonVerdict::NoHsqcPartner)}, solvent);
        const std::optional<int> open = exchanges(solvent) ? std::nullopt : std::optional{1};
        exchanging += exchanges(solvent) ? 1 : 0;
        EXPECT_EQ(found, (ExchangeableCount{.min = 1, .max = open}))
            << arise::nmr1d::to_string(solvent);
    }
    EXPECT_GT(exchanging, 0);
}

}  // namespace

#pragma once

#include "arise/core/result.hpp"
#include "arise/core/units.hpp"
#include "arise/io/correlation_csv.hpp"
#include "arise/io/peak_csv.hpp"
#include "arise/nmr1d/solvent.hpp"

#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace arise::hsqc {

/// How far a 1H peak's shift may lie from an HSQC cross peak's 1H shift, or from a solvent
/// or water signal, and still match it (inclusive). A convention the chemist chose
/// (2026-11-05), not a constant of nature. On BMRB's experimental spectra
/// (tests/data/reference/bmrb_*.tsv), 98.4% of 1H peaks on carbon have an HSQC partner
/// within 0.15 ppm, against 94.3% within 0.05 ppm (tests/hsqc/bmrb_test.cpp). The cost: an
/// exchangeable peak within this distance of a C-H or solvent peak is not suggested; on
/// nmrshiftdb2's assigned spectra that is about 9% of them (tests/hsqc/nmrshiftdb2_test.cpp),
/// and the D2O-shake list is the remedy. Lower the tolerance for crowded spectra.
inline constexpr core::Ppm kDefaultProtonTolerance{0.15};

/// What one 1H peak most likely is, from 1H vs HSQC (and the D2O shake, when given).
enum class ProtonVerdict : std::uint8_t {
    /// An HSQC cross peak lies within tolerance: the H is on a carbon.
    OnCarbon,
    /// A residual-solvent or water signal of the solvent lies within tolerance, and no HSQC
    /// peak does. The water shift varies (data/solvent_peaks.tsv cites its sources), so a
    /// water match means "possibly water", never certainly.
    SolventOrWater,
    /// The peak is missing from the D2O-shake list. This overrides every other evidence.
    Exchangeable,
    /// Neither an HSQC peak nor a solvent or water signal lies within tolerance: suggested
    /// as exchangeable. A 1H-13C HSQC correlates only H on carbon, so H on any other atom
    /// lands here; so does an H on carbon whose cross peak is missing from the HSQC list.
    /// The chemist decides.
    NoHsqcPartner,
    /// Both an HSQC peak and a solvent or water signal lie within tolerance: it may be a
    /// C-H, the solvent, water, or an exchangeable H hidden under one of them.
    Undetermined,
};

/// "on carbon", "solvent or water", "exchangeable (D2O)", "no HSQC partner",
/// "undetermined".
[[nodiscard]] std::string_view to_string(ProtonVerdict verdict) noexcept;

/// A solvent signal a 1H peak matched.
struct SolventSignal {
    core::Ppm shift;
    /// A water (H2O/HDO) signal rather than the solvent's residual peak.
    bool water = false;

    friend bool operator==(const SolventSignal&, const SolventSignal&) = default;
};

/// The verdict for one 1H peak, with its evidence.
struct ProtonClass {
    io::ObservedPeak peak;
    ProtonVerdict verdict = ProtonVerdict::NoHsqcPartner;
    /// HSQC cross peaks within tolerance, in HSQC input order.
    std::vector<io::CorrelationPeak> hsqc_partners;
    /// Solvent and water signals within tolerance, ascending.
    std::vector<SolventSignal> solvent_signals;
    /// Only when a D2O-shake list was given: no peak of it lies within tolerance.
    bool gone_after_d2o = false;

    friend bool operator==(const ProtonClass&, const ProtonClass&) = default;
};

/// Spectra for classify_protons. The 1H and HSQC lists and the solvent are required; the
/// 1H list after a D2O shake is optional (std::nullopt when not acquired).
struct ProtonSpectra {
    std::span<const io::ObservedPeak> proton;
    std::span<const io::CorrelationPeak> hsqc;
    nmr1d::Solvent solvent = nmr1d::Solvent::CDCl3;
    std::optional<std::span<const io::ObservedPeak>> after_d2o = std::nullopt;
    core::Ppm tolerance = kDefaultProtonTolerance;
};

/// The verdicts, and the D2O-shake peaks that match no original peak, which are reported,
/// not used.
struct ProtonClassification {
    /// One per 1H peak, in input order.
    std::vector<ProtonClass> peaks;
    std::vector<io::ObservedPeak> new_after_d2o;

    friend bool operator==(const ProtonClassification&, const ProtonClassification&) = default;
};

/// Classifies every 1H peak. The evidence is applied in this order:
///   1. gone after the D2O shake                     -> Exchangeable
///   2. an HSQC partner and a solvent/water signal   -> Undetermined
///   3. an HSQC partner only                         -> OnCarbon
///   4. a solvent/water signal only                  -> SolventOrWater
///   5. neither                                      -> NoHsqcPartner
/// "Within tolerance" compares 1H shifts only and is inclusive. The solvent's signals are
/// nmr1d::solvent_peaks(solvent).
///
/// Errors (InvalidArgument): an empty 1H or HSQC list, a non-finite shift, and a
/// non-finite or non-positive tolerance. A missing solvent row is passed on from
/// nmr1d::solvent_peaks.
[[nodiscard]] core::Result<ProtonClassification> classify_protons(const ProtonSpectra& spectra);

/// A suggested number of H on O, N, S or Se (exchangeable protons; chem::exchangeable_hydrogens
/// counts them in a candidate), for the chemist to confirm.
struct ExchangeableCount {
    /// Protons of the Exchangeable and NoHsqcPartner peaks.
    int min = 0;
    /// min plus the protons of the Undetermined peaks; std::nullopt (no upper bound) in a
    /// solvent carrying deuterium on O, N or S (nmr1d's exchanges_labile_h: D2O, CD3OD,
    /// TFE-d3), the chemist's decision (2026-11-01). The chemist's confirmed count may
    /// exceed this suggestion.
    std::optional<int> max = 0;

    friend bool operator==(const ExchangeableCount&, const ExchangeableCount&) = default;
};

/// The count suggested by classify_protons' verdicts. Each peak counts its integral rounded
/// to the nearest whole number, at least 1; a peak without an integral counts 1. The solvent
/// row is read with nmr1d::solvent_peaks; its error is passed on.
[[nodiscard]] core::Result<ExchangeableCount> suggest_exchangeable_count(
    std::span<const ProtonClass> classes, nmr1d::Solvent solvent);

}  // namespace arise::hsqc

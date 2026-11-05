#include "arise/hsqc/exchangeable.hpp"

#include "arise/core/error.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

namespace arise::hsqc {

namespace {

core::Error invalid(std::string message) {
    return {.code = core::ErrorCode::InvalidArgument, .message = std::move(message)};
}

bool within(core::Ppm a, core::Ppm b, core::Ppm tolerance) noexcept {
    return std::abs(a.value() - b.value()) <= tolerance.value();
}

template <typename Peaks, typename Shift>
bool all_finite(const Peaks& peaks, Shift shift) {
    return std::all_of(peaks.begin(), peaks.end(),
                       [&shift](const auto& peak) { return std::isfinite(shift(peak).value()); });
}

core::Result<core::Ppm> check(const ProtonSpectra& spectra) {
    if (spectra.proton.empty()) {
        return invalid("classify_protons: the 1H list has no peaks");
    }
    if (spectra.hsqc.empty()) {
        return invalid("classify_protons: the HSQC list has no peaks");
    }
    const auto tolerance = spectra.tolerance;
    if (!std::isfinite(tolerance.value()) || tolerance.value() <= 0.0) {
        return invalid("classify_protons: the matching tolerance must be a positive number");
    }
    const auto proton_shift = [](const io::ObservedPeak& peak) { return peak.shift; };
    if (!all_finite(spectra.proton, proton_shift)) {
        return invalid("classify_protons: a 1H shift is not a finite number");
    }
    if (!all_finite(spectra.hsqc, [](const io::CorrelationPeak& peak) { return peak.proton; })) {
        return invalid("classify_protons: an HSQC 1H shift is not a finite number");
    }
    if (spectra.after_d2o && !all_finite(*spectra.after_d2o, proton_shift)) {
        return invalid("classify_protons: a shift of the D2O-shake list is not a finite number");
    }
    return tolerance;
}

// The solvent's residual and water signals, ascending.
std::vector<SolventSignal> signals_of(const nmr1d::SolventPeaks& solvent) {
    std::vector<SolventSignal> signals;
    signals.reserve(solvent.residual.size() + solvent.water.size());
    for (const auto shift : solvent.residual) {
        signals.push_back({.shift = shift, .water = false});
    }
    for (const auto shift : solvent.water) {
        signals.push_back({.shift = shift, .water = true});
    }
    std::stable_sort(
        signals.begin(), signals.end(),
        [](const SolventSignal& a, const SolventSignal& b) { return a.shift < b.shift; });
    return signals;
}

ProtonVerdict verdict_of(const ProtonClass& found) noexcept {
    const bool on_carbon = !found.hsqc_partners.empty();
    const bool solvent = !found.solvent_signals.empty();
    if (found.gone_after_d2o) {
        return ProtonVerdict::Exchangeable;
    }
    if (on_carbon && solvent) {
        return ProtonVerdict::Undetermined;
    }
    if (on_carbon) {
        return ProtonVerdict::OnCarbon;
    }
    return solvent ? ProtonVerdict::SolventOrWater : ProtonVerdict::NoHsqcPartner;
}

// Protons one peak stands for: its integral rounded to the nearest whole number (halves
// up), at least 1; 1 without an integral.
int protons_of(const io::ObservedPeak& peak) {
    if (!peak.integral) {
        return 1;
    }
    return std::max(1, static_cast<int>(std::lround(*peak.integral)));
}

}  // namespace

std::string_view to_string(ProtonVerdict verdict) noexcept {
    switch (verdict) {
        case ProtonVerdict::OnCarbon:
            return "on carbon";
        case ProtonVerdict::SolventOrWater:
            return "solvent or water";
        case ProtonVerdict::Exchangeable:
            return "exchangeable (D2O)";
        case ProtonVerdict::NoHsqcPartner:
            return "no HSQC partner";
        case ProtonVerdict::Undetermined:
            return "undetermined";
    }
    return "unknown";
}

core::Result<ProtonClassification> classify_protons(const ProtonSpectra& spectra) {
    const auto checked = check(spectra);
    if (!checked) {
        return checked.error();
    }
    const auto tolerance = checked.value();
    const auto solvent = nmr1d::solvent_peaks(spectra.solvent);
    if (!solvent) {
        return solvent.error();
    }
    const auto signals = signals_of(solvent.value());

    ProtonClassification result;
    for (const auto& peak : spectra.proton) {
        ProtonClass found{
            .peak = peak,
            .verdict = ProtonVerdict::NoHsqcPartner,
            .hsqc_partners = {},
            .solvent_signals = {},
            .gone_after_d2o = false,
        };
        for (const auto& cross : spectra.hsqc) {
            if (within(peak.shift, cross.proton, tolerance)) {
                found.hsqc_partners.push_back(cross);
            }
        }
        for (const auto& signal : signals) {
            if (within(peak.shift, signal.shift, tolerance)) {
                found.solvent_signals.push_back(signal);
            }
        }
        if (spectra.after_d2o) {
            found.gone_after_d2o =
                std::none_of(spectra.after_d2o->begin(), spectra.after_d2o->end(),
                             [&](const io::ObservedPeak& after) {
                                 return within(peak.shift, after.shift, tolerance);
                             });
        }
        found.verdict = verdict_of(found);
        result.peaks.push_back(std::move(found));
    }
    if (spectra.after_d2o) {
        for (const auto& after : *spectra.after_d2o) {
            const bool known = std::any_of(spectra.proton.begin(), spectra.proton.end(),
                                           [&](const io::ObservedPeak& peak) {
                                               return within(peak.shift, after.shift, tolerance);
                                           });
            if (!known) {
                result.new_after_d2o.push_back(after);
            }
        }
    }
    return result;
}

core::Result<ExchangeableCount> suggest_exchangeable_count(std::span<const ProtonClass> classes,
                                                           nmr1d::Solvent solvent) {
    const auto peaks = nmr1d::solvent_peaks(solvent);
    if (!peaks) {
        return peaks.error();
    }
    int certain = 0;
    int possible = 0;
    for (const auto& found : classes) {
        switch (found.verdict) {
            case ProtonVerdict::Exchangeable:
            case ProtonVerdict::NoHsqcPartner:
                certain += protons_of(found.peak);
                break;
            case ProtonVerdict::Undetermined:
                possible += protons_of(found.peak);
                break;
            case ProtonVerdict::OnCarbon:
            case ProtonVerdict::SolventOrWater:
                break;
        }
    }
    return ExchangeableCount{
        .min = certain,
        .max = peaks.value().exchanges_labile_h ? std::nullopt : std::optional{certain + possible},
    };
}

}  // namespace arise::hsqc

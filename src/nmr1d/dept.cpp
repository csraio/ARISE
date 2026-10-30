#include "arise/nmr1d/dept.hpp"

#include "arise/core/error.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

namespace arise::nmr1d {
namespace {

constexpr std::array kAllTypes{CarbonType::C, CarbonType::CH, CarbonType::CH2, CarbonType::CH3};

// Shortest fixed-point form: 29.000000 -> "29", 50.100000 -> "50.1".
std::string format_ppm(core::Ppm shift) {
    std::string text = std::to_string(shift.value());
    text.erase(text.find_last_not_of('0') + 1);
    if (text.back() == '.') {
        text.pop_back();
    }
    return text;
}

core::Error invalid(std::string message) {
    return core::Error{.code = core::ErrorCode::InvalidArgument, .message = std::move(message)};
}

// Index of the single 13C peak within `tolerance` of `peak`.
core::Result<std::size_t> match_carbon(std::span<const io::CarbonPeak> carbon,
                                       const io::DeptPeak& peak, core::Ppm tolerance,
                                       std::string_view spectrum) {
    std::optional<std::size_t> found;
    for (std::size_t i = 0; i < carbon.size(); ++i) {
        if (std::fabs((peak.shift - carbon[i].shift).value()) > tolerance.value()) {
            continue;
        }
        if (found.has_value()) {
            return invalid(std::string{spectrum} + " peak at " + format_ppm(peak.shift) +
                           " ppm is ambiguous: 13C peaks at " + format_ppm(carbon[*found].shift) +
                           " and " + format_ppm(carbon[i].shift) + " ppm are both within " +
                           format_ppm(tolerance) + " ppm; lower the matching tolerance");
        }
        found = i;
    }
    if (!found.has_value()) {
        return invalid(std::string{spectrum} + " peak at " + format_ppm(peak.shift) +
                       " ppm has no 13C peak within " + format_ppm(tolerance) +
                       " ppm; every DEPT peak must also appear in the 13C list");
    }
    return *found;
}

// For each 13C peak, the DEPT peak matched to it, if any.
core::Result<std::vector<std::optional<io::DeptPeak>>> match_spectrum(
    std::span<const io::CarbonPeak> carbon, std::span<const io::DeptPeak> dept, core::Ppm tolerance,
    std::string_view spectrum) {
    std::vector<std::optional<io::DeptPeak>> matched(carbon.size());
    for (const auto& peak : dept) {
        auto index = match_carbon(carbon, peak, tolerance, spectrum);
        if (!index.ok()) {
            return std::move(index).error();
        }
        auto& slot = matched[index.value()];
        if (slot.has_value()) {
            return invalid(std::string{spectrum} + " peaks at " + format_ppm(slot->shift) +
                           " and " + format_ppm(peak.shift) + " ppm both match the 13C peak at " +
                           format_ppm(carbon[index.value()].shift) + " ppm");
        }
        slot = peak;
    }
    return matched;
}

// DEPT-135: CH and CH3 positive, CH2 negative, C absent.
CarbonTypeSet from_dept135(const std::optional<io::DeptPeak>& peak) {
    if (!peak.has_value()) {
        return {CarbonType::C};
    }
    return peak->intensity > 0.0 ? CarbonTypeSet{CarbonType::CH, CarbonType::CH3}
                                 : CarbonTypeSet{CarbonType::CH2};
}

// DEPT-90: only CH, positive.
CarbonTypeSet from_dept90(const std::optional<io::DeptPeak>& peak) {
    if (!peak.has_value()) {
        return {CarbonType::C, CarbonType::CH2, CarbonType::CH3};
    }
    return {CarbonType::CH};
}

}  // namespace

std::string_view to_string(CarbonType type) noexcept {
    switch (type) {
        case CarbonType::C:
            return "C";
        case CarbonType::CH:
            return "CH";
        case CarbonType::CH2:
            return "CH2";
        case CarbonType::CH3:
            return "CH3";
    }
    return "?";
}

std::string to_string(CarbonTypeSet set) {
    std::string text = "{";
    for (const auto type : kAllTypes) {
        if (set.contains(type)) {
            text += (text.size() > 1 ? ", " : "") + std::string{to_string(type)};
        }
    }
    return text + "}";
}

core::Result<std::vector<CarbonAssignment>> classify_carbons(const DeptSpectra& spectra) {
    const auto carbon = spectra.carbon;
    if (carbon.empty()) {
        return invalid("the 13C peak list is empty");
    }
    if (!std::isfinite(spectra.tolerance.value()) || spectra.tolerance.value() <= 0.0) {
        return invalid("matching tolerance must be finite and > 0, got " +
                       std::to_string(spectra.tolerance.value()) + " ppm");
    }

    std::vector<CarbonTypeSet> by_dept135(carbon.size(), CarbonTypeSet::all());
    if (spectra.dept135.has_value()) {
        auto matched = match_spectrum(carbon, *spectra.dept135, spectra.tolerance, "DEPT-135");
        if (!matched.ok()) {
            return std::move(matched).error();
        }
        for (std::size_t i = 0; i < carbon.size(); ++i) {
            by_dept135[i] = from_dept135(matched.value()[i]);
        }
    }

    std::vector<CarbonTypeSet> by_dept90(carbon.size(), CarbonTypeSet::all());
    if (spectra.dept90.has_value()) {
        auto matched = match_spectrum(carbon, *spectra.dept90, spectra.tolerance, "DEPT-90");
        if (!matched.ok()) {
            return std::move(matched).error();
        }
        for (std::size_t i = 0; i < carbon.size(); ++i) {
            const auto& peak = matched.value()[i];
            if (peak.has_value() && peak->intensity < 0.0) {
                return invalid("DEPT-90 peak at " + format_ppm(peak->shift) +
                               " ppm is negative; DEPT-90 shows only CH, which is positive "
                               "(check the phasing)");
            }
            by_dept90[i] = from_dept90(peak);
        }
    }

    std::vector<CarbonAssignment> assignments;
    assignments.reserve(carbon.size());
    for (std::size_t i = 0; i < carbon.size(); ++i) {
        const auto possible = by_dept135[i] & by_dept90[i];
        if (possible.empty()) {
            return invalid("13C peak at " + format_ppm(carbon[i].shift) + " ppm: DEPT-90 implies " +
                           to_string(by_dept90[i]) + " but DEPT-135 implies " +
                           to_string(by_dept135[i]) + "; the spectra contradict each other");
        }
        assignments.push_back(CarbonAssignment{.shift = carbon[i].shift, .possible = possible});
    }
    return assignments;
}

}  // namespace arise::nmr1d

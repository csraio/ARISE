#pragma once

#include "arise/core/result.hpp"
#include "arise/core/units.hpp"
#include "arise/io/carbon_csv.hpp"

#include <bit>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace arise::nmr1d {

/// A carbon by its number of attached hydrogens, which is what DEPT measures. The
/// enumerator value is that number. These are often nicknamed quaternary, tertiary,
/// secondary and primary, but that naming strictly counts attached carbons, which DEPT
/// cannot see (an aromatic CH has two carbon neighbours; a C=O often has two).
enum class CarbonType : std::uint8_t { C = 0, CH = 1, CH2 = 2, CH3 = 3 };

/// "C", "CH", "CH2" or "CH3".
[[nodiscard]] std::string_view to_string(CarbonType type) noexcept;

/// The carbon types still possible for one 13C peak.
class CarbonTypeSet {
public:
    /// The empty set.
    constexpr CarbonTypeSet() noexcept = default;
    constexpr CarbonTypeSet(std::initializer_list<CarbonType> types) noexcept {
        for (const auto type : types) {
            bits_ = static_cast<std::uint8_t>(bits_ | bit(type));
        }
    }

    [[nodiscard]] static constexpr CarbonTypeSet all() noexcept {
        return {CarbonType::C, CarbonType::CH, CarbonType::CH2, CarbonType::CH3};
    }

    [[nodiscard]] constexpr bool contains(CarbonType type) const noexcept {
        return (bits_ & bit(type)) != 0;
    }
    [[nodiscard]] constexpr bool empty() const noexcept { return bits_ == 0; }
    [[nodiscard]] constexpr int size() const noexcept { return std::popcount(bits_); }
    /// The type, when exactly one remains.
    [[nodiscard]] constexpr std::optional<CarbonType> resolved() const noexcept {
        for (const auto type : {CarbonType::C, CarbonType::CH, CarbonType::CH2, CarbonType::CH3}) {
            if (size() == 1 && contains(type)) {
                return type;
            }
        }
        return std::nullopt;
    }

    /// Intersection: the types allowed by both.
    [[nodiscard]] friend constexpr CarbonTypeSet operator&(CarbonTypeSet a,
                                                           CarbonTypeSet b) noexcept {
        CarbonTypeSet both;
        both.bits_ = static_cast<std::uint8_t>(a.bits_ & b.bits_);
        return both;
    }
    /// Union: the types allowed by either.
    [[nodiscard]] friend constexpr CarbonTypeSet operator|(CarbonTypeSet a,
                                                           CarbonTypeSet b) noexcept {
        CarbonTypeSet either;
        either.bits_ = static_cast<std::uint8_t>(a.bits_ | b.bits_);
        return either;
    }
    friend constexpr bool operator==(CarbonTypeSet, CarbonTypeSet) noexcept = default;

private:
    static constexpr std::uint8_t bit(CarbonType type) noexcept {
        return static_cast<std::uint8_t>(1U << static_cast<unsigned>(type));
    }

    std::uint8_t bits_ = 0;
};

/// "{CH, CH3}" style, in C, CH, CH2, CH3 order; "{}" when empty.
[[nodiscard]] std::string to_string(CarbonTypeSet set);

/// The classification of one 13C peak.
struct CarbonAssignment {
    core::Ppm shift;
    CarbonTypeSet possible;

    friend bool operator==(const CarbonAssignment&, const CarbonAssignment&) = default;
};

/// How far a DEPT peak's shift may differ from its 13C peak and still match. A
/// convention, not a constant of nature: spectra of one sample, referenced alike,
/// usually agree within about 0.01 ppm, so 0.05 ppm leaves headroom while keeping
/// carbons about 0.1 ppm apart distinct. Lower it for crowded spectra.
inline constexpr core::Ppm kDefaultMatchTolerance{0.05};

/// Spectra for classify_carbons. The 13C list is required; each DEPT is optional
/// (std::nullopt when not acquired; an empty span is a DEPT with no peaks).
struct DeptSpectra {
    std::span<const io::CarbonPeak> carbon;
    std::optional<std::span<const io::DeptPeak>> dept135 = std::nullopt;
    std::optional<std::span<const io::DeptPeak>> dept90 = std::nullopt;
    core::Ppm tolerance = kDefaultMatchTolerance;
};

/// Classifies every 13C peak by attached hydrogens, returning one assignment per
/// 13C peak in 13C order. Each DEPT peak must match exactly one 13C peak within
/// `tolerance` (inclusive), and each 13C peak at most one peak per DEPT spectrum.
///
///   DEPT-135: positive -> {CH, CH3}, negative -> {CH2}, absent -> {C}
///   DEPT-90:  positive -> {CH},                         absent -> {C, CH2, CH3}
///
/// The results of the supplied spectra are intersected. Errors (InvalidArgument):
/// empty 13C list, non-finite or non-positive tolerance, a DEPT peak with no or several
/// 13C matches, two DEPT peaks on one 13C peak, a negative DEPT-90 peak, and spectra
/// that contradict each other (an empty intersection).
///
/// Not detected: a DEPT-135 phased upside down when no DEPT-90 is given, and 13C peaks
/// that are not from the sample (a CDCl3 solvent peak is correctly reported as C).
[[nodiscard]] core::Result<std::vector<CarbonAssignment>> classify_carbons(
    const DeptSpectra& spectra);

}  // namespace arise::nmr1d

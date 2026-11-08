#pragma once

#include "arise/chem/molecule.hpp"
#include "arise/core/result.hpp"
#include "arise/filter/groups.hpp"

#include <optional>

namespace arise::filter {

/// The number of exchangeable hydrogens (H on O, N, S, Se; chem::exchangeable_hydrogens)
/// the chemist has confirmed, as an inclusive range. One exact count is min == max. No
/// upper bound (std::nullopt) when exchangeable protons may be hidden, as in D2O or CD3OD
/// (hsqc::suggest_exchangeable_count suggests a range; the chemist confirms it).
struct ConfirmedExchangeable {
    int min = 0;
    std::optional<int> max;

    friend bool operator==(const ConfirmedExchangeable&, const ConfirmedExchangeable&) = default;
};

/// Checks a candidate against the confirmed count: Eliminated, with the reason, when its
/// chem::exchangeable_hydrogens falls outside [min, max]; Active otherwise. A confirmed
/// count is a hard fact, so this runs with the confirmed groups, before DEPT tiering.
/// InvalidArgument if min is negative or max is below min.
[[nodiscard]] core::Result<Verdict> check_confirmed_exchangeable(
    const chem::Molecule& molecule, const ConfirmedExchangeable& confirmed);

}  // namespace arise::filter

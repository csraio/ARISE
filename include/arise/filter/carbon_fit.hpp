#pragma once

#include "arise/chem/molecule.hpp"
#include "arise/core/result.hpp"
#include "arise/filter/groups.hpp"
#include "arise/filter/hybridization.hpp"
#include "arise/nmr1d/dept.hpp"

#include <cstddef>
#include <span>
#include <string>
#include <vector>

namespace arise::filter {

/// How one candidate's carbons fit the chemist's 13C list (plan 3c). Soft evidence: a
/// candidate is never Eliminated here, only Active or BackBurner, and assumed overlaps only
/// order the list (the chemist's decisions, 2026-11-09).
struct CarbonFit {
    /// Active, or BackBurner with one reason per conflict.
    Verdict verdict;
    /// How many of the candidate's always-split carbon environments must coincide for the
    /// list to have as few peaks as it does: max(0, fewest - peaks). Orders, never tiers.
    int assumed_overlaps = 0;
    /// The candidate's carbon environments: always split (chem::nmr_equivalence().fewest)
    /// and at most (.most).
    int environments_fewest = 0;
    int environments_most = 0;
    /// Remarks that do not change the tier, e.g. "assumes 1 overlapping peak".
    std::vector<std::string> notes;

    friend bool operator==(const CarbonFit&, const CarbonFit&) = default;
};

/// Checks one candidate against the chemist's classified 13C peaks (nmr1d::classify_carbons):
/// - more peaks than the candidate's `most` environments: BackBurner ("13C list has N
///   peaks; the candidate has at most M carbon environments");
/// - with options.use_dept on, each peak must take a distinct `most` environment whose
///   carbon type (chem::carbon_type) is among the peak's possible types, found by a maximum
///   matching; every peak left over is a BackBurner reason naming its shift and types. With
///   use_dept off, or peaks of unknown type (no DEPT list: every type possible), only the
///   count is compared. An environment without a DEPT type (methane) accepts any peak;
/// - assumed_overlaps = max(0, fewest - peaks), with a note when positive.
/// The same options as suggest_groups, so one tick-box drives both.
/// InvalidArgument for an empty peak list.
[[nodiscard]] core::Result<CarbonFit> check_carbon_fit(
    const chem::Molecule& molecule, std::span<const nmr1d::CarbonAssignment> peaks,
    SuggestOptions options = {});

/// check_carbon_fit, then the chemist's confirmed minimum counts (plan 3c-b): the candidate's
/// `most` carbon environments must be able to supply every confirmed count at once. An
/// environment counts for the hybridization (chem::carbon_hybridization) and DEPT type
/// (chem::carbon_type) of its carbons; one whose carbons have no single hybridization
/// (charged carbons) or no DEPT type (methane) stands in for any one type. Hybridization and
/// DEPT-type counts are checked separately. Each type the environments cannot supply is a
/// BackBurner reason ("confirmed at least 3 sp2 peaks; the candidate has at most 2 sp2 carbon
/// environments"), and so is a set of counts that stand-in environments cannot all cover
/// ("together"). Maximums are never used against a candidate. Reasons follow those of the
/// peak check. InvalidArgument for a negative count, and as check_carbon_fit.
[[nodiscard]] core::Result<CarbonFit> check_carbon_fit(
    const chem::Molecule& molecule, std::span<const nmr1d::CarbonAssignment> peaks,
    const ConfirmedCarbonCounts& confirmed, SuggestOptions options = {});

/// Candidate indices in ranking order: Active before BackBurner (and Eliminated last),
/// then fewer assumed overlaps, then fewer reasons, then input order (stable).
[[nodiscard]] std::vector<std::size_t> order_by_fit(std::span<const CarbonFit> fits);

}  // namespace arise::filter

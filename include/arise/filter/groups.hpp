#pragma once

#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "arise/core/result.hpp"
#include "arise/core/units.hpp"
#include "arise/nmr1d/dept.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace arise::filter {

/// The 13C window of one hallmark of one group: a peak inside it (inclusive), from a carbon
/// whose DEPT type is among `types`, suggests the group.
struct HallmarkWindow {
    chem::Group group = chem::Group::CarboxylicAcid;
    /// As chem::hallmarks names it, plus the carbon type where measured separately
    /// ("O-C carbon (sp3)").
    std::string hallmark;
    /// The DEPT types of the carbons the window was measured on.
    nmr1d::CarbonTypeSet types;
    core::Ppm low;
    core::Ppm high;
    /// How many carbons the window was measured on.
    std::size_t n = 0;

    friend bool operator==(const HallmarkWindow&, const HallmarkWindow&) = default;
};

/// Reads hallmark windows from a table in the format of tests/data/reference/group_shifts.tsv
/// (written by tools/measure_group_shifts): '#' comment lines, then tab-separated rows of
/// group name (as chem::to_string), hallmark, DEPT types ("{C, CH}"), n, min, 1st, 5th
/// percentile, median, 95th, 99th percentile, max. A window runs from the 1st to the 99th
/// percentile (user decision 2026-11-01). The engine holds no table of its own: the caller
/// supplies it, so its source and licence stay the caller's choice.
///
/// Errors (ParseError, citing the line): a row without 11 fields, an unknown group, empty or
/// unknown DEPT types, a number that is not a plain decimal or not finite, n of 0, 1st
/// percentile above the 99th, a group and hallmark listed twice, and a table with no rows.
[[nodiscard]] core::Result<std::vector<HallmarkWindow>> parse_hallmark_windows(
    std::string_view text);

/// One peak that falls in one window of a suggested group.
struct PeakEvidence {
    /// Index of the peak in the list given to suggest_groups.
    std::size_t peak = 0;
    core::Ppm shift;
    HallmarkWindow window;

    friend bool operator==(const PeakEvidence&, const PeakEvidence&) = default;
};

/// A group the spectrum suggests, with every peak and window that suggest it.
struct GroupSuggestion {
    chem::Group group = chem::Group::CarboxylicAcid;
    /// By peak index, then window order.
    std::vector<PeakEvidence> evidence;

    friend bool operator==(const GroupSuggestion&, const GroupSuggestion&) = default;
};

/// The groups the observed 13C peaks suggest, in chem::kAllGroups order: a group is
/// suggested when any peak lies in any of its windows (bounds inclusive) and the peak's
/// possible DEPT types share a type with the window's. Suggestions are never
/// eliminations: overlapping windows suggest every group they belong to, and the chemist
/// decides. A peak with unknown DEPT type (all types possible) matches on shift alone.
[[nodiscard]] std::vector<GroupSuggestion> suggest_groups(
    std::span<const nmr1d::CarbonAssignment> peaks, std::span<const HallmarkWindow> windows);

/// Groups the chemist has confirmed as present or ruled out.
struct ConfirmedGroups {
    std::vector<chem::Group> present;
    std::vector<chem::Group> absent;
};

/// Where a candidate stands (plan 3: rank, don't delete).
enum class Tier : std::uint8_t {
    Active,      // consistent with all evidence so far
    BackBurner,  // soft conflicts; kept, ranked lower, revisited later
    Eliminated,  // impossible: violates a hard fact, such as a confirmed group
};

/// "active", "back-burner" or "eliminated".
[[nodiscard]] std::string_view to_string(Tier tier) noexcept;

/// A candidate's tier and the reasons for it (empty when Active).
struct Verdict {
    Tier tier = Tier::Active;
    std::vector<std::string> reasons;

    friend bool operator==(const Verdict&, const Verdict&) = default;
};

/// Checks a candidate against the chemist's confirmed groups: Eliminated, with one reason
/// per violated group (in the order given), if it lacks a group confirmed present or has
/// one confirmed absent (counted by chem::count_group); Active otherwise.
/// InvalidArgument if a group is confirmed both present and absent.
[[nodiscard]] core::Result<Verdict> check_confirmed_groups(const chem::Molecule& molecule,
                                                           const ConfirmedGroups& confirmed);

}  // namespace arise::filter

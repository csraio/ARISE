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

/// Hallmarks measured on fewer molecules than this use their full measured range (min to
/// max) instead of the 1st to 99th percentile: below 100 molecules a percentile rests on
/// less than one molecule (user decision 2026-11-02). A convention, not a constant of
/// nature.
inline constexpr std::size_t kMinMoleculesForPercentiles = 100;

/// Where a window's bounds come from, for showing the evidence.
enum class WindowSource : std::uint8_t {
    Percentiles,    // 1st to 99th percentile of the measured shifts
    MeasuredRange,  // min to max of the measured shifts (fewer than 100 molecules)
    Lab,            // the chemist's own table
};

/// The 13C window of one hallmark of one group: a peak inside it (inclusive), from a carbon
/// whose DEPT type is among `types`, suggests the group.
struct HallmarkWindow {
    chem::Group group = chem::Group::CarboxylicAcid;
    /// As chem::hallmarks names it, plus the carbon type where measured separately
    /// ("O-C carbon (sp3)"); for a lab window, the lab's label.
    std::string hallmark;
    /// The DEPT types of the carbons the window was measured on (a lab's choice for a lab
    /// window; all four when it gives none).
    nmr1d::CarbonTypeSet types;
    core::Ppm low;
    core::Ppm high;
    /// How many carbons, and how many distinct molecules, the window was measured on
    /// (0 for a lab window).
    std::size_t n = 0;
    std::size_t molecules = 0;
    WindowSource source = WindowSource::Percentiles;

    friend bool operator==(const HallmarkWindow&, const HallmarkWindow&) = default;
};

/// Reads hallmark windows from a table in the format of data/group_shifts.tsv (written by
/// tools/measure_group_shifts): '#' comment lines, then tab-separated rows of group name
/// (as chem::to_string), hallmark, DEPT types ("{C, CH}"), n, molecules, min, 1st, 5th
/// percentile, median, 95th, 99th percentile, max. A window runs from the 1st to the 99th
/// percentile (user decision 2026-11-01), or from min to max when the row has fewer than
/// kMinMoleculesForPercentiles molecules.
///
/// Errors (ParseError, citing the line): a row without 12 fields, an unknown group, empty or
/// unknown DEPT types, a count that is not a positive whole number, a number that is not a
/// plain decimal or not finite, a 1st percentile above the 99th or a min above the max, a
/// group and hallmark listed twice, and a table with no rows.
[[nodiscard]] core::Result<std::vector<HallmarkWindow>> parse_hallmark_windows(
    std::string_view text);

/// The windows ARISE ships with: data/group_shifts.tsv, compiled in at configure time and
/// read with parse_hallmark_windows. Contains information from nmrshiftdb2
/// (www.nmrshiftdb.org), made available under the nmrshiftdb2 Database License.
[[nodiscard]] core::Result<std::vector<HallmarkWindow>> measured_hallmark_windows();

/// Reads a lab's own windows (io::parse_window_csv: columns group, low_ppm, high_ppm and
/// optionally hallmark and dept). Group names are matched to chem::to_string ignoring
/// case. `dept` lists types separated by spaces (C, CH, CH2, CH3, any case); blank means
/// any type. A blank `hallmark` becomes the group's hallmark names joined by " or "
/// ("13C window" for a group without one). Errors cite the line: an unknown group (the
/// message lists the valid names), an unknown DEPT type, a group and hallmark listed
/// twice, and the reader's own errors.
[[nodiscard]] core::Result<std::vector<HallmarkWindow>> parse_lab_windows(std::string_view text);

/// The built-in windows with a lab's applied: every group the lab lists keeps only the
/// lab's windows; every other group keeps its built-in ones. Ordered by group
/// (chem::kAllGroups), then as given.
[[nodiscard]] std::vector<HallmarkWindow> merge_windows(std::span<const HallmarkWindow> built_in,
                                                        std::span<const HallmarkWindow> lab);

/// The built-in windows as a lab table (parse_lab_windows' format, with comment lines
/// explaining it), for the chemist to download, edit and upload: shifts to one decimal, as
/// measured. Rows the chemist deletes keep ARISE's windows.
[[nodiscard]] std::string lab_table_template();

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

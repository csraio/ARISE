#pragma once
// Group suggestions from the 1H spectrum (the 1H plan; the chemist's decisions, 2026-11-22).
// One 1H peak in a crowded range says little, so every window carries its measured
// precision, and only the windows the chemist enabled for a round suggest anything (round 1:
// methoxy O-CH3 and benzene ring CH). A group confirmed from 1H evidence is soft: it goes in
// ConfirmedGroups::likely_present (filter/groups.hpp), which back-burners, never eliminates.

#include "arise/chem/perception.hpp"
#include "arise/core/result.hpp"
#include "arise/core/units.hpp"
#include "arise/io/peak_csv.hpp"
#include "arise/nmr1d/dept.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace arise::filter {

/// The 1H window of the H on one hallmark carbon type of one group, measured from nmrshiftdb2
/// (data/group_proton_shifts.tsv, written by tools/measure_group_proton_shifts).
struct ProtonWindow {
    chem::Group group = chem::Group::CarboxylicAcid;
    std::string hallmark;  // e.g. "O-CH3 carbon", "ring carbons", "O-C carbon (sp3)"
    /// The key carbon's H count: a peak must integrate to at least this many protons
    /// (CH 1, CH2 2, CH3 3) when ProtonOptions::use_counts is on.
    nmr1d::CarbonType type = nmr1d::CarbonType::CH;
    core::Ppm low;
    core::Ppm high;
    std::size_t n = 0;          // carbons measured
    std::size_t molecules = 0;  // molecules measured
    /// Of the measured molecules with an H of this count inside the window, the share that
    /// have the group; and the share of all measured molecules that have it.
    double precision = 0;
    double base_rate = 0;
    /// Whether the chemist enabled the window for round 1.
    bool round1 = false;

    friend bool operator==(const ProtonWindow&, const ProtonWindow&) = default;
};

/// Reads windows in the format of data/group_proton_shifts.tsv: '#' comment lines, then
/// tab-separated rows of group name (as chem::to_string), hallmark, H count (CH, CH2 or CH3),
/// n, molecules, min, 1st, 5th percentile, median, 95th, 99th percentile, max, window low,
/// window high, hits, true hits, precision, molecules with the group, molecules, base rate,
/// round ("1" or "-"). ParseError naming the line for anything else: an unknown group or H
/// count, a malformed number, low above high, a duplicate (group, hallmark, H count).
[[nodiscard]] core::Result<std::vector<ProtonWindow>> parse_proton_windows(std::string_view text);

/// The windows measured from nmrshiftdb2, compiled into the engine.
[[nodiscard]] core::Result<std::vector<ProtonWindow>> measured_proton_windows();

/// The number of protons a 1H peak's integral stands for: the integral rounded to the nearest
/// whole number (halves up), at least 1; none without an integral.
[[nodiscard]] std::optional<int> proton_count(const io::ObservedPeak& peak);

/// One 1H peak inside one window.
struct ProtonEvidence {
    /// Index of the peak in the list given to suggest_groups_from_protons.
    std::size_t peak = 0;
    core::Ppm shift;
    std::optional<int> protons;  // proton_count of the peak
    ProtonWindow window;

    friend bool operator==(const ProtonEvidence&, const ProtonEvidence&) = default;
};

/// A group the 1H spectrum suggests, with every peak and window that suggest it.
struct ProtonSuggestion {
    chem::Group group = chem::Group::CarboxylicAcid;
    /// By peak index, then window order.
    std::vector<ProtonEvidence> evidence;

    friend bool operator==(const ProtonSuggestion&, const ProtonSuggestion&) = default;
};

/// The chemist's choices for suggest_groups_from_protons.
struct ProtonOptions {
    /// Whether a peak's proton count must reach the window's H count (on by default). Off,
    /// every peak matches on shift alone; a peak without an integral always does.
    bool use_counts = true;
};

/// The groups the 1H peaks suggest, in chem::kAllGroups order, from the round-1 windows only:
/// a group is suggested when any peak lies in any of its enabled windows (bounds inclusive)
/// and, with `options.use_counts`, the peak's proton count is at least the window's H count.
/// Suggestions are never eliminations; the chemist decides.
[[nodiscard]] std::vector<ProtonSuggestion> suggest_groups_from_protons(
    std::span<const io::ObservedPeak> peaks, std::span<const ProtonWindow> windows,
    ProtonOptions options = {});

}  // namespace arise::filter

// The suggestion windows checked on data they were not measured from (plan 3b): the
// windows come from nmrshiftdb2 records with even IDs (data/group_shifts.tsv); here every
// record with an odd ID is the hold-out. A window runs from the 1st to the 99th percentile
// (or over the full measured range for hallmarks from fewer than 100 molecules), so at
// least about 98% of hold-out hallmark carbons should fall in
// it. The pass mark is that rate less three binomial standard errors, sqrt(0.98 * 0.02 / n)
// (a stated convention: it fails a real shortfall without failing small groups on chance).
// n counts molecules (nmrshiftdb2 IDs), not carbons or spectra: the carbons of one molecule
// (a ring's six, an anhydride's two) and its repeat spectra (the extract lists 3,349
// molecules more than once, 488 spectra exactly twice) are not independent.
// DEPT types of the hold-out carbons come from their structures, as a complete DEPT
// (DEPT-90 and DEPT-135) would give them.
#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "arise/filter/groups.hpp"
#include "arise/nmr1d/dept.hpp"
#include "group_shifts.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using arise::chem::Group;

constexpr double kExpectedRate = 0.98;

// About 30 s natively, far slower under the wasm checks; the logic is portable and runs in
// every native configuration. A run-time skip keeps the code compiled under wasm too.
#ifdef __EMSCRIPTEN__
constexpr bool kSkipSlowChecks = true;
#else
constexpr bool kSkipSlowChecks = false;
#endif

double pass_mark(std::size_t n) {
    return kExpectedRate -
           (3.0 * std::sqrt(kExpectedRate * (1.0 - kExpectedRate) / static_cast<double>(n)));
}

struct Tally {
    std::size_t n = 0;
    std::size_t hits = 0;
    std::set<long long> molecules;  // nmrshiftdb2 IDs, the independent units
};

struct HoldoutTallies {
    std::size_t records = 0;
    // Per (group, hallmark): hallmark carbons whose own peak alone suggests the group by
    // that hallmark's window.
    std::map<std::pair<Group, std::string>, Tally> carbons;
    // Per group: hold-out spectra of molecules with the group (and an assigned hallmark
    // carbon) whose full peak list suggests it.
    std::map<Group, Tally> molecules;
};

const std::vector<arise::filter::HallmarkWindow>& windows() {
    static const auto parsed = arise::filter::parse_hallmark_windows(
        arise::test::join_chunks(arise::test::kGroupShiftsChunks));
    return parsed.value();
}

bool suggests(const std::vector<arise::filter::GroupSuggestion>& suggestions, Group group,
              const std::string& hallmark) {
    for (const auto& suggestion : suggestions) {
        if (suggestion.group != group) {
            continue;
        }
        for (const auto& evidence : suggestion.evidence) {
            if (hallmark.empty() || evidence.window.hallmark == hallmark) {
                return true;
            }
        }
    }
    return false;
}

void tally_record(const arise::tools::ShiftRecord& record, const arise::chem::Molecule& molecule,
                  HoldoutTallies& tallies) {
    std::vector<arise::nmr1d::CarbonAssignment> peaks;
    constexpr auto kNoPeak = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> peak_of(record.shifts.size(), kNoPeak);
    for (std::size_t atom = 0; atom < record.shifts.size(); ++atom) {
        const auto& shift = record.shifts[atom];
        const auto type = arise::chem::carbon_type(molecule, atom);
        if (shift && type) {
            peak_of[atom] = peaks.size();
            peaks.push_back({
                .shift = arise::core::Ppm{*shift},
                .possible = arise::nmr1d::CarbonTypeSet{*type},
            });
        }
    }
    const auto all = arise::filter::suggest_groups(peaks, windows());
    std::set<Group> present;
    for (const auto& key : arise::tools::key_carbons(molecule)) {
        if (peak_of[key.atom] == kNoPeak) {
            continue;  // carbon not assigned in the record
        }
        const auto group = arise::chem::kAllGroups[key.group_index];
        present.insert(group);
        const auto one = arise::filter::suggest_groups(
            std::span{peaks}.subspan(peak_of[key.atom], 1), windows());
        auto& tally = tallies.carbons[{group, key.hallmark}];
        tally.molecules.insert(record.id);
        ++tally.n;
        if (suggests(one, group, key.hallmark)) {
            ++tally.hits;
        }
    }
    for (const auto group : present) {
        auto& tally = tallies.molecules[group];
        tally.molecules.insert(record.id);
        ++tally.n;
        if (suggests(all, group, "")) {
            ++tally.hits;
        }
    }
}

const HoldoutTallies& holdout() {
    static const HoldoutTallies tallies = [] {
        HoldoutTallies result;
        const auto extract = arise::test::join_chunks(arise::test::kNmrshiftdbAllChunks);
        for (const auto& record : arise::tools::parse_shift_extract(extract)) {
            if (record.id % 2 == 0) {
                continue;  // measured half
            }
            const auto molecule = arise::chem::parse_smiles(record.smiles);
            if (molecule.ok() && molecule.value().atoms().size() == record.shifts.size()) {
                ++result.records;
                tally_record(record, molecule.value(), result);
            }
        }
        return result;
    }();
    return tallies;
}

// Shortfalls the chemist signs off: tests/data/reference/known_differences.tsv rows with
// source "holdout" (molecule = group name, check = hallmark, or "molecules" for the
// per-molecule rate). A listed row must still fall short, so it cannot linger.
std::set<std::pair<std::string, std::string>> listed_shortfalls() {
    std::set<std::pair<std::string, std::string>> listed;
    const auto text = arise::test::join_chunks(arise::test::kKnownDifferencesChunks);
    std::string_view rest = text;
    while (!rest.empty()) {
        const auto end = rest.find('\n');
        const auto line = rest.substr(0, end);
        rest = end == std::string_view::npos ? std::string_view{} : rest.substr(end + 1);
        if (!line.starts_with("holdout\t")) {
            continue;
        }
        const auto group = line.substr(8, line.find('\t', 8) - 8);
        const auto check_start = 8 + group.size() + 1;
        const auto check = line.substr(check_start, line.find('\t', check_start) - check_start);
        listed.insert({std::string{group}, std::string{check}});
    }
    return listed;
}

// Group name and hallmark (or "molecules" for the per-molecule rate), as listed.
std::pair<std::string, std::string> shortfall_id(Group group) {
    return {std::string{arise::chem::to_string(group)}, "molecules"};
}

std::pair<std::string, std::string> shortfall_id(const std::pair<Group, std::string>& key) {
    return {std::string{arise::chem::to_string(key.first)}, key.second};
}

// Prints the rate and holds it to the pass mark (or, if listed, to falling short).
void expect_rate(const std::pair<std::string, std::string>& id, const Tally& tally,
                 const std::string& what, bool listed) {
    const double rate = static_cast<double>(tally.hits) / static_cast<double>(tally.n);
    const std::string name = id.first + " / " + id.second;
    const double mark = pass_mark(tally.molecules.size());
    std::cout << "[ holdout  ] " << what << " " << name << ": " << tally.hits << "/" << tally.n
              << " = " << rate << " (" << tally.molecules.size() << " molecules, pass mark " << mark
              << ")\n";
    if (listed) {
        EXPECT_LT(rate, mark) << name << " is listed as a shortfall but now passes";
    } else {
        EXPECT_GE(rate, mark) << name;
    }
}

// Every rate reaches its pass mark, except the listed shortfalls, which must still fall
// short; each listed row of this kind ("molecules" or a hallmark) must have been checked.
template <typename Key>
void expect_rates(const std::map<Key, Tally>& tallies, const std::string& what) {
    auto listed = listed_shortfalls();
    for (const auto& [key, tally] : tallies) {
        const auto id = shortfall_id(key);
        expect_rate(id, tally, what, listed.erase(id) == 1);
    }
    for (const auto& [group, check] : listed) {
        EXPECT_NE(check == "molecules", what == "molecule")
            << group << " / " << check << " is listed but was not checked";
    }
}

TEST(GroupSuggestionHoldout, HallmarkCarbonsFallInTheirWindowsAtTheMeasuredRate) {
    if (kSkipSlowChecks) {
        GTEST_SKIP() << "about 30 s natively, far slower under wasm checks; checked natively";
    }
    ASSERT_GT(holdout().records, 15000U);
    // Every measured hallmark has hold-out carbons, and every hold-out hallmark a window.
    std::set<std::pair<Group, std::string>> measured;
    for (const auto& window : windows()) {
        measured.insert({window.group, window.hallmark});
    }
    std::set<std::pair<Group, std::string>> held_out;
    for (const auto& [key, tally] : holdout().carbons) {
        held_out.insert(key);
    }
    EXPECT_EQ(held_out, measured);
    expect_rates(holdout().carbons, "carbon");
}

TEST(GroupSuggestionHoldout, MoleculesWithAGroupHaveItSuggestedAtTheMeasuredRate) {
    if (kSkipSlowChecks) {
        GTEST_SKIP() << "about 30 s natively, far slower under wasm checks; checked natively";
    }
    EXPECT_GT(holdout().molecules.size(), 20U);
    expect_rates(holdout().molecules, "molecule");
}

}  // namespace

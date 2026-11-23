// Hold-out check of the round-1 1H windows (the 1H plan): the windows were measured on the
// even-ID nmrshiftdb2 records (data/group_proton_shifts.tsv); here they are applied to the
// odd-ID half, which they were not measured on.
//   - Recall: the share of the group's key carbons (with that H count) whose 1H shift falls
//     in the window must reach 98% less three binomial standard errors over the molecules,
//     as for the 13C windows (tests/filter/holdout_test.cpp).
//   - Precision: of the odd-ID molecules with an H of that count inside the window, the share
//     with the group must reach the even-ID precision less three binomial standard errors.

#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "arise/filter/proton_groups.hpp"
#include "group_shifts.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

using arise::filter::ProtonWindow;

constexpr double kExpectedRecall = 0.98;

double pass_mark(double rate, std::size_t n) {
    return rate - (3.0 * std::sqrt(rate * (1.0 - rate) / static_cast<double>(n)));
}

bool inside(double shift, const ProtonWindow& window) {
    return window.low.value() <= shift && shift <= window.high.value();
}

// The hold-out counts of one window over the odd-ID records.
struct HoldOut {
    std::size_t key_carbons = 0;
    std::size_t in_window = 0;
    std::set<long long> molecules;  // with a key carbon of the window
    std::size_t hits = 0;           // molecules with an H of the window's count inside it
    std::size_t true_hits = 0;      // of those, with the group
};

// Adds one odd-ID record to `counts`.
void add_record(const arise::tools::ShiftRecord& record, const arise::chem::Molecule& molecule,
                const ProtonWindow& window, HoldOut& counts) {
    bool hit = false;
    for (std::size_t atom = 0; atom < record.shifts.size(); ++atom) {
        const auto& shift = record.shifts[atom];
        hit =
            hit || (shift.has_value() && arise::chem::carbon_type(molecule, atom) == window.type &&
                    inside(*shift, window));
    }
    const bool has = arise::chem::count_group(molecule, window.group) > 0;
    counts.hits += hit ? 1 : 0;
    counts.true_hits += hit && has ? 1 : 0;
    for (const auto& key : arise::tools::key_carbons(molecule)) {
        const auto& shift = record.shifts[key.atom];
        if (arise::chem::kAllGroups[key.group_index] == window.group &&
            key.hallmark == window.hallmark &&
            arise::chem::carbon_type(molecule, key.atom) == window.type && shift.has_value()) {
            ++counts.key_carbons;
            counts.in_window += inside(*shift, window) ? 1 : 0;
            counts.molecules.insert(record.id);
        }
    }
}

HoldOut hold_out(const std::vector<arise::tools::ShiftRecord>& records,
                 const ProtonWindow& window) {
    HoldOut counts;
    for (const auto& record : records) {
        if (record.id % 2 == 0) {
            continue;  // measured half
        }
        const auto molecule = arise::chem::parse_smiles(record.smiles);
        if (molecule.ok() && molecule.value().atoms().size() == record.shifts.size()) {
            add_record(record, molecule.value(), window, counts);
        }
    }
    return counts;
}

// Checks one round-1 window's recall and precision on the odd-ID half, and prints them.
void expect_window_holds(const std::vector<arise::tools::ShiftRecord>& records,
                         const ProtonWindow& window) {
    const auto counts = hold_out(records, window);
    const auto label = std::string{arise::chem::to_string(window.group)} + " / " + window.hallmark +
                       " / " + std::string{arise::nmr1d::to_string(window.type)};
    if (counts.key_carbons == 0 || counts.hits == 0) {
        ADD_FAILURE() << label << ": nothing to check on the odd half";
        return;
    }
    const double recall =
        static_cast<double>(counts.in_window) / static_cast<double>(counts.key_carbons);
    const double recall_mark = pass_mark(kExpectedRecall, counts.molecules.size());
    const double precision =
        static_cast<double>(counts.true_hits) / static_cast<double>(counts.hits);
    const double precision_mark = pass_mark(window.precision, counts.hits);
    std::cout << "[ holdout  ] 1H " << label << ": recall " << counts.in_window << "/"
              << counts.key_carbons << " = " << recall << " (pass mark " << recall_mark
              << "), precision " << counts.true_hits << "/" << counts.hits << " = " << precision
              << " (measured " << window.precision << ", pass mark " << precision_mark << ")\n";
    EXPECT_GE(recall, recall_mark) << label;
    EXPECT_GE(precision, precision_mark) << label;
}

TEST(GroupProtonHoldout, RoundOneWindowsHoldOnTheOddHalf) {
#ifdef __EMSCRIPTEN__
    // Skipped at run time, not compiled out, so the helpers stay used under wasm.
    GTEST_SKIP() << "measured offline; far slower under wasm checks";
#endif
    const auto windows = arise::filter::measured_proton_windows();
    ASSERT_TRUE(windows.ok()) << windows.error().message;
    const auto records = arise::tools::parse_shift_extract(
        arise::test::join_chunks(arise::test::kNmrshiftdb1hChunks));
    std::size_t checked = 0;
    for (const auto& window : windows.value()) {
        if (window.round1) {
            ++checked;
            expect_window_holds(records, window);
        }
    }
    EXPECT_EQ(checked, 2U);  // methoxy and benzene ring
}

}  // namespace

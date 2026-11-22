#include "arise/filter/carbon_fit.hpp"

#include "arise/chem/perception.hpp"
#include "arise/core/error.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <optional>
#include <tuple>
#include <utility>

namespace arise::filter {

namespace {

// Shortest fixed-point form: 29.000000 -> "29", 45.500000 -> "45.5".
std::string format_ppm(core::Ppm shift) {
    std::string text = std::to_string(shift.value());
    text.erase(text.find_last_not_of('0') + 1);
    if (text.back() == '.') {
        text.pop_back();
    }
    return text;
}

// The DEPT type of each `most` carbon environment, in order of first appearance; empty for
// an environment with no DEPT type (methane). `classes` come from chem::nmr_equivalence,
// numbered 0, 1, 2, ... below the atom count, so they index a per-class array.
std::vector<std::optional<nmr1d::CarbonType>> environment_types(
    const chem::Molecule& molecule, const std::vector<std::size_t>& classes) {
    std::vector<bool> seen(classes.size(), false);
    std::vector<std::optional<nmr1d::CarbonType>> types;
    for (std::size_t atom = 0; atom < classes.size(); ++atom) {
        if (molecule.atoms()[atom].element == chem::Element::C && !seen[classes[atom]]) {
            seen[classes[atom]] = true;
            types.push_back(chem::carbon_type(molecule, atom));
        }
    }
    return types;
}

// Maximum matching of peaks to environments (augmenting paths, peaks and environments in
// a fixed order). `owner[e]` is the peak holding environment e.
class Matching {
public:
    Matching(std::span<const nmr1d::CarbonAssignment> peaks,
             const std::vector<std::optional<nmr1d::CarbonType>>& environments, bool use_dept)
        : peaks_{peaks},
          environments_{environments},
          use_dept_{use_dept},
          owner_(environments.size()) {}

    // Whether peak `p` found an environment.
    bool place(std::size_t p) {
        std::vector<bool> visited(environments_.size(), false);
        return augment(p, visited);
    }

private:
    [[nodiscard]] bool fits(std::size_t p, std::size_t e) const {
        const auto& type = environments_[e];
        return !use_dept_ || !type.has_value() || peaks_[p].possible.contains(*type);
    }

    bool augment(std::size_t p, std::vector<bool>& visited) {
        for (std::size_t e = 0; e < environments_.size(); ++e) {
            if (visited[e] || !fits(p, e)) {
                continue;
            }
            visited[e] = true;
            const auto holder = owner_[e];
            if (!holder.has_value() || augment(*holder, visited)) {
                owner_[e] = p;
                return true;
            }
        }
        return false;
    }

    std::span<const nmr1d::CarbonAssignment> peaks_;
    std::span<const std::optional<nmr1d::CarbonType>> environments_;
    bool use_dept_;
    std::vector<std::optional<std::size_t>> owner_;
};

// The peaks in a fixed order (by shift, then by possible types), so the result does not
// depend on the order the chemist listed them in.
std::vector<nmr1d::CarbonAssignment> sorted(std::span<const nmr1d::CarbonAssignment> peaks) {
    std::vector<nmr1d::CarbonAssignment> ordered(peaks.begin(), peaks.end());
    std::stable_sort(ordered.begin(), ordered.end(),
                     [](const nmr1d::CarbonAssignment& a, const nmr1d::CarbonAssignment& b) {
                         if (a.shift != b.shift) {
                             return a.shift < b.shift;
                         }
                         return nmr1d::to_string(a.possible) < nmr1d::to_string(b.possible);
                     });
    return ordered;
}

// Per `most` carbon environment (in order of first appearance), the index of its kind, or
// std::nullopt for an environment that stands in for any kind. `kind_of(atom)` gives a
// carbon's kind; an environment whose carbons disagree, or have none, stands in for any.
template <typename KindOf>
std::vector<std::optional<std::size_t>> environment_kinds(const chem::Molecule& molecule,
                                                          const std::vector<std::size_t>& classes,
                                                          KindOf kind_of) {
    // Per class: its index in `kinds`, or kNone before its first carbon.
    constexpr auto kNone = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> index(classes.size(), kNone);
    std::vector<std::optional<std::size_t>> kinds;
    std::vector<bool> mixed;
    for (std::size_t atom = 0; atom < classes.size(); ++atom) {
        if (molecule.atoms()[atom].element != chem::Element::C) {
            continue;
        }
        const std::optional<std::size_t> kind = kind_of(atom);
        auto& e = index[classes[atom]];
        if (e == kNone) {
            e = kinds.size();
            kinds.push_back(kind);
            mixed.push_back(false);
        } else if (kinds[e] != kind) {
            mixed[e] = true;
        }
    }
    for (std::size_t e = 0; e < kinds.size(); ++e) {
        if (mixed[e]) {
            kinds[e] = std::nullopt;
        }
    }
    return kinds;
}

// Adds a BackBurner reason for each confirmed minimum (`minimums[k]`, named `names[k]`) the
// environments of `kinds` cannot supply, counting each stand-in for any one kind, and one for
// minimums that stand-ins cannot all cover together.
void check_minimums(const std::vector<std::optional<std::size_t>>& kinds,
                    std::span<const std::optional<int>> minimums,
                    std::span<const std::string_view> names, std::vector<std::string>& reasons) {
    std::vector<int> counts(minimums.size(), 0);
    int stand_ins = 0;
    for (const auto& kind : kinds) {
        if (kind) {
            ++counts[*kind];
        } else {
            ++stand_ins;
        }
    }
    int deficit = 0;
    bool short_alone = false;
    std::string confirmed;
    for (std::size_t k = 0; k < minimums.size(); ++k) {
        const auto minimum = minimums[k];
        if (!minimum) {
            continue;
        }
        const int wanted = *minimum;
        confirmed += (confirmed.empty() ? "" : ", ") + std::string{names[k]};
        deficit += std::max(0, wanted - counts[k]);
        if (counts[k] + stand_ins < wanted) {
            short_alone = true;
            reasons.push_back("confirmed at least " + std::to_string(wanted) + " " +
                              std::string{names[k]} + " peaks; the candidate has at most " +
                              std::to_string(counts[k] + stand_ins) + " " + std::string{names[k]} +
                              " carbon environments");
        }
    }
    if (!short_alone && deficit > stand_ins) {
        reasons.push_back("the confirmed " + confirmed +
                          " counts together need more carbon environments than the candidate "
                          "has");
    }
}

bool any_negative(const ConfirmedCarbonCounts& c) noexcept {
    const std::array counts{c.sp, c.sp2, c.sp3, c.c, c.ch, c.ch2, c.ch3};
    return std::any_of(counts.begin(), counts.end(),
                       [](const std::optional<int>& count) { return count.value_or(0) < 0; });
}

int tier_rank(Tier tier) noexcept {
    switch (tier) {
        case Tier::Active:
            return 0;
        case Tier::BackBurner:
            return 1;
        case Tier::Eliminated:
            return 2;
    }
    return 3;
}

core::Error no_peaks() {
    return {
        .code = core::ErrorCode::InvalidArgument,
        .message = "check_carbon_fit: the 13C list has no peaks",
    };
}

// check_carbon_fit on a non-empty peak list, given the candidate's equivalence classes
// (computed once by the caller, so the confirmed-count overload does not repeat them).
CarbonFit fit_carbons(const chem::Molecule& molecule,
                      std::span<const nmr1d::CarbonAssignment> peaks,
                      const chem::NmrEquivalence& equivalence, SuggestOptions options) {
    const auto environments = environment_types(molecule, equivalence.most);
    CarbonFit fit;
    fit.environments_fewest = chem::carbon_environments(molecule, equivalence.fewest);
    fit.environments_most = static_cast<int>(environments.size());
    const auto count = static_cast<int>(peaks.size());
    fit.assumed_overlaps = std::max(0, fit.environments_fewest - count);
    if (fit.assumed_overlaps > 0) {
        fit.notes.push_back(
            "assumes " + std::to_string(fit.assumed_overlaps) +
            (fit.assumed_overlaps == 1 ? " overlapping peak" : " overlapping peaks"));
    }
    if (count > fit.environments_most) {
        fit.verdict.tier = Tier::BackBurner;
        fit.verdict.reasons.push_back(
            "13C list has " + std::to_string(count) + " peaks; the candidate has at most " +
            std::to_string(fit.environments_most) + " carbon environments");
        return fit;
    }
    const auto ordered = sorted(peaks);
    Matching matching{ordered, environments, options.use_dept};
    for (std::size_t p = 0; p < ordered.size(); ++p) {
        if (!matching.place(p)) {
            fit.verdict.reasons.push_back(
                "no carbon environment of type " + nmr1d::to_string(ordered[p].possible) +
                " is left for the peak at " + format_ppm(ordered[p].shift) + " ppm");
        }
    }
    if (!fit.verdict.reasons.empty()) {
        fit.verdict.tier = Tier::BackBurner;
    }
    return fit;
}

}  // namespace

core::Result<CarbonFit> check_carbon_fit(const chem::Molecule& molecule,
                                         std::span<const nmr1d::CarbonAssignment> peaks,
                                         SuggestOptions options) {
    if (peaks.empty()) {
        return no_peaks();
    }
    return fit_carbons(molecule, peaks, chem::nmr_equivalence(molecule), options);
}

core::Result<CarbonFit> check_carbon_fit(const chem::Molecule& molecule,
                                         std::span<const nmr1d::CarbonAssignment> peaks,
                                         const ConfirmedCarbonCounts& confirmed,
                                         SuggestOptions options) {
    if (any_negative(confirmed)) {
        return core::Error{
            .code = core::ErrorCode::InvalidArgument,
            .message = "check_carbon_fit: a confirmed carbon count is negative",
        };
    }
    if (peaks.empty()) {
        return no_peaks();
    }
    const auto equivalence = chem::nmr_equivalence(molecule);
    auto fit = fit_carbons(molecule, peaks, equivalence, options);
    const auto& classes = equivalence.most;

    // Hybridization kinds: 0 sp, 1 sp2, 2 sp3.
    const auto hybridizations = environment_kinds(
        molecule, classes, [&molecule](std::size_t atom) -> std::optional<std::size_t> {
            const auto h = chem::carbon_hybridization(molecule, atom);
            if (!h) {
                return std::nullopt;
            }
            return static_cast<std::size_t>(*h) - 1;
        });
    const std::array hybridization_minimums{confirmed.sp, confirmed.sp2, confirmed.sp3};
    constexpr std::array<std::string_view, 3> kHybridizationNames{"sp", "sp2", "sp3"};
    check_minimums(hybridizations, hybridization_minimums, kHybridizationNames,
                   fit.verdict.reasons);

    // DEPT kinds: 0 C, 1 CH, 2 CH2, 3 CH3.
    const auto types = environment_kinds(
        molecule, classes, [&molecule](std::size_t atom) -> std::optional<std::size_t> {
            const auto type = chem::carbon_type(molecule, atom);
            if (!type) {
                return std::nullopt;
            }
            switch (*type) {
                case nmr1d::CarbonType::C:
                    return 0;
                case nmr1d::CarbonType::CH:
                    return 1;
                case nmr1d::CarbonType::CH2:
                    return 2;
                case nmr1d::CarbonType::CH3:
                    return 3;
            }
            return std::nullopt;
        });
    const std::array type_minimums{confirmed.c, confirmed.ch, confirmed.ch2, confirmed.ch3};
    constexpr std::array<std::string_view, 4> kTypeNames{"C", "CH", "CH2", "CH3"};
    check_minimums(types, type_minimums, kTypeNames, fit.verdict.reasons);

    if (!fit.verdict.reasons.empty()) {
        fit.verdict.tier = Tier::BackBurner;
    }
    return fit;
}

std::vector<std::size_t> order_by_fit(std::span<const CarbonFit> fits) {
    std::vector<std::size_t> order(fits.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(), [&fits](std::size_t a, std::size_t b) {
        const auto key = [&fits](std::size_t i) {
            return std::tuple{tier_rank(fits[i].verdict.tier), fits[i].assumed_overlaps,
                              fits[i].verdict.reasons.size()};
        };
        return key(a) < key(b);
    });
    return order;
}

}  // namespace arise::filter

#include "arise/filter/carbon_fit.hpp"

#include "arise/chem/perception.hpp"
#include "arise/core/error.hpp"

#include <algorithm>
#include <map>
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
// an environment with no DEPT type (methane).
std::vector<std::optional<nmr1d::CarbonType>> environment_types(
    const chem::Molecule& molecule, const std::vector<std::size_t>& classes) {
    std::map<std::size_t, std::size_t> seen;
    std::vector<std::optional<nmr1d::CarbonType>> types;
    for (std::size_t atom = 0; atom < classes.size(); ++atom) {
        if (molecule.atoms()[atom].element == chem::Element::C &&
            seen.emplace(classes[atom], types.size()).second) {
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

}  // namespace

core::Result<CarbonFit> check_carbon_fit(const chem::Molecule& molecule,
                                         std::span<const nmr1d::CarbonAssignment> peaks,
                                         SuggestOptions options) {
    if (peaks.empty()) {
        return core::Error{
            .code = core::ErrorCode::InvalidArgument,
            .message = "check_carbon_fit: the 13C list has no peaks",
        };
    }
    const auto equivalence = chem::nmr_equivalence(molecule);
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

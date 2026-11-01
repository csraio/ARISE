#include "arise/chem/perception.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <utility>

namespace arise::chem {
namespace {

constexpr int kAromaticBondCode = 4;  // distinct from bond orders 1-3

using Ring = std::array<std::size_t, 6>;

bool is_ring_element(Element element) noexcept {
    return element == Element::C || element == Element::N;
}

// Index into molecule.bonds() of the bond between a and b, if any.
std::optional<std::size_t> find_bond(const Molecule& molecule, std::size_t a,
                                     std::size_t b) noexcept {
    for (const auto bond : molecule.bonds_of(a)) {
        if (molecule.neighbor(bond, a) == b) {
            return bond;
        }
    }
    return std::nullopt;
}

BondOrder order_of(const Molecule& molecule, std::size_t bond) noexcept {
    return molecule.bonds()[bond].order;
}

// Every simple 6-cycle through C/N atoms, each listed once in walking order, starting
// at its smallest atom.
class SixRingFinder {
public:
    explicit SixRingFinder(const Molecule& molecule) : molecule_{&molecule} {}

    std::vector<Ring> find() {
        for (std::size_t start = 0; start < molecule_->atoms().size(); ++start) {
            if (is_ring_element(molecule_->atoms()[start].element)) {
                path_ = {start};
                extend(start);
            }
        }
        return std::move(rings_);
    }

private:
    void extend(std::size_t start) {
        const std::size_t last = path_.back();
        for (const auto bond : molecule_->bonds_of(last)) {
            const std::size_t next = molecule_->neighbor(bond, last);
            if (path_.size() == 6) {
                if (next == start) {
                    record();
                }
                continue;
            }
            if (next <= start || !is_ring_element(molecule_->atoms()[next].element) ||
                std::find(path_.begin(), path_.end(), next) != path_.end()) {
                continue;
            }
            path_.push_back(next);
            extend(start);
            path_.pop_back();
        }
    }

    // Each cycle is walked in both directions; keep the one whose second atom is smaller.
    void record() {
        if (path_[1] > path_[5]) {
            return;
        }
        Ring ring{};
        std::copy(path_.begin(), path_.end(), ring.begin());
        rings_.push_back(ring);
    }

    const Molecule* molecule_;  // not null
    std::vector<std::size_t> path_;
    std::vector<Ring> rings_;
};

// The single double-bond partner of `atom`, or nothing when the atom has no double bond,
// several, or a triple bond.
std::optional<std::size_t> sole_double_partner(const Molecule& molecule, std::size_t atom) {
    std::optional<std::size_t> partner;
    for (const auto bond : molecule.bonds_of(atom)) {
        const auto order = order_of(molecule, bond);
        if (order == BondOrder::Triple || (order == BondOrder::Double && partner.has_value())) {
            return std::nullopt;
        }
        if (order == BondOrder::Double) {
            partner = molecule.neighbor(bond, atom);
        }
    }
    return partner;
}

// A ring is aromatic when every ring atom has exactly one double bond, and that bond
// either lies in the ring or belongs to an aromatic ring already found (fused systems).
bool is_aromatic(const Molecule& molecule, const Ring& ring,
                 const std::set<std::size_t>& aromatic_bonds) {
    for (std::size_t i = 0; i < 6; ++i) {
        const std::size_t atom = ring[i];
        const auto partner = sole_double_partner(molecule, atom);
        if (!partner.has_value()) {
            return false;
        }
        const bool in_ring = *partner == ring[(i + 1) % 6] || *partner == ring[(i + 5) % 6];
        const auto bond = find_bond(molecule, atom, *partner);
        const bool in_aromatic_ring = bond.has_value() && aromatic_bonds.contains(*bond);
        if (!in_ring && !in_aromatic_ring) {
            return false;
        }
    }
    return true;
}

struct Aromaticity {
    std::vector<Ring> rings;      // in walking order
    std::set<std::size_t> bonds;  // bond indices inside aromatic rings
    std::vector<bool> atoms;      // per atom: in an aromatic ring
};

Aromaticity perceive_aromaticity(const Molecule& molecule) {
    Aromaticity result;
    result.atoms.assign(molecule.atoms().size(), false);
    auto pending = SixRingFinder{molecule}.find();
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto it = pending.begin(); it != pending.end();) {
            if (!is_aromatic(molecule, *it, result.bonds)) {
                ++it;
                continue;
            }
            for (std::size_t i = 0; i < 6; ++i) {
                if (const auto bond = find_bond(molecule, (*it)[i], (*it)[(i + 1) % 6])) {
                    result.bonds.insert(*bond);
                }
                result.atoms[(*it)[i]] = true;
            }
            result.rings.push_back(*it);
            it = pending.erase(it);
            changed = true;
        }
    }
    return result;
}

// Dense ranks of `signatures`: equal signatures get equal ranks, ordered by value, so
// the ranks do not depend on atom order.
std::vector<int> rank(const std::vector<std::vector<int>>& signatures) {
    std::map<std::vector<int>, int> ranks;
    for (const auto& signature : signatures) {
        ranks.emplace(signature, 0);
    }
    int next = 0;
    for (auto& [signature, value] : ranks) {
        value = next++;
    }
    std::vector<int> result;
    result.reserve(signatures.size());
    for (const auto& signature : signatures) {
        result.push_back(ranks.at(signature));
    }
    return result;
}

std::size_t distinct(const std::vector<int>& labels) {
    return std::set<int>(labels.begin(), labels.end()).size();
}

bool is_carbonyl_carbon(const Molecule& molecule, std::size_t atom) {
    if (molecule.atoms()[atom].element != Element::C) {
        return false;
    }
    const auto bonds = molecule.bonds_of(atom);
    return std::any_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        return order_of(molecule, bond) == BondOrder::Double &&
               molecule.atoms()[molecule.neighbor(bond, atom)].element == Element::O;
    });
}

bool has_only_single_bonds(const Molecule& molecule, std::size_t atom) {
    const auto bonds = molecule.bonds_of(atom);
    return std::all_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        return order_of(molecule, bond) == BondOrder::Single;
    });
}

// Neighbours of `atom` reached by single bonds.
std::vector<std::size_t> single_neighbors(const Molecule& molecule, std::size_t atom) {
    std::vector<std::size_t> result;
    for (const auto bond : molecule.bonds_of(atom)) {
        if (order_of(molecule, bond) == BondOrder::Single) {
            result.push_back(molecule.neighbor(bond, atom));
        }
    }
    return result;
}

// Counts the carbonyl group on carbon `atom` (which is a carbonyl carbon).
bool carbonyl_matches(const Molecule& molecule, std::size_t atom, Group group) {
    const auto& atoms = molecule.atoms();
    const auto singles = single_neighbors(molecule, atom);
    const auto has_single = [&](auto predicate) {
        return std::any_of(singles.begin(), singles.end(), predicate);
    };
    const auto is_element = [&](std::size_t n, Element e) { return atoms[n].element == e; };
    switch (group) {
        case Group::CarboxylicAcid:
            return has_single([&](std::size_t n) {
                return is_element(n, Element::O) && molecule.hydrogen_count(n) > 0;
            });
        case Group::Aldehyde:
            return molecule.hydrogen_count(atom) > 0 && !has_single([&](std::size_t n) {
                       return is_element(n, Element::O) || is_element(n, Element::N) ||
                              is_element(n, Element::S);
                   });
        case Group::Ketone:
            return molecule.bonds_of(atom).size() == 3 && singles.size() == 2 &&
                   is_element(singles[0], Element::C) && is_element(singles[1], Element::C);
        case Group::Ester:
            return has_single([&](std::size_t o) {
                if (!is_element(o, Element::O) || molecule.hydrogen_count(o) > 0) {
                    return false;
                }
                const auto others = single_neighbors(molecule, o);
                return std::any_of(others.begin(), others.end(), [&](std::size_t c) {
                    return c != atom && is_element(c, Element::C) &&
                           !is_carbonyl_carbon(molecule, c);
                });
            });
        case Group::Amide:
            return has_single([&](std::size_t n) { return is_element(n, Element::N); });
        default:
            return false;
    }
}

bool oxygen_matches(const Molecule& molecule, std::size_t atom, Group group,
                    const Aromaticity& aromaticity) {
    const auto& atoms = molecule.atoms();
    const auto neighbors = single_neighbors(molecule, atom);
    const bool hydroxyl = molecule.hydrogen_count(atom) > 0;
    const auto heavy = std::count_if(neighbors.begin(), neighbors.end(),
                                     [&](std::size_t n) { return atoms[n].element != Element::H; });
    if (molecule.bonds_of(atom).size() != neighbors.size()) {
        return false;  // a C=O oxygen
    }
    switch (group) {
        case Group::Alcohol:
        case Group::Phenol: {
            if (!hydroxyl || heavy != 1) {
                return false;
            }
            const auto carbon =
                *std::find_if(neighbors.begin(), neighbors.end(),
                              [&](std::size_t n) { return atoms[n].element != Element::H; });
            if (atoms[carbon].element != Element::C) {
                return false;
            }
            if (group == Group::Phenol) {
                return aromaticity.atoms[carbon];
            }
            return has_only_single_bonds(molecule, carbon);
        }
        case Group::Ether:
            return !hydroxyl && neighbors.size() == 2 &&
                   std::all_of(neighbors.begin(), neighbors.end(), [&](std::size_t n) {
                       return atoms[n].element == Element::C && !is_carbonyl_carbon(molecule, n);
                   });
        default:
            return false;
    }
}

bool amine_nitrogen(const Molecule& molecule, std::size_t atom) {
    if (!has_only_single_bonds(molecule, atom)) {
        return false;
    }
    const auto neighbors = single_neighbors(molecule, atom);
    return std::none_of(neighbors.begin(), neighbors.end(),
                        [&](std::size_t n) { return is_carbonyl_carbon(molecule, n); });
}

int count_bonds(const Molecule& molecule, BondOrder order, Element a, Element b,
                const std::set<std::size_t>& excluded) {
    int count = 0;
    const auto& atoms = molecule.atoms();
    for (std::size_t i = 0; i < molecule.bonds().size(); ++i) {
        const auto& bond = molecule.bonds()[i];
        const auto ea = atoms[bond.first].element;
        const auto eb = atoms[bond.second].element;
        const bool matches = (ea == a && eb == b) || (ea == b && eb == a);
        if (bond.order == order && matches && !excluded.contains(i)) {
            ++count;
        }
    }
    return count;
}

// ---- Symmetry refinement ----------------------------------------------------------------

bool is_heavy(const Molecule& molecule, std::size_t atom) {
    return molecule.atoms()[atom].element != Element::H;
}

// Starting invariants: element, isotope, charge, total hydrogens, heavy neighbours.
// Hydrogen atoms enter only through hydrogen_count, never as neighbours, so a methyl
// written with explicit [H] atoms matches one with implicit hydrogens.
std::vector<int> initial_labels(const Molecule& molecule) {
    const auto& atoms = molecule.atoms();
    std::vector<std::vector<int>> signatures(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto bonds = molecule.bonds_of(i);
        const auto heavy_degree = std::count_if(bonds.begin(), bonds.end(), [&](std::size_t b) {
            return is_heavy(molecule, molecule.neighbor(b, i));
        });
        signatures[i] = {
            static_cast<int>(atoms[i].element),
            atoms[i].isotope,
            atoms[i].charge,
            molecule.hydrogen_count(i),
            static_cast<int>(heavy_degree),
        };
    }
    return rank(signatures);
}

// Refines `labels` on neighbour labels until the number of classes stops growing. Bonds
// of aromatic rings count as one kind, so Kekule benzene stays symmetric.
std::vector<int> refine(const Molecule& molecule, const Aromaticity& aromaticity,
                        std::vector<int> labels) {
    const std::size_t n = molecule.atoms().size();
    std::vector<std::vector<int>> signatures(n);
    for (std::size_t round = 0; round <= n; ++round) {
        for (std::size_t i = 0; i < n; ++i) {
            std::vector<std::pair<int, int>> around;
            for (const auto bond : molecule.bonds_of(i)) {
                const auto other = molecule.neighbor(bond, i);
                if (!is_heavy(molecule, other)) {
                    continue;
                }
                const int code = aromaticity.bonds.contains(bond)
                                     ? kAromaticBondCode
                                     : static_cast<int>(order_of(molecule, bond));
                around.emplace_back(code, labels[other]);
            }
            std::sort(around.begin(), around.end());
            signatures[i] = {labels[i]};
            for (const auto& [code, label] : around) {
                signatures[i].push_back(code);
                signatures[i].push_back(label);
            }
        }
        auto refined = rank(signatures);
        const bool stable = distinct(refined) == distinct(labels);
        labels = std::move(refined);
        if (stable) {
            break;
        }
    }
    return labels;
}

std::vector<std::size_t> first_appearance(const std::vector<int>& labels) {
    std::map<int, std::size_t> first_seen;
    std::vector<std::size_t> classes;
    classes.reserve(labels.size());
    for (const int label : labels) {
        const auto [it, inserted] = first_seen.emplace(label, first_seen.size());
        classes.push_back(it->second);
    }
    return classes;
}

// ---- NMR equivalence rules ----------------------------------------------------------------

// An atom to separate from its class. Picks of one class normally share their new label,
// so symmetric copies stay paired; `own_label` gives the atom a label of its own.
struct Pick {
    std::size_t atom;
    bool own_label = false;
};

// Heavy neighbours of `atom`, except `skip`.
std::vector<std::size_t> heavy_neighbors(const Molecule& molecule, std::size_t atom,
                                         std::optional<std::size_t> skip = std::nullopt) {
    std::vector<std::size_t> result;
    for (const auto bond : molecule.bonds_of(atom)) {
        const auto other = molecule.neighbor(bond, atom);
        if (is_heavy(molecule, other) && other != skip) {
            result.push_back(other);
        }
    }
    return result;
}

bool has_double_bond_to(const Molecule& molecule, std::size_t atom, Element element) {
    const auto bonds = molecule.bonds_of(atom);
    return std::any_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        return order_of(molecule, bond) == BondOrder::Double &&
               molecule.atoms()[molecule.neighbor(bond, atom)].element == element;
    });
}

// The lower-index member of a pair of identical groups, if `groups` is exactly such a pair.
std::optional<std::size_t> identical_pair(const std::vector<std::size_t>& groups,
                                          const std::vector<int>& labels) {
    if (groups.size() != 2 || labels[groups[0]] != labels[groups[1]]) {
        return std::nullopt;
    }
    return std::min(groups[0], groups[1]);
}

// Whether the substituents of `atom` other than `skip` are two *different* groups. A
// nitrogen's lone pair counts as a group, so a C=N nitrogen always qualifies.
bool two_different_groups(const Molecule& molecule, std::size_t atom, std::size_t skip,
                          const std::vector<int>& labels) {
    const auto heavy = heavy_neighbors(molecule, atom, skip);
    const int hydrogens = molecule.hydrogen_count(atom);
    if (molecule.atoms()[atom].element == Element::N) {
        return heavy.size() + static_cast<std::size_t>(hydrogens) == 1;
    }
    if (heavy.size() + static_cast<std::size_t>(hydrogens) != 2) {
        return false;
    }
    if (heavy.size() == 2) {
        return labels[heavy[0]] != labels[heavy[1]];
    }
    return heavy.size() == 1;  // one group and one H; two H are identical
}

// E/Z: identical groups on one end of a non-aromatic C=C or C=N, different on the other.
std::vector<Pick> ez_splits(const Molecule& molecule, const Aromaticity& aromaticity,
                            const std::vector<int>& labels) {
    std::vector<Pick> picks;
    const auto& atoms = molecule.atoms();
    for (std::size_t b = 0; b < molecule.bonds().size(); ++b) {
        const auto& bond = molecule.bonds()[b];
        if (bond.order != BondOrder::Double || aromaticity.bonds.contains(b)) {
            continue;
        }
        for (const auto& [end, other] :
             {std::pair{bond.first, bond.second}, std::pair{bond.second, bond.first}}) {
            const auto other_element = atoms[other].element;
            if (atoms[end].element != Element::C ||
                (other_element != Element::C && other_element != Element::N)) {
                continue;
            }
            const auto pair = identical_pair(heavy_neighbors(molecule, end, other), labels);
            if (pair.has_value() && molecule.hydrogen_count(end) == 0 &&
                two_different_groups(molecule, other, end, labels)) {
                picks.push_back(Pick{.atom = *pair});
            }
        }
    }
    return picks;
}

enum class Barrier : std::uint8_t { Slow, Borderline };

// The C=O/C=S (or vinylogous) partner that slows rotation about N's bond to it, with how
// slow: amides and thioamides are slow; carbamates, ureas and enaminones borderline.
std::optional<std::pair<std::size_t, Barrier>> rotation_partner(const Molecule& molecule,
                                                                std::size_t nitrogen) {
    const auto& atoms = molecule.atoms();
    for (const auto k : heavy_neighbors(molecule, nitrogen)) {
        if (atoms[k].element != Element::C) {
            continue;
        }
        if (has_double_bond_to(molecule, k, Element::O) ||
            has_double_bond_to(molecule, k, Element::S)) {
            const auto others = heavy_neighbors(molecule, k, nitrogen);
            const bool hetero = std::any_of(others.begin(), others.end(), [&](std::size_t o) {
                const auto e = atoms[o].element;
                return (e == Element::N || e == Element::O) &&
                       find_bond(molecule, k, o).has_value() &&
                       order_of(molecule, *find_bond(molecule, k, o)) == BondOrder::Single;
            });
            return std::pair{k, hetero ? Barrier::Borderline : Barrier::Slow};
        }
        // Enaminone: N-C=C-C(=O).
        for (const auto bond : molecule.bonds_of(k)) {
            const auto c2 = molecule.neighbor(bond, k);
            if (order_of(molecule, bond) != BondOrder::Double || atoms[c2].element != Element::C) {
                continue;
            }
            for (const auto c3 : heavy_neighbors(molecule, c2, k)) {
                if (atoms[c3].element == Element::C &&
                    has_double_bond_to(molecule, c3, Element::O)) {
                    return std::pair{k, Barrier::Borderline};
                }
            }
        }
    }
    return std::nullopt;
}

// Identical groups on an amide-type nitrogen, for the given barrier class.
std::vector<Pick> rotation_splits(const Molecule& molecule, const std::vector<int>& labels,
                                  Barrier barrier) {
    std::vector<Pick> picks;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        if (molecule.atoms()[i].element != Element::N || !has_only_single_bonds(molecule, i) ||
            molecule.hydrogen_count(i) != 0) {
            continue;
        }
        const auto partner = rotation_partner(molecule, i);
        if (!partner.has_value() || partner->second != barrier) {
            continue;
        }
        if (const auto pair =
                identical_pair(heavy_neighbors(molecule, i, partner->first), labels)) {
            picks.push_back(Pick{.atom = *pair});
        }
    }
    return picks;
}

// Whether all substituent groups of `atom` (heavy neighbours by label, plus each H) differ.
bool all_groups_differ(const Molecule& molecule, std::size_t atom, const std::vector<int>& labels) {
    const auto heavy = heavy_neighbors(molecule, atom);
    std::set<int> seen;
    for (const auto n : heavy) {
        if (!seen.insert(labels[n]).second) {
            return false;
        }
    }
    return molecule.hydrogen_count(atom) <= 1;
}

bool is_stereocentre(const Molecule& molecule, std::size_t atom, const std::vector<int>& labels) {
    const auto element = molecule.atoms()[atom].element;
    const auto groups = heavy_neighbors(molecule, atom).size() +
                        static_cast<std::size_t>(molecule.hydrogen_count(atom));
    if (element == Element::C) {
        return has_only_single_bonds(molecule, atom) && groups == 4 &&
               all_groups_differ(molecule, atom, labels);
    }
    if (element == Element::S) {  // sulfoxide: S(=O) with two different groups
        return groups == 3 && has_double_bond_to(molecule, atom, Element::O) &&
               molecule.bonds_of(atom).size() == 3 && all_groups_differ(molecule, atom, labels);
    }
    if (element == Element::P) {  // three different groups, with or without P=O/P=S
        return groups >= 3 && all_groups_differ(molecule, atom, labels);
    }
    return false;
}

// The lower-index member of the identical pair on a prochiral sp3 carbon: one with two
// identical groups and two other, different ones (Caabc, at most one H).
std::optional<std::size_t> prochiral_pair(const Molecule& molecule, std::size_t atom,
                                          const std::vector<int>& labels) {
    if (molecule.atoms()[atom].element != Element::C || !has_only_single_bonds(molecule, atom)) {
        return std::nullopt;
    }
    const auto heavy = heavy_neighbors(molecule, atom);
    const int hydrogens = molecule.hydrogen_count(atom);
    if (heavy.size() + static_cast<std::size_t>(hydrogens) != 4 || hydrogens > 1) {
        return std::nullopt;
    }
    std::map<int, std::vector<std::size_t>> by_label;
    for (const auto h : heavy) {
        by_label[labels[h]].push_back(h);
    }
    std::optional<std::size_t> pair;
    for (const auto& [label, members] : by_label) {
        if (members.size() > 2 || (members.size() == 2 && pair.has_value())) {
            return std::nullopt;  // three identical groups, or two identical pairs
        }
        if (members.size() == 2) {
            pair = std::min(members[0], members[1]);
        }
    }
    const bool rest_differ = by_label.size() + static_cast<std::size_t>(hydrogens) == 3;
    return rest_differ ? pair : std::nullopt;
}

// Diastereotopic splits, which may or may not show: when the molecule has a stereocentre
// or at least two prochiral centres (these make each other's groups diastereotopic), one
// of each prochiral pair, and every stereocentre on its own (unknown configuration can
// break a symmetry). Over-splitting only widens the range, the safe direction.
std::vector<Pick> stereo_splits(const Molecule& molecule, const std::vector<int>& labels) {
    std::vector<Pick> stereocentres;
    std::vector<Pick> prochiral;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        if (is_stereocentre(molecule, i, labels)) {
            stereocentres.push_back(Pick{.atom = i, .own_label = true});
        } else if (const auto pair = prochiral_pair(molecule, i, labels)) {
            prochiral.push_back(Pick{.atom = *pair});
        }
    }
    if (stereocentres.empty() && prochiral.size() < 2) {
        return {};
    }
    stereocentres.insert(stereocentres.end(), prochiral.begin(), prochiral.end());
    return stereocentres;
}

// Gives each class of picked atoms one new label and refines, repeating until the rules
// pick nothing new. Picked atoms of one class share a label, so symmetric copies (both
// "cis" methyls of a symmetric molecule) stay equivalent to each other.
template <typename Rules>
std::vector<int> split_until_stable(const Molecule& molecule, const Aromaticity& aromaticity,
                                    std::vector<int> labels, Rules rules) {
    for (std::size_t round = 0; round <= molecule.atoms().size(); ++round) {
        const auto picks = rules(labels);
        if (picks.empty()) {
            break;
        }
        const std::size_t before = distinct(labels);
        int next = *std::max_element(labels.begin(), labels.end()) + 1;
        std::map<int, int> new_label;
        for (const auto& pick : picks) {
            if (pick.own_label) {
                labels[pick.atom] = next++;
                continue;
            }
            const auto [it, inserted] = new_label.emplace(labels[pick.atom], next);
            next += inserted ? 1 : 0;
            labels[pick.atom] = it->second;
        }
        labels = refine(molecule, aromaticity, std::move(labels));
        if (distinct(labels) == before) {
            break;
        }
    }
    return labels;
}

}  // namespace

std::optional<nmr1d::CarbonType> carbon_type(const Molecule& molecule, std::size_t atom) noexcept {
    if (molecule.atoms()[atom].element != Element::C) {
        return std::nullopt;
    }
    switch (molecule.hydrogen_count(atom)) {
        case 0:
            return nmr1d::CarbonType::C;
        case 1:
            return nmr1d::CarbonType::CH;
        case 2:
            return nmr1d::CarbonType::CH2;
        case 3:
            return nmr1d::CarbonType::CH3;
        default:
            return std::nullopt;
    }
}

std::vector<std::array<std::size_t, 6>> aromatic_rings(const Molecule& molecule) {
    auto rings = perceive_aromaticity(molecule).rings;
    for (auto& ring : rings) {
        std::sort(ring.begin(), ring.end());
    }
    std::sort(rings.begin(), rings.end());
    return rings;
}

std::vector<std::size_t> symmetry_classes(const Molecule& molecule) {
    const auto aromaticity = perceive_aromaticity(molecule);
    return first_appearance(refine(molecule, aromaticity, initial_labels(molecule)));
}

NmrEquivalence nmr_equivalence(const Molecule& molecule) {
    const auto aromaticity = perceive_aromaticity(molecule);
    const auto topological = refine(molecule, aromaticity, initial_labels(molecule));
    const auto hard = [&](const std::vector<int>& labels) {
        auto picks = ez_splits(molecule, aromaticity, labels);
        auto rotation = rotation_splits(molecule, labels, Barrier::Slow);
        picks.insert(picks.end(), rotation.begin(), rotation.end());
        return picks;
    };
    const auto all = [&](const std::vector<int>& labels) {
        auto picks = hard(labels);
        auto borderline = rotation_splits(molecule, labels, Barrier::Borderline);
        picks.insert(picks.end(), borderline.begin(), borderline.end());
        auto stereo = stereo_splits(molecule, labels);
        picks.insert(picks.end(), stereo.begin(), stereo.end());
        return picks;
    };
    const auto fewest = split_until_stable(molecule, aromaticity, topological, hard);
    const auto most = split_until_stable(molecule, aromaticity, fewest, all);
    return NmrEquivalence{.fewest = first_appearance(fewest), .most = first_appearance(most)};
}

int carbon_environments(const Molecule& molecule, const std::vector<std::size_t>& classes) {
    std::set<std::size_t> carbon_classes;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        if (molecule.atoms()[i].element == Element::C) {
            carbon_classes.insert(classes[i]);
        }
    }
    return static_cast<int>(carbon_classes.size());
}

std::string_view to_string(Group group) noexcept {
    switch (group) {
        case Group::CarboxylicAcid:
            return "carboxylic acid";
        case Group::Aldehyde:
            return "aldehyde";
        case Group::Ketone:
            return "ketone";
        case Group::Ester:
            return "ester";
        case Group::Amide:
            return "amide";
        case Group::Alcohol:
            return "alcohol";
        case Group::Phenol:
            return "phenol";
        case Group::Ether:
            return "ether";
        case Group::Amine:
            return "amine";
        case Group::Nitrile:
            return "nitrile";
        case Group::Alkene:
            return "alkene";
        case Group::BenzeneRing:
            return "benzene ring";
    }
    return "?";
}

int count_group(const Molecule& molecule, Group group) {
    const auto aromaticity = perceive_aromaticity(molecule);
    const auto& atoms = molecule.atoms();
    switch (group) {
        case Group::Nitrile:
            return count_bonds(molecule, BondOrder::Triple, Element::C, Element::N, {});
        case Group::Alkene:
            return count_bonds(molecule, BondOrder::Double, Element::C, Element::C,
                               aromaticity.bonds);
        case Group::BenzeneRing:
            return static_cast<int>(std::count_if(
                aromaticity.rings.begin(), aromaticity.rings.end(), [&](const Ring& ring) {
                    return std::all_of(ring.begin(), ring.end(), [&](std::size_t a) {
                        return atoms[a].element == Element::C;
                    });
                }));
        default:
            break;
    }
    int count = 0;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto element = atoms[i].element;
        if (element == Element::C && is_carbonyl_carbon(molecule, i)) {
            count += carbonyl_matches(molecule, i, group) ? 1 : 0;
        } else if (element == Element::O) {
            count += oxygen_matches(molecule, i, group, aromaticity) ? 1 : 0;
        } else if (element == Element::N && group == Group::Amine) {
            count += amine_nitrogen(molecule, i) ? 1 : 0;
        }
    }
    return count;
}

}  // namespace arise::chem

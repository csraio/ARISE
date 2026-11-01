// ARISE's chemistry checked against independent references (plan 3a-ref):
// - RDKit, a separate implementation, over every cleaned PubChem candidate of C4H8O2,
//   C8H8O2 and C10H12O2 (hydrogens, symmetry classes, RDKit's published fragments);
// - ChEBI's expert-curated chemical classes (functional groups);
// - experimentally assigned 13C spectra from nmrshiftdb2 (equivalent carbons).
// No group definition here is ARISE's own: every disagreement is either an ARISE bug or a
// named row in docs/chem_group_definitions.md ("Known differences"), signed off there.
// The reference files are built by scripts/build_reference_sets.py; see their headers.
#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using arise::chem::Group;

// The embedded reference files, joined once (they outlive every view taken into them).
const std::string& rdkit_oracle_text() {
    static const std::string text = arise::test::join_chunks(arise::test::kRdkitOracleChunks);
    return text;
}

const std::string& chebi_text() {
    static const std::string text = arise::test::join_chunks(arise::test::kChebiClassesChunks);
    return text;
}

const std::string& known_differences_text() {
    static const std::string text = arise::test::join_chunks(arise::test::kKnownDifferencesChunks);
    return text;
}

const std::string& nmrshiftdb_text() {
    static const std::string text =
        arise::test::join_chunks(arise::test::kNmrshiftdb13cSampleChunks);
    return text;
}

// Data lines of a reference file (comment lines start with '#').
std::vector<std::string_view> data_lines(std::string_view text) {
    std::vector<std::string_view> lines;
    while (!text.empty()) {
        const auto end = text.find('\n');
        const auto line = text.substr(0, end);
        text = end == std::string_view::npos ? std::string_view{} : text.substr(end + 1);
        if (!line.empty() && line.front() != '#') {
            lines.push_back(line);
        }
    }
    return lines;
}

std::vector<std::string_view> split(std::string_view text, char separator) {
    std::vector<std::string_view> parts;
    while (true) {
        const auto end = text.find(separator);
        parts.push_back(text.substr(0, end));
        if (end == std::string_view::npos) {
            return parts;
        }
        text = text.substr(end + 1);
    }
}

std::vector<int> integers(std::string_view text) {
    std::vector<int> values;
    for (const auto part : split(text, ' ')) {
        int value = 0;
        bool negative = false;
        for (const char c : part) {
            if (c == '-') {
                negative = true;
            } else {
                value = (value * 10) + (c - '0');
            }
        }
        values.push_back(negative ? -value : value);
    }
    return values;
}

// One line of pubchem_rdkit_oracle.tsv.
struct OracleEntry {
    std::string_view smiles;
    std::vector<int> hydrogens;
    std::vector<int> classes;
    std::vector<int> fragments;  // order of the file's "# fragments:" line
    bool aromaticity_differs = false;
};

std::vector<OracleEntry> oracle() {
    std::vector<OracleEntry> entries;
    for (const auto line : data_lines(rdkit_oracle_text())) {
        const auto fields = split(line, '\t');
        if (fields.size() != 5) {
            ADD_FAILURE() << "malformed oracle line: " << line;
            continue;
        }
        entries.push_back(OracleEntry{
            .smiles = fields[0],
            .hydrogens = integers(fields[1]),
            .classes = integers(fields[2]),
            .fragments = integers(fields[3]),
            .aromaticity_differs = fields[4] == "1",
        });
    }
    return entries;
}

std::string join(const std::vector<int>& values) {
    std::string text;
    for (const int value : values) {
        text += (text.empty() ? "" : " ") + std::to_string(value);
    }
    return text;
}

// Reports up to `shown` collected disagreements and fails if there are any.
void expect_no_disagreements(const std::vector<std::string>& problems, std::size_t checked) {
    constexpr std::size_t kShown = 15;
    std::string listing;
    for (std::size_t i = 0; i < problems.size() && i < kShown; ++i) {
        listing += "\n  " + problems[i];
    }
    EXPECT_TRUE(problems.empty()) << problems.size() << " of " << checked
                                  << " molecules disagree:" << listing;
}

// ---- RDKit oracle -------------------------------------------------------------------

TEST(Reference, OracleCoversThousandsOfRealCandidates) {
    // C4H8O2 122 + C8H8O2 and C10H12O2 (radical-free, cleaned): several thousand.
    EXPECT_GT(oracle().size(), 4000U);
}

TEST(Reference, EveryOracleMoleculeParses) {
    std::vector<std::string> problems;
    const auto entries = oracle();
    for (const auto& entry : entries) {
        const auto molecule = arise::chem::parse_smiles(entry.smiles);
        if (!molecule.ok()) {
            problems.push_back(std::string{entry.smiles} + ": " + molecule.error().message);
        }
    }
    expect_no_disagreements(problems, entries.size());
}

TEST(Reference, HydrogenCountsMatchRdkit) {
    std::vector<std::string> problems;
    const auto entries = oracle();
    for (const auto& entry : entries) {
        const auto molecule = arise::chem::parse_smiles(entry.smiles);
        if (!molecule.ok()) {
            continue;  // reported by EveryOracleMoleculeParses
        }
        std::vector<int> ours;
        for (std::size_t i = 0; i < molecule.value().atoms().size(); ++i) {
            ours.push_back(molecule.value().hydrogen_count(i));
        }
        if (ours != entry.hydrogens) {
            problems.push_back(std::string{entry.smiles} + ": ARISE " + join(ours) + " vs RDKit " +
                               join(entry.hydrogens));
        }
    }
    expect_no_disagreements(problems, entries.size());
}

TEST(Reference, SymmetryClassesMatchRdkit) {
    // Rings RDKit calls aromatic but ARISE's model does not (azulene, azoles with two
    // heteroatoms, quinoid tautomers, arynes) are excluded: ARISE reads their Kekule bonds
    // literally by design (documented in perception.hpp).
    std::vector<std::string> problems;
    std::size_t checked = 0;
    for (const auto& entry : oracle()) {
        const auto molecule = arise::chem::parse_smiles(entry.smiles);
        if (entry.aromaticity_differs || !molecule.ok()) {
            continue;
        }
        ++checked;
        std::vector<int> ours;
        for (const auto cls : arise::chem::symmetry_classes(molecule.value())) {
            ours.push_back(static_cast<int>(cls));
        }
        if (ours != entry.classes) {
            problems.push_back(std::string{entry.smiles} + ": ARISE " + join(ours) + " vs RDKit " +
                               join(entry.classes));
        }
    }
    expect_no_disagreements(problems, checked);
}

// Column names from a "# <prefix> a b c" header line.
std::vector<std::string_view> header_names(std::string_view text, std::string_view prefix) {
    const auto start = text.find(prefix);
    if (start == std::string_view::npos) {
        return {};
    }
    text = text.substr(start + prefix.size());
    return split(text.substr(0, text.find('\n')), ' ');
}

// Signed-off differences from tests/data/reference/known_differences.tsv (columns: source,
// molecule, check, reason; docs/chem_group_definitions.md summarises them). Each row must
// still be a disagreement, so fixed cases cannot linger in the list.
class KnownDifferences {
public:
    explicit KnownDifferences(std::string_view source) {
        for (const auto line : data_lines(known_differences_text())) {
            const auto fields = split(line, '\t');
            if (fields.size() != 4) {
                ADD_FAILURE() << "malformed known difference: " << line;
            } else if (fields[0] == source) {
                unused_.insert({std::string{fields[1]}, std::string{fields[2]}});
            }
        }
        all_ = unused_;
    }

    // Whether {molecule, check} is listed; marks it as used.
    bool contains(const std::string& molecule, std::string_view check) {
        const std::pair<std::string, std::string> key{molecule, std::string{check}};
        unused_.erase(key);
        return all_.contains(key);
    }

    // Listed differences that no longer disagree.
    [[nodiscard]] std::vector<std::string> stale() const {
        std::vector<std::string> rows;
        for (const auto& [molecule, check] : unused_) {
            rows.push_back(molecule);
            rows.back().append(" ").append(check).append(" is listed but now agrees");
        }
        return rows;
    }

private:
    std::set<std::pair<std::string, std::string>> all_;
    std::set<std::pair<std::string, std::string>> unused_;
};

// RDKit fragments whose published definition (RDKit's FragmentDescriptors.csv) matches the
// sum of these ARISE groups. Fragments without an equivalent are not compared; ChEBI covers
// those groups.
struct FragmentCheck {
    std::string_view fragment;
    std::vector<Group> groups;
};

const std::vector<FragmentCheck>& fragment_checks() {
    static const std::vector<FragmentCheck> checks{
        {.fragment = "fr_ketone", .groups = {Group::Ketone}},
        {.fragment = "fr_aldehyde", .groups = {Group::Aldehyde}},
        {.fragment = "fr_COO2", .groups = {Group::CarboxylicAcid}},
        {.fragment = "fr_ester", .groups = {Group::CarboxylicEster}},
        {.fragment = "fr_Al_OH", .groups = {Group::Alcohol, Group::Enol}},
        {.fragment = "fr_phenol", .groups = {Group::Phenol}},
        {.fragment = "fr_nitrile", .groups = {Group::Nitrile}},
        {.fragment = "fr_urea", .groups = {Group::Urea}},
        {.fragment = "fr_benzene", .groups = {Group::BenzeneRing}},
        {.fragment = "fr_pyridine", .groups = {Group::PyridineRing}},
        {.fragment = "fr_furan", .groups = {Group::FuranRing}},
        {.fragment = "fr_thiophene", .groups = {Group::ThiopheneRing}},
    };
    return checks;
}

// The fragment checks one oracle molecule fails (and that are not signed-off differences),
// as " <fragment> <ARISE>/<RDKit>" items.
std::string fragment_differences(const arise::chem::Molecule& molecule, const OracleEntry& entry,
                                 const std::vector<std::string_view>& names,
                                 KnownDifferences& known) {
    std::string differences;
    for (const auto& check : fragment_checks()) {
        const auto column = std::find(names.begin(), names.end(), check.fragment);
        if (column == names.end()) {
            ADD_FAILURE() << check.fragment << " missing from the oracle";
            continue;
        }
        const int theirs = entry.fragments[static_cast<std::size_t>(column - names.begin())];
        int ours = 0;
        for (const auto group : check.groups) {
            ours += arise::chem::count_group(molecule, group);
        }
        if (ours != theirs && !known.contains(std::string{entry.smiles}, check.fragment)) {
            differences.append(" ")
                .append(check.fragment)
                .append(" ")
                .append(std::to_string(ours))
                .append("/")
                .append(std::to_string(theirs));
        }
    }
    return differences;
}

TEST(Reference, GroupCountsMatchRdkitFragments) {
    const auto names = header_names(rdkit_oracle_text(), "# fragments: ");
    ASSERT_FALSE(names.empty()) << "oracle has no '# fragments:' line";
    KnownDifferences known{"rdkit"};
    std::vector<std::string> problems;
    std::size_t checked = 0;
    for (const auto& entry : oracle()) {
        const auto molecule = arise::chem::parse_smiles(entry.smiles);
        if (entry.aromaticity_differs || !molecule.ok()) {
            continue;
        }
        ASSERT_EQ(entry.fragments.size(), names.size()) << entry.smiles;
        ++checked;
        const auto differences = fragment_differences(molecule.value(), entry, names, known);
        if (!differences.empty()) {
            problems.push_back(std::string{entry.smiles} + " (ARISE/RDKit):" + differences);
        }
    }
    const auto stale = known.stale();
    problems.insert(problems.end(), stale.begin(), stale.end());
    expect_no_disagreements(problems, checked);
}

// ---- ChEBI curated classes -------------------------------------------------------------

// docs/chem_group_definitions.md: each ChEBI class and the ARISE groups that satisfy it.
// Where ChEBI's class contains several ARISE groups, finding any one of them counts (user
// decision 2026-10-31): ChEBI files carbamate and carbonate esters under carboxylic esters,
// and pyridones under pyridines. ChEBI's ring classes are skeleton classes that include
// non-aromatic rings, so they are checked only where RDKit sees the aromatic ring
// (`rdkit_ring`, a column of the file's "# rings:" line).
struct ChebiCheck {
    int chebi = 0;
    std::vector<Group> groups;
    std::string_view rdkit_ring;
};

const std::vector<ChebiCheck>& chebi_checks() {
    static const std::vector<ChebiCheck> checks{
        {.chebi = 33575, .groups = {Group::CarboxylicAcid}, .rdkit_ring = {}},
        {.chebi = 17478, .groups = {Group::Aldehyde}, .rdkit_ring = {}},
        {.chebi = 17087, .groups = {Group::Ketone}, .rdkit_ring = {}},
        {
            .chebi = 33308,
            .groups = {Group::CarboxylicEster, Group::CarbamateEster, Group::CarbonateEster},
            .rdkit_ring = {},
        },
        {.chebi = 46722, .groups = {Group::CarbonateEster}, .rdkit_ring = {}},
        {.chebi = 35873, .groups = {Group::CarboxylicAnhydride}, .rdkit_ring = {}},
        {.chebi = 37622, .groups = {Group::Carboxamide}, .rdkit_ring = {}},
        {.chebi = 23003, .groups = {Group::CarbamateEster}, .rdkit_ring = {}},
        {.chebi = 47857, .groups = {Group::Urea}, .rdkit_ring = {}},
        {.chebi = 30879, .groups = {Group::Alcohol}, .rdkit_ring = {}},
        {.chebi = 33823, .groups = {Group::Enol}, .rdkit_ring = {}},
        {.chebi = 33853, .groups = {Group::Phenol}, .rdkit_ring = {}},
        {.chebi = 38183, .groups = {Group::Pyridone}, .rdkit_ring = {}},
        {.chebi = 25698, .groups = {Group::Ether}, .rdkit_ring = {}},
        {.chebi = 32952, .groups = {Group::Amine}, .rdkit_ring = {}},
        {.chebi = 35358, .groups = {Group::Sulfonamide}, .rdkit_ring = {}},
        {.chebi = 24631, .groups = {Group::Hydrazine}, .rdkit_ring = {}},
        {.chebi = 24709, .groups = {Group::Hydroxylamine}, .rdkit_ring = {}},
        {.chebi = 24783, .groups = {Group::Imine}, .rdkit_ring = {}},
        {.chebi = 25750, .groups = {Group::Oxime}, .rdkit_ring = {}},
        {.chebi = 18379, .groups = {Group::Nitrile}, .rdkit_ring = {}},
        {.chebi = 78840, .groups = {Group::Alkene}, .rdkit_ring = {}},
        {.chebi = 22712, .groups = {Group::BenzeneRing}, .rdkit_ring = "fr_benzene"},
        {
            .chebi = 26421,
            .groups = {Group::PyridineRing, Group::Pyridone},
            .rdkit_ring = "fr_pyridine",
        },
        {.chebi = 24129, .groups = {Group::FuranRing}, .rdkit_ring = "fr_furan"},
        {.chebi = 26455, .groups = {Group::PyrroleRing}, .rdkit_ring = "n1cccc1"},
        {.chebi = 26961, .groups = {Group::ThiopheneRing}, .rdkit_ring = "fr_thiophene"},
    };
    return checks;
}

TEST(Reference, ChebiChecksCoverEveryGroup) {
    std::set<int> classes;
    std::set<Group> groups;
    for (const auto& check : chebi_checks()) {
        classes.insert(check.chebi);
        groups.insert(check.groups.begin(), check.groups.end());
    }
    EXPECT_EQ(classes.size(), chebi_checks().size());  // each class once
    EXPECT_EQ(groups.size(), arise::chem::kAllGroups.size());
}

// Whether RDKit sees the aromatic ring `check` needs (always true for non-ring classes).
bool rdkit_sees_ring(const ChebiCheck& check, const std::vector<std::string_view>& rings,
                     const std::vector<int>& ring_counts) {
    if (check.rdkit_ring.empty()) {
        return true;
    }
    const auto column = std::find(rings.begin(), rings.end(), check.rdkit_ring);
    EXPECT_NE(column, rings.end()) << check.rdkit_ring;
    return column != rings.end() &&
           ring_counts[static_cast<std::size_t>(column - rings.begin())] > 0;
}

// The ChEBI classes on one line of chebi_classes.tsv that ARISE does not find (and that are
// not signed-off differences), as " CHEBI:<id> (<group>)" items.
std::string missing_chebi_classes(const arise::chem::Molecule& molecule, const std::string& id,
                                  std::string_view classes, const std::vector<int>& ring_counts,
                                  const std::vector<std::string_view>& rings,
                                  KnownDifferences& known) {
    std::string missing;
    for (const int cls : integers(classes)) {
        const auto check = std::find_if(chebi_checks().begin(), chebi_checks().end(),
                                        [&](const ChebiCheck& c) { return c.chebi == cls; });
        if (check == chebi_checks().end()) {
            ADD_FAILURE() << "unknown ChEBI class " << cls;
            continue;
        }
        if (!rdkit_sees_ring(*check, rings, ring_counts)) {
            continue;  // a non-aromatic member of a ChEBI skeleton class
        }
        const bool found =
            std::any_of(check->groups.begin(), check->groups.end(),
                        [&](Group group) { return arise::chem::count_group(molecule, group) > 0; });
        if (!found && !known.contains(id, std::to_string(cls))) {
            missing.append(" CHEBI:")
                .append(std::to_string(cls))
                .append(" (")
                .append(arise::chem::to_string(check->groups.front()))
                .append(")");
        }
    }
    return missing;
}

TEST(Reference, FindsEveryGroupChebiAsserts) {
    // ChEBI's annotation is incomplete, so only asserted classes are checked: when ChEBI
    // places a molecule in a class, ARISE must find one of the class's groups.
    const auto rings = header_names(chebi_text(), "# rings: ");
    ASSERT_FALSE(rings.empty()) << "ChEBI file has no '# rings:' line";
    KnownDifferences known{"chebi"};
    std::vector<std::string> problems;
    std::size_t checked = 0;
    for (const auto line : data_lines(chebi_text())) {
        const auto fields = split(line, '\t');
        if (fields.size() != 4) {
            ADD_FAILURE() << "malformed ChEBI line: " << line;
            continue;
        }
        const std::string id{fields[0]};
        const auto ring_counts = integers(fields[3]);
        ASSERT_EQ(ring_counts.size(), rings.size()) << line;
        const auto molecule = arise::chem::parse_smiles(fields[1]);
        ++checked;
        const std::string label = "CHEBI:" + id + " " + std::string{fields[1]};
        if (!molecule.ok()) {
            problems.push_back(label + ": " + molecule.error().message);
            continue;
        }
        const auto missing =
            missing_chebi_classes(molecule.value(), id, fields[2], ring_counts, rings, known);
        if (!missing.empty()) {
            problems.push_back(label);
            problems.back().append(": ChEBI says").append(missing);
        }
    }
    const auto stale = known.stale();
    problems.insert(problems.end(), stale.begin(), stale.end());
    EXPECT_GT(checked, 1000U);
    expect_no_disagreements(problems, checked);
}

// ---- nmrshiftdb2 experimental 13C assignments ------------------------------------------

// One line of nmrshiftdb2_13c_sample.tsv: per atom, the assigned peak, -1 for a carbon
// the database does not list, -2 for a non-carbon.
struct Assigned {
    std::string_view id;
    std::string_view smiles;
    std::vector<int> peaks;
};

std::vector<Assigned> assigned_spectra() {
    std::vector<Assigned> spectra;
    for (const auto line : data_lines(nmrshiftdb_text())) {
        const auto fields = split(line, '\t');
        if (fields.size() != 3) {
            ADD_FAILURE() << "malformed nmrshiftdb2 line: " << line;
            continue;
        }
        std::vector<int> peaks;
        for (const auto mark : split(fields[2], ' ')) {
            if (mark == ".") {
                peaks.push_back(-2);  // not carbon
            } else if (mark == "-") {
                peaks.push_back(-1);  // carbon not listed
            } else {
                peaks.push_back(integers(mark).front());
            }
        }
        spectra.push_back(Assigned{.id = fields[0], .smiles = fields[1], .peaks = peaks});
    }
    return spectra;
}

TEST(Reference, EquivalentCarbonsNeverHaveDifferentExperimentalShifts) {
    // If ARISE calls two listed carbons equivalent even after every split NMR may show
    // (nmr_equivalence().most), they must sit on the same peak: equivalent nuclei have
    // identical shifts. The converse is not required: carbons can overlap by accident.
    std::vector<std::string> problems;
    const auto spectra = assigned_spectra();
    for (const auto& spectrum : spectra) {
        const auto molecule = arise::chem::parse_smiles(spectrum.smiles);
        if (!molecule.ok()) {
            problems.push_back(std::string{spectrum.id} + ": " + molecule.error().message);
            continue;
        }
        const auto classes = arise::chem::nmr_equivalence(molecule.value()).most;
        std::map<std::size_t, int> peak_of_class;
        for (std::size_t i = 0; i < spectrum.peaks.size(); ++i) {
            if (spectrum.peaks[i] < 0) {
                continue;
            }
            const auto [it, inserted] = peak_of_class.emplace(classes[i], spectrum.peaks[i]);
            if (!inserted && it->second != spectrum.peaks[i]) {
                problems.push_back(std::string{spectrum.id} + " " + std::string{spectrum.smiles} +
                                   ": atom " + std::to_string(i) +
                                   " is equivalent to an atom on peak " +
                                   std::to_string(it->second) + " but sits on peak " +
                                   std::to_string(spectrum.peaks[i]));
                break;
            }
        }
    }
    expect_no_disagreements(problems, spectra.size());
}

TEST(Reference, ReportsUnlistedCarbonsWithoutAListedEquivalent) {
    // Databases list one atom per carbon environment. A carbon that is not listed should
    // therefore be equivalent to a listed one; if ARISE gives it a class of its own,
    // either ARISE failed to see a symmetry or the database left a peak unassigned.
    std::vector<std::string> problems;
    const auto spectra = assigned_spectra();
    for (const auto& spectrum : spectra) {
        const auto molecule = arise::chem::parse_smiles(spectrum.smiles);
        if (!molecule.ok()) {
            continue;  // reported by the test above
        }
        const auto classes = arise::chem::nmr_equivalence(molecule.value()).fewest;
        std::map<std::size_t, bool> listed;
        for (std::size_t i = 0; i < spectrum.peaks.size(); ++i) {
            if (spectrum.peaks[i] >= 0) {
                listed[classes[i]] = true;
            }
        }
        for (std::size_t i = 0; i < spectrum.peaks.size(); ++i) {
            if (spectrum.peaks[i] == -1 && !listed[classes[i]]) {
                problems.push_back(std::string{spectrum.id} + " " + std::string{spectrum.smiles} +
                                   ": carbon " + std::to_string(i) + " has no listed equivalent");
                break;
            }
        }
    }
    // Report-only: unassigned peaks in the database make this noisy (plan 3a-eq).
    std::cout << "[ reference ] " << problems.size() << " of " << spectra.size()
              << " spectra have an unlisted carbon without a listed equivalent\n";
    for (std::size_t i = 0; i < problems.size() && i < 10; ++i) {
        std::cout << "  " << problems[i] << "\n";
    }
}

}  // namespace

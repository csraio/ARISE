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
#include "arise/core/formula.hpp"
#include "arise/io/carbon_csv.hpp"
#include "arise/io/number.hpp"
#include "arise/io/peak_csv.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cctype>
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

const std::string& group_definitions_text() {
    static const std::string text = arise::test::join_chunks(arise::test::kGroupDefinitionsChunks);
    return text;
}

const std::string& measurements_text() {
    static const std::string text = arise::test::join_chunks(arise::test::kMeasurementsChunks);
    return text;
}

const std::string& spec_oracle_text() {
    static const std::string text = arise::test::join_chunks(arise::test::kSpecOracleChunks);
    return text;
}

const std::string& valence_probes_text() {
    static const std::string text = arise::test::join_chunks(arise::test::kValenceProbesChunks);
    return text;
}

const std::string& dou_checks_text() {
    static const std::string text = arise::test::join_chunks(arise::test::kDouChecksChunks);
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

std::string_view trim(std::string_view text) {
    const auto first = text.find_first_not_of(' ');
    if (first == std::string_view::npos) {
        return {};
    }
    return text.substr(first, text.find_last_not_of(' ') - first + 1);
}

std::string lowercase(std::string_view text) {
    std::string lower{text};
    for (auto& c : lower) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return lower;
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

// Oracle lines of `text` (pubchem_rdkit_oracle.tsv, or the accepted rows of
// spec_rdkit_oracle.tsv; rows marked 'rejected' are skipped).
std::vector<OracleEntry> oracle(std::string_view text) {
    std::vector<OracleEntry> entries;
    for (const auto line : data_lines(text)) {
        if (line.ends_with("\trejected")) {
            continue;
        }
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

std::vector<OracleEntry> oracle() {
    return oracle(rdkit_oracle_text());
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
    constexpr std::size_t kShown = 40;
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

// The signed-off group table of docs/chem_group_definitions.md: row number (the Group's
// position in kAllGroups, from 1), group name and the "CHEBI:<id> <name>" cell.
struct DefinitionRow {
    std::string_view group;
    int chebi = 0;
    std::string_view chebi_name;
};

std::map<int, DefinitionRow> definition_rows() {
    std::map<int, DefinitionRow> rows;
    for (const auto line : split(group_definitions_text(), '\n')) {
        const auto cells = split(line, '|');
        if (cells.size() != 9) {
            continue;  // not a row of the 7-column groups table
        }
        const auto number = integers(trim(cells[1]));
        const auto chebi = trim(cells[6]);
        if (number.size() != 1 || !chebi.starts_with("CHEBI:")) {
            continue;  // header or separator row
        }
        const auto space = chebi.find(' ');
        const auto id = integers(chebi.substr(6, space - 6));
        rows[number.front()] = {
            .group = trim(cells[2]),
            .chebi = id.size() == 1 ? id.front() : -1,
            .chebi_name = chebi.substr(space + 1),
        };
    }
    return rows;
}

// ChEBI's own name of each class, as the builder read it from the ChEBI ontology
// ("# class <id>: <name>" header lines of chebi_classes.tsv).
std::map<int, std::string_view> chebi_names() {
    constexpr std::string_view kPrefix = "# class ";
    std::map<int, std::string_view> names;
    for (const auto line : split(chebi_text(), '\n')) {
        const auto colon = line.find(": ");
        if (!line.starts_with(kPrefix) || colon == std::string_view::npos) {
            continue;
        }
        const auto id = integers(line.substr(kPrefix.size(), colon - kPrefix.size()));
        if (id.size() == 1) {
            names[id.front()] = line.substr(colon + 2);
        }
    }
    return names;
}

// What is wrong with the ChEBI class of `group` (row `number` of the doc table); empty if
// nothing. The class ID is typed twice by hand (the C++ table and the doc), and the doc
// pairs it with a name, so a typo in any of them disagrees with ChEBI's ontology.
std::string chebi_class_problem(Group group, const DefinitionRow& row,
                                const std::map<int, std::string_view>& names) {
    const auto check = std::find_if(chebi_checks().begin(), chebi_checks().end(),
                                    [&](const ChebiCheck& c) { return c.groups.front() == group; });
    std::string problem;
    if (!lowercase(row.group).starts_with(arise::chem::to_string(group))) {
        problem.append(" doc row names ").append(row.group);
    }
    if (check == chebi_checks().end() || check->chebi != row.chebi) {
        problem.append(" C++ and doc ChEBI IDs differ");
    }
    const auto name = names.find(row.chebi);
    if (name == names.end() || name->second != row.chebi_name) {
        problem.append(" doc says CHEBI:")
            .append(std::to_string(row.chebi))
            .append(" is '")
            .append(row.chebi_name)
            .append("', ChEBI says '")
            .append(name == names.end() ? std::string_view{"(not listed)"} : name->second)
            .append("'");
    }
    return problem;
}

TEST(Reference, ChebiClassIdsMatchTheSignedOffDefinitionsAndChebiNames) {
    const auto rows = definition_rows();
    const auto names = chebi_names();
    ASSERT_EQ(rows.size(), arise::chem::kAllGroups.size());
    EXPECT_EQ(names.size(), chebi_checks().size());
    for (std::size_t i = 0; i < arise::chem::kAllGroups.size(); ++i) {
        const auto group = arise::chem::kAllGroups[i];
        const auto row = rows.find(static_cast<int>(i) + 1);
        ASSERT_NE(row, rows.end()) << "doc has no row " << i + 1;
        EXPECT_EQ(chebi_class_problem(group, row->second, names), "")
            << arise::chem::to_string(group);
    }
}

// Whether `text` has the shape of a standard InChIKey: 14 letters, '-', 10 letters, '-', 1.
bool is_inchikey(std::string_view text) {
    constexpr std::string_view kShape = "AAAAAAAAAAAAAA-AAAAAAAAAA-A";
    if (text.size() != kShape.size()) {
        return false;
    }
    for (std::size_t i = 0; i < text.size(); ++i) {
        const bool letter = text[i] >= 'A' && text[i] <= 'Z';
        if (kShape[i] == '-' ? text[i] != '-' : !letter) {
            return false;
        }
    }
    return true;
}

// Ring checks that are not one of RDKit's published fr_* fragments: patterns ARISE's
// checks wrote themselves.
std::set<std::string_view> unpublished_ring_patterns() {
    std::set<std::string_view> patterns;
    for (const auto& check : chebi_checks()) {
        if (!check.rdkit_ring.empty() && !check.rdkit_ring.starts_with("fr_")) {
            patterns.insert(check.rdkit_ring);
        }
    }
    return patterns;
}

// What is wrong with a 'pubchem' row (columns: source, InChI, pattern, reason naming the
// PubChem CID and InChIKey); empty if nothing.
std::string_view pubchem_row_problem(const std::vector<std::string_view>& fields,
                                     const std::set<std::string_view>& patterns) {
    if (!patterns.contains(fields[2])) {
        return "no check uses the pattern";
    }
    if (!fields[1].starts_with("InChI=1S/")) {
        return "molecule is not a standard InChI";
    }
    const auto key = fields[3].find("InChIKey ");
    if (key == std::string_view::npos || !is_inchikey(fields[3].substr(key + 9, 27))) {
        return "reason names no InChIKey";
    }
    if (fields[3].find("PubChem CID ") == std::string_view::npos) {
        return "reason names no PubChem CID";
    }
    return {};
}

TEST(Reference, RingPatternsWithoutAnRdkitFragmentArePubChemCompounds) {
    // A ring pattern ARISE's checks wrote themselves must be tied in known_differences.tsv to
    // the PubChem compound it encodes: its InChI and InChIKey, as PubChem gives them
    // (scripts/build_reference_sets.py re-checks both against PubChem and RDKit).
    const auto patterns = unpublished_ring_patterns();
    EXPECT_EQ(patterns, (std::set<std::string_view>{"n1cccc1"}));
    std::multiset<std::string_view> backed;
    for (const auto line : data_lines(known_differences_text())) {
        const auto fields = split(line, '\t');
        if (fields.size() == 4 && fields[0] == "pubchem") {
            EXPECT_EQ(pubchem_row_problem(fields, patterns), "") << line;
            backed.insert(fields[2]);
        }
    }
    EXPECT_EQ(backed, (std::multiset<std::string_view>{patterns.begin(), patterns.end()}));
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

// ---- Provenance: hand-typed values held to machine-read sources ------------------------

// Whether ARISE accepts `smiles` as a candidate structure (parses, closed shell, neutral).
bool arise_accepts(std::string_view smiles) {
    const auto molecule = arise::chem::parse_smiles(smiles);
    return molecule.ok() && !arise::chem::structure_issue(molecule.value()).has_value();
}

// Measurement rows of measurements.tsv whose first column is `kind`.
std::vector<std::vector<std::string_view>> measurements(std::string_view kind) {
    std::vector<std::vector<std::string_view>> rows;
    for (const auto line : data_lines(measurements_text())) {
        auto fields = split(line, '\t');
        if (fields.front() == kind) {
            rows.push_back(std::move(fields));
        }
    }
    return rows;
}

double decimal(std::string_view text) {
    const auto parsed = arise::io::parse_decimal(text);
    EXPECT_TRUE(parsed.ok()) << text;
    return parsed.ok() ? parsed.value() : 0.0;
}

// Whether a measured shift range ("shift-range", nucleus, low, high) lies inside ARISE's
// plausibility window for that nucleus; empty if it does.
std::string window_problem(const std::vector<std::string_view>& range) {
    const bool carbon = range[1] == "13C";
    if (!carbon && range[1] != "1H") {
        return "unknown nucleus " + std::string{range[1]};
    }
    const double low =
        carbon ? arise::io::kMinCarbonShift.value() : arise::io::kMinProtonShift.value();
    const double high =
        carbon ? arise::io::kMaxCarbonShift.value() : arise::io::kMaxProtonShift.value();
    if (decimal(range[2]) < low || decimal(range[3]) > high) {
        return std::string{range[1]} + " window " + std::to_string(low) + ".." +
               std::to_string(high) + " misses measured " + std::string{range[2]} + ".." +
               std::string{range[3]};
    }
    return {};
}

TEST(Provenance, ShiftWindowsContainEveryShiftListedInNmrshiftdb) {
    // The plausibility windows are conventions; they must at least admit every shift
    // nmrshiftdb2 lists for an in-scope molecule (measured by the builder).
    const auto ranges = measurements("shift-range");
    ASSERT_EQ(ranges.size(), 2U);
    for (const auto& range : ranges) {
        EXPECT_EQ(window_problem(range), "");
    }
}

TEST(Provenance, DocumentedOximeShiftsAreTheMeasuredOnes) {
    // docs/chem_group_definitions.md quotes the C=N measurement that decided "no oxime
    // ether group"; its table rows must be exactly the builder's output.
    const auto rows = measurements("cn-shift");
    ASSERT_EQ(rows.size(), 2U);
    for (const auto& row : rows) {
        std::string expected = "| ";
        expected.append(row[1]).append(" | ").append(row[2]).append(" | ").append(row[3]);
        expected.append(" | ").append(row[4]).append("–").append(row[5]).append(" |");
        EXPECT_NE(group_definitions_text().find(expected), std::string::npos)
            << "doc lacks the measured row: " << expected;
    }
}

// "1234567" -> "1,234,567", as the doc writes counts.
std::string with_thousands(std::size_t value) {
    std::string digits = std::to_string(value);
    for (auto i = static_cast<std::ptrdiff_t>(digits.size()) - 3; i > 0; i -= 3) {
        digits.insert(static_cast<std::size_t>(i), ",");
    }
    return digits;
}

std::string_view source_label(std::string_view source) {
    if (source == "rdkit") {
        return "RDKit";
    }
    if (source == "chebi") {
        return "ChEBI";
    }
    if (source == "spec") {
        return "chem tests";
    }
    return source == "pubchem" ? std::string_view{"PubChem"} : source;
}

// Rows per "source category" of known_differences.tsv, as the doc's table rows
// "| category | count | Source |".
std::set<std::string> known_difference_rows() {
    std::map<std::pair<std::string_view, std::string_view>, int> counts;
    for (const auto line : data_lines(known_differences_text())) {
        const auto fields = split(line, '\t');
        if (fields.size() == 4) {
            ++counts[{fields[0], fields[3].substr(0, fields[3].find(':'))}];
        }
    }
    std::set<std::string> rows;
    for (const auto& [key, count] : counts) {
        std::string row = "| ";
        row.append(key.second).append(" | ").append(std::to_string(count)).append(" | ");
        rows.insert(row.append(source_label(key.first)).append(" |"));
    }
    return rows;
}

// The doc's rows of the known-differences category table (first three cells).
std::set<std::string> documented_difference_rows() {
    const std::string_view text = group_definitions_text();
    const auto start = text.find("## Known differences");
    const auto end = text.find("\n## ", start + 1);
    std::set<std::string> rows;
    for (const auto line : split(text.substr(start, end - start), '\n')) {
        const auto cells = split(line, '|');
        if (cells.size() == 6 && !integers(trim(cells[2])).empty() &&
            trim(cells[2]).find_first_not_of("0123456789") == std::string_view::npos) {
            std::string row = "| ";
            row.append(trim(cells[1])).append(" | ").append(trim(cells[2])).append(" | ");
            rows.insert(row.append(trim(cells[3])).append(" |"));
        }
    }
    return rows;
}

TEST(Provenance, DocumentedCountsAreTheDataFilesCounts) {
    EXPECT_EQ(documented_difference_rows(), known_difference_rows());
    const auto chebi = data_lines(chebi_text()).size();
    const auto candidates = oracle().size();
    const std::string_view text = group_definitions_text();
    EXPECT_NE(text.find(with_thousands(chebi) + " molecules"), std::string_view::npos)
        << "doc should say " << with_thousands(chebi) << " molecules";
    EXPECT_NE(text.find(with_thousands(candidates) + " PubChem candidates"), std::string_view::npos)
        << "doc should say " << with_thousands(candidates) << " PubChem candidates";
}

// What disagrees between ARISE and RDKit on one valence probe (empty if nothing).
std::string valence_problem(const std::vector<std::string_view>& fields) {
    const bool rdkit_accepts = fields[1] == "ok";
    if (arise_accepts(fields[0]) != rdkit_accepts) {
        return rdkit_accepts ? "accept" : "reject";
    }
    if (!rdkit_accepts) {
        return {};
    }
    const auto molecule = arise::chem::parse_smiles(fields[0]);
    const auto expected = integers(fields[2]);
    return molecule.value().hydrogen_count(0) == expected.front() ? std::string{} : "hydrogens";
}

TEST(Provenance, ValencesMatchRdkitOnRuleGeneratedProbes) {
    // ARISE's normal valences (implicit hydrogens; radical and unusual-valence rejection)
    // against RDKit on every element with 0-6 neighbours, generated by rule in the builder.
    KnownDifferences known{"valence"};
    std::vector<std::string> problems;
    std::size_t checked = 0;
    for (const auto line : data_lines(valence_probes_text())) {
        const auto fields = split(line, '\t');
        ASSERT_EQ(fields.size(), 3U) << line;
        ++checked;
        const auto problem = valence_problem(fields);
        if (!problem.empty() && !known.contains(std::string{fields[0]}, problem)) {
            problems.push_back(std::string{fields[0]} + ": RDKit " + std::string{fields[1]} + " " +
                               std::string{fields[2]} + ", differs in " + problem);
        }
    }
    const auto stale = known.stale();
    problems.insert(problems.end(), stale.begin(), stale.end());
    EXPECT_GT(checked, 300U);
    expect_no_disagreements(problems, checked);
}

TEST(Provenance, DegreesOfUnsaturationMatchRdkitRingsAndPiBonds) {
    KnownDifferences known{"dou"};
    std::vector<std::string> problems;
    std::size_t checked = 0;
    for (const auto line : data_lines(dou_checks_text())) {
        const auto fields = split(line, '\t');
        ASSERT_EQ(fields.size(), 2U) << line;
        const auto formula = arise::core::Formula::parse(fields[0]);
        ASSERT_TRUE(formula.ok()) << line;
        ++checked;
        const auto dou = formula.value().degrees_of_unsaturation();
        const bool same = dou.ok() && dou.value() == integers(fields[1]).front();
        if (!same && !known.contains(std::string{fields[0]}, fields[1])) {
            problems.push_back(std::string{line} + ": ARISE " +
                               (dou.ok() ? std::to_string(dou.value()) : dou.error().message));
        }
    }
    const auto stale = known.stale();
    problems.insert(problems.end(), stale.begin(), stale.end());
    EXPECT_GT(checked, 1000U);
    expect_no_disagreements(problems, checked);
}

// What disagrees between ARISE and RDKit on one spec SMILES RDKit accepts: "accept",
// "hydrogens", "symmetry" or the fragments that differ (empty if nothing).
std::string spec_problem(const OracleEntry& entry, const std::vector<std::string_view>& names,
                         KnownDifferences& known) {
    if (!arise_accepts(entry.smiles)) {
        return "accept";
    }
    const auto molecule = arise::chem::parse_smiles(entry.smiles).value();
    std::vector<int> hydrogens;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        hydrogens.push_back(molecule.hydrogen_count(i));
    }
    if (hydrogens != entry.hydrogens) {
        return "hydrogens";
    }
    if (entry.aromaticity_differs) {
        return {};
    }
    std::vector<int> classes;
    for (const auto cls : arise::chem::symmetry_classes(molecule)) {
        classes.push_back(static_cast<int>(cls));
    }
    if (classes != entry.classes) {
        return "symmetry";
    }
    return fragment_differences(molecule, entry, names, known);
}

// The disagreement on one line of spec_rdkit_oracle.tsv that is not a signed-off
// difference; empty if none.
std::string spec_line_problem(std::string_view line, const std::vector<std::string_view>& names,
                              KnownDifferences& known) {
    const auto fields = split(line, '\t');
    if (fields.size() == 2 && fields[1] == "rejected") {
        const bool listed = known.contains(std::string{fields[0]}, "reject");
        return arise_accepts(fields[0]) && !listed
                   ? std::string{fields[0]} + ": RDKit rejects, ARISE accepts"
                   : std::string{};
    }
    const auto entries = oracle(line);
    if (entries.size() != 1) {
        return "malformed spec line: " + std::string{line};
    }
    const auto& entry = entries.front();
    const auto problem = spec_problem(entry, names, known);
    if (problem.starts_with(' ')) {  // fragment differences, already filtered
        return std::string{entry.smiles} + " (ARISE/RDKit):" + problem;
    }
    if (problem.empty() || known.contains(std::string{entry.smiles}, problem)) {
        return {};
    }
    return std::string{entry.smiles} + ": differs in " + problem;
}

// The SMILES literals passed to the spec helpers in the hand-written chem tests (the same
// helpers as SPEC_HELPERS in the builder).
std::set<std::string> spec_literals() {
    static constexpr std::array kHelpers{
        "parse_smiles(\"",       "parse_ok(\"",    "parse(\"",
        "expect_error(\"",       "groups(\"",      "environments(\"",
        "carbon_class_sizes(\"", "topological(\"",
    };
    const auto sources = arise::test::join_chunks(arise::test::kSpecSourcesChunks);
    std::set<std::string> literals;
    for (const std::string_view helper : kHelpers) {
        for (auto at = sources.find(helper); at != std::string::npos;
             at = sources.find(helper, at + 1)) {
            const bool word_start =
                at == 0 || (std::isalnum(static_cast<unsigned char>(sources[at - 1])) == 0 &&
                            sources[at - 1] != '_');
            const auto begin = at + helper.size();
            if (word_start) {
                std::string literal;  // C++ escapes: "\\" is one backslash in the SMILES
                for (auto i = begin; i < sources.size() && sources[i] != '"'; ++i) {
                    i += sources[i] == '\\' ? 1 : 0;
                    literal += sources[i];
                }
                literals.insert(literal);
            }
        }
    }
    return literals;
}

TEST(Provenance, SpecOracleCoversEverySmilesInTheHandWrittenTests) {
    // A new hand-written case must be re-checked by RDKit: rerun the builder.
    std::set<std::string> covered;
    for (const auto line : data_lines(spec_oracle_text())) {
        covered.insert(std::string{split(line, '\t').front()});
    }
    for (const auto& literal : spec_literals()) {
        EXPECT_TRUE(covered.contains(literal))
            << literal << " is not in spec_rdkit_oracle.tsv; run scripts/build_reference_sets.py";
    }
    EXPECT_GT(spec_literals().size(), 150U);
}

TEST(Provenance, HandWrittenSpecMoleculesMatchRdkit) {
    // Every SMILES literal in the hand-written chem tests (found by pattern by the builder):
    // the tests assert ARISE's output, so ARISE agreeing with RDKit checks the hand-typed
    // expectations for hydrogens, symmetry and the mapped fragments.
    const auto names = header_names(spec_oracle_text(), "# fragments: ");
    ASSERT_FALSE(names.empty()) << "spec oracle has no '# fragments:' line";
    KnownDifferences known{"spec"};
    std::vector<std::string> problems;
    std::size_t checked = 0;
    for (const auto line : data_lines(spec_oracle_text())) {
        ++checked;
        auto problem = spec_line_problem(line, names, known);
        if (!problem.empty()) {
            problems.push_back(std::move(problem));
        }
    }
    const auto stale = known.stale();
    problems.insert(problems.end(), stale.begin(), stale.end());
    EXPECT_GT(checked, 150U);
    expect_no_disagreements(problems, checked);
}

}  // namespace

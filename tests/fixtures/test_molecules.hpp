#pragma once
// The molecules and SMILES-reader strings the hand-written tests use, by name (the
// chemist's rule, 2026-11-21: no test contains a SMILES literal). The data come from
// scripts/build_reference_sets.py:
//   * test_compounds.tsv: PubChem compounds by CID (name, CID, formula, PubChem SMILES,
//     the SMILES tests read, RDKit hybridizations);
//   * rule_compounds.tsv: molecules made by a named rule from a PubChem compound;
//   * smiles_syntax_cases.tsv: notation variants and deliberately broken strings, each
//     checked against RDKit.
// A name that is not listed fails the calling test.

#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <string_view>

namespace arise::test {

namespace detail {

// name -> the text of `column` (0-based) of each data line of `text`.
inline std::map<std::string, std::string, std::less<>> table(std::string_view text,
                                                             std::size_t column) {
    std::map<std::string, std::string, std::less<>> rows;
    while (!text.empty()) {
        const auto end = text.find('\n');
        const auto line = text.substr(0, end);
        text = end == std::string_view::npos ? std::string_view{} : text.substr(end + 1);
        if (line.empty() || line.front() == '#') {
            continue;
        }
        std::size_t start = 0;
        std::string_view field;
        std::string_view name;
        for (std::size_t i = 0; i <= column; ++i) {
            const auto tab = line.find('\t', start);
            field = line.substr(start, tab == std::string_view::npos ? tab : tab - start);
            if (i == 0) {
                name = field;
            }
            start = tab == std::string_view::npos ? line.size() : tab + 1;
        }
        rows.emplace(std::string{name}, std::string{field});
    }
    return rows;
}

inline const std::map<std::string, std::string, std::less<>>& compounds() {
    static const auto rows = [] {
        auto all = table(join_chunks(kTestCompoundsChunks), 4);
        for (auto& [name, smiles] : table(join_chunks(kRuleCompoundsChunks), 4)) {
            all.emplace(name, smiles);
        }
        return all;
    }();
    return rows;
}

inline const std::map<std::string, std::string, std::less<>>& syntax_cases() {
    static const auto rows = table(join_chunks(kSyntaxCasesChunks), 1);
    return rows;
}

}  // namespace detail

/// The SMILES of the test compound `name` (test_compounds.tsv or rule_compounds.tsv).
/// Fails the calling test, and returns an empty string, when `name` is not listed.
inline std::string compound_smiles(std::string_view name) {
    const auto& rows = detail::compounds();
    const auto found = rows.find(name);
    if (found == rows.end()) {
        ADD_FAILURE() << "no test compound named '" << name
                      << "' (add it to TEST_COMPOUNDS or RULE_COMPOUNDS in "
                         "scripts/build_reference_sets.py)";
        return {};
    }
    return found->second;
}

/// The PubChem CID of the test compound `name` (test_compounds.tsv). Fails the calling
/// test, and returns 0, when `name` is not listed.
inline long long compound_cid(std::string_view name) {
    static const auto rows = detail::table(join_chunks(kTestCompoundsChunks), 1);
    const auto found = rows.find(name);
    if (found == rows.end()) {
        ADD_FAILURE() << "no PubChem test compound named '" << name << "'";
        return 0;
    }
    return std::stoll(found->second);
}

/// The string of the SMILES-reader case `name` (smiles_syntax_cases.tsv). Fails the
/// calling test, and returns an empty string, when `name` is not listed.
inline std::string syntax_case(std::string_view name) {
    const auto& rows = detail::syntax_cases();
    const auto found = rows.find(name);
    if (found == rows.end()) {
        ADD_FAILURE() << "no SMILES syntax case named '" << name
                      << "' (add it to SYNTAX_CASES in scripts/build_reference_sets.py)";
        return {};
    }
    return found->second;
}

}  // namespace arise::test

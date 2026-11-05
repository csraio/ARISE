#pragma once
// Reader for the BMRB reference cases (tests/data/reference/bmrb_cases.tsv and
// bmrb_peaks.tsv, built by scripts/build_bmrb_reference.py; see BMRB_NOTICE.md). Shared by
// the tests that check ARISE against them. A malformed line fails the calling test.

#include "arise/io/number.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace arise::test {

// One peak of a BMRB peak list. Atoms are indices into the case's SMILES atoms: in the 1H
// dimension the atom each assigned H is bonded to, in the 13C dimension the carbon.
struct BmrbPeak {
    std::string list;  // "1H", "13C", "DEPT135", "DEPT90" or "HSQC"
    std::optional<double> h_shift;
    std::optional<double> c_shift;
    std::string phase;  // "positive", "negative" or "-"
    std::vector<std::size_t> h_atoms;
    std::vector<std::size_t> c_atoms;
};

// One sample of a BMRB entry, with its structure and every peak of its lists.
struct BmrbCase {
    std::string entry;
    std::string sample;
    std::string solvent;  // as nmr1d::to_string(Solvent)
    std::string compound;
    std::string smiles;
    std::vector<std::string> atom_names;  // BMRB's, per SMILES atom
    std::vector<int> hydrogens;           // bonded H per SMILES atom, from the entry's bonds
    std::vector<BmrbPeak> peaks;          // in file order: by list, then by shift

    [[nodiscard]] std::string id() const { return entry + ":" + sample; }

    // The peaks of one list, in file order.
    [[nodiscard]] std::vector<BmrbPeak> list(std::string_view name) const {
        std::vector<BmrbPeak> found;
        for (const auto& peak : peaks) {
            if (peak.list == name) {
                found.push_back(peak);
            }
        }
        return found;
    }
};

namespace bmrb_detail {

inline std::vector<std::string_view> split(std::string_view text, char separator) {
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

inline std::vector<std::string_view> data_lines(std::string_view text) {
    std::vector<std::string_view> lines;
    for (const auto line : split(text, '\n')) {
        if (!line.empty() && line.front() != '#') {
            lines.push_back(line);
        }
    }
    return lines;
}

inline std::optional<double> shift(std::string_view text, std::string_view line) {
    if (text == "-") {
        return std::nullopt;
    }
    const auto value = io::parse_decimal(text);
    if (!value) {
        ADD_FAILURE() << "bad shift '" << text << "' in: " << line;
        return std::nullopt;
    }
    return value.value();
}

inline std::vector<std::size_t> atoms(std::string_view text, std::size_t atom_count,
                                      std::string_view line) {
    std::vector<std::size_t> found;
    if (text == "-") {
        return found;
    }
    for (const auto item : split(text, ',')) {
        std::size_t value = 0;
        bool digits = !item.empty();
        for (const char c : item) {
            digits = digits && c >= '0' && c <= '9';
            value = (value * 10) + static_cast<std::size_t>(c - '0');
        }
        if (!digits || value >= atom_count) {
            ADD_FAILURE() << "bad atom '" << item << "' in: " << line;
            return {};
        }
        found.push_back(value);
    }
    return found;
}

inline std::vector<BmrbCase> read_cases() {
    std::vector<BmrbCase> cases;
    std::map<std::string, std::size_t, std::less<>> index;
    const std::string case_text = join_chunks(kBmrbCasesChunks);
    for (const auto line : data_lines(case_text)) {
        const auto fields = split(line, '\t');
        if (fields.size() != 7) {
            ADD_FAILURE() << "malformed BMRB case: " << line;
            continue;
        }
        BmrbCase found{
            .entry = std::string{fields[0]},
            .sample = std::string{fields[1]},
            .solvent = std::string{fields[2]},
            .compound = std::string{fields[3]},
            .smiles = std::string{fields[4]},
            .atom_names = {},
            .hydrogens = {},
            .peaks = {},
        };
        for (const auto name : split(fields[5], ' ')) {
            found.atom_names.emplace_back(name);
        }
        for (const auto count : split(fields[6], ' ')) {
            found.hydrogens.push_back(
                count.size() == 1 && count[0] >= '0' && count[0] <= '9' ? count[0] - '0' : -1);
        }
        if (found.hydrogens.size() != found.atom_names.size()) {
            ADD_FAILURE() << "atom names and H counts differ in length: " << line;
            continue;
        }
        index.emplace(found.id(), cases.size());
        cases.push_back(std::move(found));
    }
    const std::string peak_text = join_chunks(kBmrbPeaksChunks);
    for (const auto line : data_lines(peak_text)) {
        const auto fields = split(line, '\t');
        const auto found = fields.size() == 10
                               ? index.find(std::string{fields[0]} + ":" + std::string{fields[1]})
                               : index.end();
        if (found == index.end()) {
            ADD_FAILURE() << "BMRB peak without its case: " << line;
            continue;
        }
        auto& owner = cases[found->second];
        const auto atom_count = owner.atom_names.size();
        owner.peaks.push_back({
            .list = std::string{fields[2]},
            .h_shift = shift(fields[4], line),
            .c_shift = shift(fields[5], line),
            .phase = std::string{fields[6]},
            .h_atoms = atoms(fields[8], atom_count, line),
            .c_atoms = atoms(fields[9], atom_count, line),
        });
    }
    return cases;
}

}  // namespace bmrb_detail

// Every BMRB case, read once per test binary.
inline const std::vector<BmrbCase>& bmrb_cases() {
    static const std::vector<BmrbCase> cases = bmrb_detail::read_cases();
    return cases;
}

}  // namespace arise::test

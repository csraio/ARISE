#pragma once

#include "arise/core/result.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace arise::io {

/// One candidate structure, as written in the input. The SMILES is not yet checked as
/// chemistry; that is `chem`'s job, which can cite `line` in its own errors.
struct Candidate {
    std::string smiles;
    /// The rest of the line after the SMILES, trimmed; empty when there is none.
    std::string name;
    /// 1-based physical line in the input.
    std::size_t line = 0;

    friend bool operator==(const Candidate&, const Candidate&) = default;
};

/// Parses a candidate list in SMILES-file (.smi) form, as written by Sherlock,
/// NMRfilter or RDKit:
/// - Shared line rules: optional byte-order mark, LF or CRLF, blank lines and lines
///   starting with `#` are skipped (a SMILES never starts with `#`, so `C#N` is data).
/// - Each line holds a SMILES, then optionally spaces or tabs and a name.
/// - SMILES may only use ASCII letters, digits and ( ) [ ] = # $ : / \ . + - @ % *.
///   Anything else is a ParseError naming the character (a comma suggests a CSV).
/// - A first line whose SMILES is the word `SMILES` is rejected as a header.
/// - Candidates are returned in file order, duplicates included. At least one.
[[nodiscard]] core::Result<std::vector<Candidate>> parse_candidate_list(std::string_view text);

}  // namespace arise::io

#pragma once

#include "arise/core/result.hpp"

#include <functional>
#include <map>
#include <string>
#include <string_view>

namespace arise::core {

/// A molecular formula: how many atoms of each element.
class Formula {
public:
    /// Parses what a chemist types: element symbols with optional counts, in any order,
    /// repeated elements summed ("CH3COOH" is C2H4O2). Surrounding whitespace is
    /// ignored. Only elements found in organic NMR samples are accepted: C H B N O F Si
    /// P S Cl Se Br I. Errors (ParseError): empty input, a lowercase start, unknown
    /// symbols (with a hint when "Co" probably meant "CO"), zero or leading-digit counts,
    /// brackets, charges and inner whitespace. Counts above 999 are OutOfRange.
    [[nodiscard]] static Result<Formula> parse(std::string_view text);

    /// Number of atoms of `element` ("C", "Cl"); 0 when absent.
    [[nodiscard]] int count(std::string_view element) const noexcept;

    /// Hill notation: C, then H, then the rest alphabetically; without carbon, all
    /// elements alphabetically. Counts of 1 are omitted ("CH4O").
    [[nodiscard]] std::string hill() const;

    /// Rings plus pi bonds, from standard valences (C, Si 4; N, P, B 3; O, S, Se 2;
    /// H and halogens 1): DoU = 1 + sum(n * (valence - 2)) / 2. P and S count with their
    /// lowest valence, so P=O and S=O add to the total, as is conventional.
    /// InvalidArgument when the result is negative or not a whole number (an odd
    /// electron count, e.g. CH3).
    [[nodiscard]] Result<int> degrees_of_unsaturation() const;

    friend bool operator==(const Formula&, const Formula&) = default;

private:
    std::map<std::string, int, std::less<>> counts_;
};

}  // namespace arise::core

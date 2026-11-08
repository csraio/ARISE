// filter::check_confirmed_exchangeable: a confirmed count of exchangeable hydrogens (H on
// O, N, S, Se) eliminates every candidate outside it, with the reason. Hand cases, then an
// external check over every cleaned PubChem candidate of the RDKit oracle
// (tests/data/reference/pubchem_rdkit_oracle.tsv): with the candidate's count taken from
// RDKit's per-atom hydrogens, confirming that count keeps it and any other count eliminates
// it.

#include "arise/filter/exchangeable.hpp"

#include "arise/chem/molecule.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::chem::Element;
using arise::filter::check_confirmed_exchangeable;
using arise::filter::ConfirmedExchangeable;
using arise::filter::Tier;
using arise::filter::Verdict;

arise::chem::Molecule molecule(std::string_view smiles) {
    auto parsed = arise::chem::parse_smiles(smiles);
    EXPECT_TRUE(parsed) << smiles;
    return parsed ? std::move(parsed).value() : arise::chem::Molecule{};
}

Verdict check(std::string_view smiles, ConfirmedExchangeable confirmed) {
    const auto result = check_confirmed_exchangeable(molecule(smiles), confirmed);
    EXPECT_TRUE(result) << smiles << ": " << (result ? "" : result.error().message);
    return result ? result.value() : Verdict{.tier = Tier::BackBurner, .reasons = {}};
}

Verdict active() {
    return {.tier = Tier::Active, .reasons = {}};
}

TEST(ConfirmedExchangeable, AnExactCountKeepsMatchingCandidatesAndEliminatesTheRest) {
    // C2H6O: ethanol has one OH, dimethyl ether none.
    const ConfirmedExchangeable one{.min = 1, .max = 1};
    EXPECT_EQ(check("CCO", one), active());
    const Verdict ether{
        .tier = Tier::Eliminated,
        .reasons = {"confirmed exchangeable H (on O, N, S, Se): 1; the candidate has 0"},
    };
    EXPECT_EQ(check("COC", one), ether);
}

TEST(ConfirmedExchangeable, ARangeKeepsEveryCountInsideIt) {
    const ConfirmedExchangeable one_to_three{.min = 1, .max = 3};
    EXPECT_EQ(check("CCO", one_to_three).tier, Tier::Active);       // 1
    EXPECT_EQ(check("NCC(=O)O", one_to_three).tier, Tier::Active);  // 3
    const auto none = check("COC", one_to_three);
    EXPECT_EQ(none.tier, Tier::Eliminated);
    EXPECT_EQ(none.reasons,
              std::vector<std::string>{
                  "confirmed exchangeable H (on O, N, S, Se): 1 to 3; the candidate has 0"});
    EXPECT_EQ(check("NC(CO)C(=O)O", one_to_three).tier, Tier::Eliminated);  // serine: 4
}

TEST(ConfirmedExchangeable, NoUpperBoundKeepsEveryCountFromTheMinimumUp) {
    const ConfirmedExchangeable at_least_one{.min = 1, .max = std::nullopt};
    EXPECT_EQ(check("NC(CO)C(=O)O", at_least_one).tier, Tier::Active);
    const auto none = check("COC", at_least_one);
    EXPECT_EQ(none.reasons,
              std::vector<std::string>{
                  "confirmed exchangeable H (on O, N, S, Se): at least 1; the candidate has 0"});
}

TEST(ConfirmedExchangeable, ZeroConfirmedEliminatesEveryOhNhShOrSeH) {
    const ConfirmedExchangeable zero{.min = 0, .max = 0};
    EXPECT_EQ(check("COC", zero), active());
    EXPECT_EQ(check("CCO", zero).tier, Tier::Eliminated);
    EXPECT_EQ(check("CS", zero).tier, Tier::Eliminated);
    EXPECT_EQ(check("C[SeH]", zero).tier, Tier::Eliminated);
    EXPECT_EQ(check("CP", zero), active());  // P-H does not exchange
}

TEST(ConfirmedExchangeable, ZeroToNothingKeepsEveryCandidate) {
    const ConfirmedExchangeable anything{.min = 0, .max = std::nullopt};
    for (const auto* smiles : {"COC", "CCO", "NC(CO)C(=O)O"}) {
        EXPECT_EQ(check(smiles, anything), active()) << smiles;
    }
}

TEST(ConfirmedExchangeable, ANegativeMinimumOrAMaximumBelowItIsAnError) {
    const auto ethanol = molecule("CCO");
    for (const auto confirmed : {
             ConfirmedExchangeable{.min = -1, .max = 1},
             ConfirmedExchangeable{.min = 2, .max = 1},
             ConfirmedExchangeable{.min = -1, .max = std::nullopt},
         }) {
        const auto result = check_confirmed_exchangeable(ethanol, confirmed);
        ASSERT_FALSE(result);
        EXPECT_EQ(result.error().code, arise::core::ErrorCode::InvalidArgument);
    }
}

// ---- external: RDKit's hydrogens on every PubChem oracle candidate ---------------------

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

// RDKit's H on O, N, S and Se (per-atom hydrogens from the oracle; the element of each
// atom from the SMILES), or -1 when the line does not fit the molecule.
int rdkit_exchangeable(const arise::chem::Molecule& parsed, std::string_view hydrogens) {
    const auto counts = split(hydrogens, ' ');
    if (counts.size() != parsed.atoms().size()) {
        return -1;
    }
    int total = 0;
    for (std::size_t atom = 0; atom < counts.size(); ++atom) {
        const auto element = parsed.atoms()[atom].element;
        if (element == Element::O || element == Element::N || element == Element::S ||
            element == Element::Se) {
            total += std::stoi(std::string{counts[atom]});
        }
    }
    return total;
}

// "" when confirming the candidate's RDKit count keeps it and confirming more eliminates it;
// otherwise what went wrong.
std::string oracle_difference(std::string_view line) {
    const auto fields = split(line, '\t');
    if (fields.size() < 2) {
        return "malformed oracle line";
    }
    const auto parsed = molecule(fields[0]);
    const int count = rdkit_exchangeable(parsed, fields[1]);
    if (count < 0) {
        return "H counts do not fit the SMILES";
    }
    const auto keep = check_confirmed_exchangeable(parsed, {.min = count, .max = count});
    const auto above = check_confirmed_exchangeable(parsed, {.min = count + 1, .max = {}});
    if (!keep || keep.value().tier != Tier::Active) {
        return "confirming its RDKit count " + std::to_string(count) + " does not keep it";
    }
    if (!above || above.value().tier != Tier::Eliminated) {
        return "confirming more than " + std::to_string(count) + " does not eliminate it";
    }
    return "";
}

TEST(ConfirmedExchangeable, OracleCandidatesSurviveExactlyTheirRdkitCount) {
    const std::string text = arise::test::join_chunks(arise::test::kRdkitOracleChunks);
    std::size_t checked = 0;
    for (const auto line : split(text, '\n')) {
        if (!line.empty() && line.front() != '#') {
            EXPECT_EQ(oracle_difference(line), "") << line.substr(0, line.find('\t'));
            ++checked;
        }
    }
    EXPECT_GT(checked, 4000U);
}

}  // namespace

// bench/snapshot (plan 3d): the record of results the speed work must not change.

#include "snapshot.hpp"

#include "pubchem_fixture.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::bench::fingerprint;
using arise::bench::split_searches;

arise::bench::BenchmarkCase benchmark_case() {
    const auto parsed = arise::bench::parse_benchmark_case(
        arise::test::join_chunks(arise::test::kBenchmarkSpectrumChunks));
    EXPECT_TRUE(parsed);
    return parsed ? parsed.value() : arise::bench::BenchmarkCase{};
}

TEST(SnapshotFingerprint, MatchesPublishedFnv1aTestVectors) {
    // 64-bit FNV-1a test vectors of the IETF draft "The FNV Non-Cryptographic Hash
    // Algorithm" (draft-eastlake-fnv), checked against an independent Python computation.
    EXPECT_EQ(fingerprint(""), "cbf29ce484222325");
    EXPECT_EQ(fingerprint("a"), "af63dc4c8601ec8c");
    EXPECT_EQ(fingerprint("foobar"), "85944171f73967e8");
}

TEST(SnapshotSearches, SplitsByFormulaHeaders) {
    const std::string all = "ignored\n## A1\nrow 1\nrow 2\n## B2\nrow 3\n";
    const auto searches = split_searches(all);
    ASSERT_EQ(searches.size(), 2U);
    EXPECT_EQ(searches[0].first, "A1");
    EXPECT_EQ(searches[0].second, "row 1\nrow 2\n");
    EXPECT_EQ(searches[1].first, "B2");
    EXPECT_EQ(searches[1].second, "row 3\n");
    EXPECT_TRUE(split_searches("").empty());
    EXPECT_TRUE(split_searches("no headers\n").empty());
}

TEST(SnapshotSearches, FindsEveryRecordedSearch) {
    const auto all = arise::test::join_chunks(arise::test::kPubchemSearchesChunks);
    const auto searches = split_searches(all);
    EXPECT_EQ(searches.size(), 8U);
    const auto has = [&](std::string_view formula) {
        return std::any_of(searches.begin(), searches.end(),
                           [&](const auto& s) { return s.first == formula; });
    };
    EXPECT_TRUE(has("C4H8O2"));
    EXPECT_TRUE(has(benchmark_case().formula));
}

// The number of tab-separated fields of each result line of the C4H8O2 search, read as if
// it were the search of `as_case`.
std::vector<long> field_counts(const arise::bench::BenchmarkCase& as_case) {
    const std::vector<arise::bench::Search> searches{{"C4H8O2", arise::test::kPubChemC4H8O2}};
    const auto results = arise::bench::results(searches, as_case);
    EXPECT_TRUE(results && results.value().size() == 1);
    std::vector<long> counts;
    for (const auto& line : results ? results.value()[0].lines : std::vector<std::string>{}) {
        counts.push_back(std::count(line.begin(), line.end(), '\t') + 1);
    }
    return counts;
}

TEST(SnapshotResults, OneLinePerCandidateAndFitsOnlyForTheCasesSearch) {
    const auto imported = arise::bench::clean(arise::test::kPubChemC4H8O2);
    ASSERT_TRUE(imported);
    const auto candidates = arise::bench::parse(imported.value()).molecules.size();
    const auto plain = field_counts(benchmark_case());
    EXPECT_EQ(plain, std::vector<long>(candidates, 3));  // SMILES, fewest, most
    auto as_case = benchmark_case();
    as_case.formula = "C4H8O2";  // now the case's search: one fit field per scenario
    EXPECT_EQ(field_counts(as_case), std::vector<long>(candidates, 6));
}

TEST(SnapshotResults, RejectsAnUnreadableSearch) {
    const std::vector<arise::bench::Search> searches{{"X", "not a PubChem CSV\n"}};
    EXPECT_FALSE(arise::bench::results(searches, benchmark_case()));
}

TEST(SnapshotBody, OneHeaderPerSearchAndOneFingerprintPerLine) {
    const std::vector<arise::bench::SearchResults> all{
        {.formula = "A1", .lines = {"x", "y"}},
        {.formula = "B2", .lines = {}},
    };
    EXPECT_EQ(arise::bench::snapshot(all),
              "## A1 2\n" + fingerprint("x") + "\n" + fingerprint("y") + "\n## B2 0\n");
}

#ifndef __EMSCRIPTEN__
// The first candidate whose fingerprint differs from the committed snapshot's, with its
// recomputed result line; empty when all match.
std::string first_difference(const std::vector<arise::bench::SearchResults>& all,
                             const std::string& committed) {
    const auto expected = arise::bench::snapshot(all);
    const auto body = committed.substr(committed.find("\n## ") + 1);
    if (body == expected) {
        return {};
    }
    std::size_t line = 0;
    for (const auto& search : all) {
        ++line;  // the "## " line
        for (const auto& result : search.lines) {
            const auto at = body.find(fingerprint(result));
            if (at == std::string::npos) {
                return search.formula + " candidate line " + std::to_string(line) + ": " + result;
            }
            ++line;
        }
    }
    return "the searches or their order changed";
}
#endif

TEST(BenchmarkReference, ResultsMatchTheEquivalenceSnapshot) {
    // tests/golden/equivalence_snapshot.tsv, written by bench/arise_snapshot before the speed
    // work: nmr_equivalence and the carbon fit must give exactly the same results.
#ifdef __EMSCRIPTEN__
    GTEST_SKIP() << "every recorded search; far slower under wasm checks";
#else
    const auto all = arise::test::join_chunks(arise::test::kPubchemSearchesChunks);
    const auto results = arise::bench::results(split_searches(all), benchmark_case());
    ASSERT_TRUE(results);
    const auto committed = arise::test::join_chunks(arise::test::kEquivalenceSnapshotChunks);
    EXPECT_TRUE(committed.starts_with(arise::bench::snapshot_header()))
        << "the snapshot's header changed; regenerate on purpose with bench/arise_snapshot";
    EXPECT_EQ(first_difference(results.value(), committed), "");
#endif
}

}  // namespace

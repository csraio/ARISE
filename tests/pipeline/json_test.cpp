// pipeline JSON: the exact text the page reads. Values are illustrative, not measured
// data; only the layout matters here. The page's own tests (Node) check that JavaScript's
// JSON.parse reads it.

#include "arise/pipeline/json.hpp"

#include "arise/chem/perception.hpp"
#include "arise/core/error.hpp"
#include "arise/core/units.hpp"
#include "arise/filter/groups.hpp"
#include "arise/io/pubchem.hpp"
#include "arise/nmr1d/dept.hpp"
#include "arise/pipeline/pipeline.hpp"
#include "test_molecules.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>

namespace {

using arise::core::Ppm;
using arise::nmr1d::CarbonType;
using arise::nmr1d::CarbonTypeSet;
using arise::pipeline::json_string;
using arise::pipeline::to_json;
using arise::pipeline::url_json;

TEST(PipelineJson, StringsEscapeQuotesBackslashesAndControlCharacters) {
    EXPECT_EQ(json_string("plain"), "\"plain\"");
    EXPECT_EQ(json_string(""), "\"\"");
    EXPECT_EQ(json_string("x\"y\\z"), "\"x\\\"y\\\\z\"");
    EXPECT_EQ(json_string("line\nnext\r\ttab"), "\"line\\nnext\\r\\ttab\"");
    EXPECT_EQ(json_string(std::string{"\x01\x1f", 2}), "\"\\u0001\\u001f\"");
    EXPECT_EQ(json_string(std::string{"nul\0x", 5}), "\"nul\\u0000x\"");
    // UTF-8 passes through unchanged ("delta 170" with a Greek delta).
    EXPECT_EQ(json_string("\xce\xb4 170"), "\"\xce\xb4 170\"");
    // A slash needs no escaping; a backslash (also a SMILES stereo bond) does.
    EXPECT_EQ(json_string("x/y\\z"), "\"x/y\\\\z\"");
}

TEST(PipelineJson, ErrorNamesItsCode) {
    const arise::core::Error error{
        .code = arise::core::ErrorCode::ParseError,
        .message = "line 2: \"x\"",
    };
    EXPECT_EQ(to_json(error), "{\"error\":{\"code\":\"" +
                                  std::string{arise::core::to_string(error.code)} +
                                  "\",\"message\":\"line 2: \\\"x\\\"\"}}");
}

TEST(PipelineJson, Url) {
    EXPECT_EQ(url_json("https://example.org/a?b=c"), "{\"url\":\"https://example.org/a?b=c\"}");
}

TEST(PipelineJson, PeaksListTheirPossibleTypesInOrder) {
    const std::vector<arise::nmr1d::CarbonAssignment> peaks{
        {.shift = Ppm{29.5}, .possible = CarbonTypeSet{CarbonType::CH2}},
        {
            .shift = Ppm{170.0},
            .possible = CarbonTypeSet{CarbonType::CH3} | CarbonTypeSet{CarbonType::C},
        },
        {.shift = Ppm{-1.25}, .possible = CarbonTypeSet::all()},
    };
    EXPECT_EQ(to_json(peaks),
              "[{\"shift\":29.5,\"types\":[\"CH2\"]},"
              "{\"shift\":170,\"types\":[\"C\",\"CH3\"]},"
              "{\"shift\":-1.25,\"types\":[\"C\",\"CH\",\"CH2\",\"CH3\"]}]");
    EXPECT_EQ(to_json(std::vector<arise::nmr1d::CarbonAssignment>{}), "[]");
}

TEST(PipelineJson, ShiftsUseAtMostSixDecimals) {
    const std::vector<arise::nmr1d::CarbonAssignment> peaks{
        {.shift = Ppm{0.1234567}, .possible = CarbonTypeSet{CarbonType::C}},
        {.shift = Ppm{0.0}, .possible = CarbonTypeSet{CarbonType::C}},
        {.shift = Ppm{-0.0000001}, .possible = CarbonTypeSet{CarbonType::C}},
    };
    // A value that rounds to zero is written "0", never "-0".
    EXPECT_EQ(to_json(peaks),
              "[{\"shift\":0.123457,\"types\":[\"C\"]},"
              "{\"shift\":0,\"types\":[\"C\"]},{\"shift\":0,\"types\":[\"C\"]}]");
}

TEST(PipelineJson, Suggestions) {
    arise::filter::HallmarkWindow window;
    window.group = arise::chem::Group::Ketone;
    window.hallmark = "carbonyl carbon";
    window.types = CarbonTypeSet{CarbonType::C};
    window.low = Ppm{195.5};
    window.high = Ppm{215.25};
    window.n = 812;
    window.molecules = 790;
    window.source = arise::filter::WindowSource::Percentiles;
    auto lab = window;
    lab.hallmark = "our \"ketone\"";
    lab.source = arise::filter::WindowSource::Lab;
    lab.n = 0;
    lab.molecules = 0;
    auto few = window;
    few.source = arise::filter::WindowSource::MeasuredRange;
    const std::vector<arise::filter::GroupSuggestion> suggestions{
        {
            .group = arise::chem::Group::Ketone,
            .evidence =
                {
                    {.peak = 3, .shift = Ppm{208.5}, .window = window},
                    {.peak = 3, .shift = Ppm{208.5}, .window = lab},
                    {.peak = 4, .shift = Ppm{199}, .window = few},
                },
        },
    };
    const std::string ketone{arise::chem::to_string(arise::chem::Group::Ketone)};
    const std::string head =
        R"({"peak":3,"shift":208.5,"hallmark":"carbonyl carbon","types":["C"],"low":195.5,"high":215.25,)";
    EXPECT_EQ(to_json(suggestions),
              R"([{"group":")" + ketone + R"(","evidence":[)" + head +
                  R"("n":812,"molecules":790,"source":"percentiles"},)" +
                  R"({"peak":3,"shift":208.5,"hallmark":"our \"ketone\"","types":["C"],)" +
                  R"("low":195.5,"high":215.25,"n":0,"molecules":0,"source":"lab"},)" +
                  R"({"peak":4,"shift":199,"hallmark":"carbonyl carbon","types":["C"],)" +
                  R"("low":195.5,"high":215.25,"n":812,"molecules":790,)" +
                  R"("source":"measured range"}]}])");
    EXPECT_EQ(to_json(std::vector<arise::filter::GroupSuggestion>{}), "[]");
}

TEST(PipelineJson, Ranking) {
    using arise::test::compound_cid;
    using arise::test::compound_smiles;
    const auto acetone = compound_smiles("Acetone");
    const auto oxide = compound_smiles("Propylene Oxide");
    const auto mixture_smiles = acetone + "." + oxide;
    arise::pipeline::Ranking ranking;
    ranking.formula = "C3H6O";
    ranking.rows = 5;
    ranking.removed = {
        {.cid = 7, .smiles = mixture_smiles, .reason = arise::io::RemovalReason::Mixture},
        {.cid = 9, .smiles = acetone, .reason = arise::io::RemovalReason::DuplicateConnectivity},
    };
    arise::pipeline::RankedCandidate first;
    first.cid = static_cast<std::uint64_t>(compound_cid("Acetone"));
    first.smiles = acetone;
    first.verdict = {.tier = arise::filter::Tier::Active, .reasons = {}};
    first.assumed_overlaps = 1;
    first.environments_fewest = 2;
    first.environments_most = 2;
    first.notes = {"assumes 1 overlapping peak"};
    arise::pipeline::RankedCandidate second;
    second.cid = static_cast<std::uint64_t>(compound_cid("Propylene Oxide"));
    second.smiles = oxide;
    second.verdict = {.tier = arise::filter::Tier::Eliminated, .reasons = {"no \"ketone\""}};
    ranking.candidates = {first, second};

    const std::string mixture{arise::io::to_string(arise::io::RemovalReason::Mixture)};
    const std::string duplicate{
        arise::io::to_string(arise::io::RemovalReason::DuplicateConnectivity)};
    const std::string active{arise::filter::to_string(arise::filter::Tier::Active)};
    const std::string eliminated{arise::filter::to_string(arise::filter::Tier::Eliminated)};
    // None of these SMILES contains a character JSON escapes.
    EXPECT_EQ(to_json(ranking),
              R"({"formula":"C3H6O","rows":5,"removed":[{"cid":7,"smiles":")" + mixture_smiles +
                  R"(","reason":")" + mixture + R"("},{"cid":9,"smiles":")" + acetone +
                  R"(","reason":")" + duplicate + R"("}],"candidates":[{"cid":)" +
                  std::to_string(first.cid) + R"(,"smiles":")" + acetone + R"(","tier":")" +
                  active + R"(","reasons":[],"notes":["assumes 1 overlapping peak"],)" +
                  R"("assumed_overlaps":1,"environments_fewest":2,"environments_most":2},)" +
                  R"({"cid":)" + std::to_string(second.cid) + R"(,"smiles":")" + oxide +
                  R"(","tier":")" + eliminated + R"(","reasons":["no \"ketone\""],"notes":[],)" +
                  R"("assumed_overlaps":0,"environments_fewest":0,"environments_most":0}]})");
}

TEST(PipelineJson, EmptyRanking) {
    arise::pipeline::Ranking ranking;
    ranking.formula = "C2H6";
    EXPECT_EQ(to_json(ranking), R"({"formula":"C2H6","rows":0,"removed":[],"candidates":[]})");
}

}  // namespace

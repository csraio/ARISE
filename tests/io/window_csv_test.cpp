#include "arise/io/window_csv.hpp"

#include "arise/core/error.hpp"
#include "arise/core/units.hpp"

#include <gtest/gtest.h>

#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::core::Ppm;
using arise::io::parse_window_csv;
using arise::io::WindowRow;

void expect_error(std::string_view text, ErrorCode code,
                  const std::vector<std::string>& fragments) {
    const auto result = parse_window_csv(text);
    ASSERT_FALSE(result.ok()) << "accepted:\n" << text;
    EXPECT_EQ(result.error().code, code) << result.error().message;
    for (const auto& fragment : fragments) {
        EXPECT_NE(result.error().message.find(fragment), std::string::npos)
            << "message '" << result.error().message << "' lacks '" << fragment << "'";
    }
}

TEST(WindowCsv, ParsesTheRequiredColumnsInFileOrder) {
    const auto result =
        parse_window_csv("group,low_ppm,high_ppm\nketone,190,220\nnitrile,110,125\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<WindowRow> expected{
        {
            .line = 2,
            .group = "ketone",
            .hallmark = "",
            .dept = "",
            .low = Ppm{190},
            .high = Ppm{220},
        },
        {
            .line = 3,
            .group = "nitrile",
            .hallmark = "",
            .dept = "",
            .low = Ppm{110},
            .high = Ppm{125},
        },
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(WindowCsv, ReadsOptionalColumnsInAnyOrderAndKeepsTextAsWritten) {
    const auto result = parse_window_csv(
        "\xEF\xBB\xBF# our lab\r\nhigh_ppm, dept ,group,hallmark,low_ppm\r\n"
        "220,C,Ketone,C=O carbon,190\r\n125,,nitrile,,110\r\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    ASSERT_EQ(result.value().size(), 2U);
    EXPECT_EQ(result.value()[0].line, 3U);
    EXPECT_EQ(result.value()[0].group, "Ketone");
    EXPECT_EQ(result.value()[0].hallmark, "C=O carbon");
    EXPECT_EQ(result.value()[0].dept, "C");
    EXPECT_EQ(result.value()[0].low, Ppm{190});
    EXPECT_EQ(result.value()[0].high, Ppm{220});
    EXPECT_EQ(result.value()[1].hallmark, "");
    EXPECT_EQ(result.value()[1].dept, "");
}

TEST(WindowCsv, AcceptsAOnePointWindow) {
    const auto result = parse_window_csv("group,low_ppm,high_ppm\nketone,200,200\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_EQ(result.value()[0].low, result.value()[0].high);
}

TEST(WindowCsv, RejectsMalformedFiles) {
    expect_error("", ErrorCode::ParseError, {"missing header"});
    expect_error("group,low_ppm\nketone,190\n", ErrorCode::ParseError,
                 {"missing required column 'high_ppm'"});
    expect_error("group,low_ppm,high_ppm,colour\nketone,190,220,red\n", ErrorCode::ParseError,
                 {"unknown column 'colour'"});
    expect_error("group,low_ppm,high_ppm\n", ErrorCode::ParseError, {"line 1", "no rows"});
    expect_error("group,low_ppm,high_ppm\nketone,190\n", ErrorCode::ParseError,
                 {"line 2", "expected 3 fields"});
    expect_error("group,low_ppm,high_ppm\n,190,220\n", ErrorCode::ParseError,
                 {"line 2", "empty group"});
    expect_error("group,low_ppm,high_ppm\nketone,abc,220\n", ErrorCode::ParseError,
                 {"line 2", "low_ppm"});
    expect_error("group,low_ppm,high_ppm\nketone,190,\n", ErrorCode::ParseError,
                 {"line 2", "empty high_ppm"});
    expect_error("group,low_ppm,high_ppm\n\"ketone\",190,220\n", ErrorCode::ParseError,
                 {"line 2", "quoted"});
}

TEST(WindowCsv, RejectsImplausibleOrReversedWindows) {
    expect_error("group,low_ppm,high_ppm\nketone,190,2200\n", ErrorCode::OutOfRange,
                 {"line 2", "high_ppm", "2200"});
    expect_error("group,low_ppm,high_ppm\nketone,-90,20\n", ErrorCode::OutOfRange,
                 {"line 2", "low_ppm", "-90"});
    expect_error("group,low_ppm,high_ppm\nketone,220,190\n", ErrorCode::OutOfRange,
                 {"line 2", "low_ppm 220 is above high_ppm 190"});
    expect_error("group,low_ppm,high_ppm\nketone,200.5,200\n", ErrorCode::OutOfRange,
                 {"line 2", "low_ppm 200.5 is above high_ppm 200"});
}

}  // namespace

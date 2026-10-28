#include "arise/core/error.hpp"
#include "arise/core/result.hpp"

#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <utility>

namespace {

using arise::core::Error;
using arise::core::ErrorCode;
using arise::core::Result;
using arise::core::to_string;

Result<int> parse_positive(int input) {
    if (input <= 0) {
        return Error{ErrorCode::InvalidArgument, "input must be positive"};
    }
    return input;
}

TEST(Result, HoldsValueOnSuccess) {
    const Result<int> result = parse_positive(42);
    ASSERT_TRUE(result.ok());
    EXPECT_TRUE(static_cast<bool>(result));
    EXPECT_EQ(result.value(), 42);
}

TEST(Result, HoldsErrorOnFailure) {
    const Result<int> result = parse_positive(-1);
    ASSERT_FALSE(result.ok());
    EXPECT_FALSE(static_cast<bool>(result));
    EXPECT_EQ(result.error().code, ErrorCode::InvalidArgument);
    EXPECT_EQ(result.error().message, "input must be positive");
}

TEST(Result, ValueOrReturnsValueOnSuccess) {
    EXPECT_EQ(parse_positive(7).value_or(0), 7);
}

TEST(Result, ValueOrReturnsFallbackOnFailure) {
    EXPECT_EQ(parse_positive(-7).value_or(99), 99);
}

TEST(Result, ValueIsMutableThroughNonConstAccess) {
    Result<std::string> result{std::string{"abc"}};
    ASSERT_TRUE(result.ok());
    result.value() += "def";
    EXPECT_EQ(result.value(), "abcdef");
}

TEST(Result, SupportsMoveOnlyTypes) {
    Result<std::unique_ptr<int>> result{std::make_unique<int>(5)};
    ASSERT_TRUE(result.ok());
    std::unique_ptr<int> owned = std::move(result).value();
    ASSERT_NE(owned, nullptr);
    EXPECT_EQ(*owned, 5);
}

TEST(Result, MoveOnlyValueOrUsesFallbackOnFailure) {
    Result<std::unique_ptr<int>> result{Error{ErrorCode::NotFound, "missing"}};
    std::unique_ptr<int> owned = std::move(result).value_or(std::make_unique<int>(9));
    ASSERT_NE(owned, nullptr);
    EXPECT_EQ(*owned, 9);
}

TEST(Result, IsCopyableWhenValueIsCopyable) {
    const Result<std::string> original{std::string{"copy me"}};
    const Result<std::string> copy = original;  // NOLINT(performance-unnecessary-copy-initialization)
    ASSERT_TRUE(copy.ok());
    EXPECT_EQ(copy.value(), "copy me");
    EXPECT_EQ(original.value(), "copy me");
}

TEST(Result, ErrorCanBeMovedOut) {
    Result<int> result{Error{ErrorCode::ParseError, "line 12: unexpected token"}};
    ASSERT_FALSE(result.ok());
    const Error error = std::move(result).error();
    EXPECT_EQ(error.code, ErrorCode::ParseError);
    EXPECT_EQ(error.message, "line 12: unexpected token");
}

// ---- Result<void> ------------------------------------------------------------

Result<void> validate(bool should_pass) {
    if (!should_pass) {
        return Error{ErrorCode::OutOfRange, "value out of range"};
    }
    return {};
}

TEST(ResultVoid, DefaultIsSuccess) {
    const Result<void> result{};
    EXPECT_TRUE(result.ok());
}

TEST(ResultVoid, ReportsSuccess) {
    EXPECT_TRUE(validate(true).ok());
}

TEST(ResultVoid, ReportsFailureWithError) {
    const Result<void> result = validate(false);
    ASSERT_FALSE(result.ok());
    EXPECT_FALSE(static_cast<bool>(result));
    EXPECT_EQ(result.error().code, ErrorCode::OutOfRange);
    EXPECT_EQ(result.error().message, "value out of range");
}

// ---- ErrorCode names ---------------------------------------------------------

constexpr std::array kAllCodes{
    ErrorCode::InvalidArgument, ErrorCode::ParseError,     ErrorCode::NotFound,
    ErrorCode::OutOfRange,      ErrorCode::ChemistryError, ErrorCode::Internal,
};

TEST(ErrorCode, EveryCodeHasANonEmptyName) {
    for (const ErrorCode code : kAllCodes) {
        EXPECT_FALSE(to_string(code).empty());
    }
}

TEST(ErrorCode, NamesAreUnique) {
    std::set<std::string_view> names;
    for (const ErrorCode code : kAllCodes) {
        names.insert(to_string(code));
    }
    EXPECT_EQ(names.size(), kAllCodes.size());
}

TEST(ErrorCode, KnownNames) {
    EXPECT_EQ(to_string(ErrorCode::InvalidArgument), "InvalidArgument");
    EXPECT_EQ(to_string(ErrorCode::ParseError), "ParseError");
    EXPECT_EQ(to_string(ErrorCode::ChemistryError), "ChemistryError");
}

}  // namespace

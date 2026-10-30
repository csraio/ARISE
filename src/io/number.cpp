#include "arise/io/number.hpp"

#include "arise/core/error.hpp"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string>

namespace arise::io {
namespace {

bool is_digit(char c) noexcept {
    return c >= '0' && c <= '9';
}

bool is_sign(char c) noexcept {
    return c == '+' || c == '-';
}

// Number of consecutive ASCII digits in `text` starting at `pos`.
std::size_t digit_run(std::string_view text, std::size_t pos) noexcept {
    std::size_t end = pos;
    while (end < text.size() && is_digit(text[end])) {
        ++end;
    }
    return end - pos;
}

// [+-]? (digits ('.' digits?)? | '.' digits) ([eE] [+-]? digits)?
bool is_plain_decimal(std::string_view text) noexcept {
    std::size_t pos = 0;
    if (pos < text.size() && is_sign(text[pos])) {
        ++pos;
    }
    const std::size_t integer_digits = digit_run(text, pos);
    pos += integer_digits;
    std::size_t fraction_digits = 0;
    if (pos < text.size() && text[pos] == '.') {
        ++pos;
        fraction_digits = digit_run(text, pos);
        pos += fraction_digits;
    }
    if (integer_digits + fraction_digits == 0) {
        return false;
    }
    if (pos < text.size() && (text[pos] == 'e' || text[pos] == 'E')) {
        ++pos;
        if (pos < text.size() && is_sign(text[pos])) {
            ++pos;
        }
        const std::size_t exponent_digits = digit_run(text, pos);
        if (exponent_digits == 0) {
            return false;
        }
        pos += exponent_digits;
    }
    return pos == text.size();
}

core::Error parse_error(std::string_view text, std::string_view reason) {
    return core::Error{
        .code = core::ErrorCode::ParseError,
        .message = "'" + std::string{text} + "' " + std::string{reason},
    };
}

}  // namespace

core::Result<double> parse_decimal(std::string_view text) {
    if (!is_plain_decimal(text)) {
        return parse_error(text, "is not a plain decimal number");
    }
    // strtod needs NUL termination, and the grammar above leaves it nothing to
    // reinterpret (no hex, inf or nan). It is correctly rounded in glibc, Apple's libc
    // and musl (Emscripten). Under a locale whose decimal point is not '.', it stops
    // early, which the length check turns into an error instead of a wrong value.
    const std::string terminated{text};
    char* end = nullptr;
    const double value = std::strtod(terminated.c_str(), &end);
    const auto consumed = static_cast<std::size_t>(end - terminated.c_str());
    if (consumed != terminated.size()) {
        return parse_error(text, "could not be converted (the C locale must use '.')");
    }
    if (!std::isfinite(value)) {
        return parse_error(text, "is outside the range of a double");
    }
    return value;
}

}  // namespace arise::io

#pragma once

#include <cstdint>
#include <string>
#include <string_view>

namespace arise::core {

/// Category of a recoverable failure. Exhaustive: `to_string` must handle every value.
enum class ErrorCode : std::uint8_t {
    InvalidArgument,
    ParseError,
    NotFound,
    OutOfRange,
    ChemistryError,
    Internal,
};

/// A recoverable failure: a machine-readable code plus a human-readable message.
struct Error {
    ErrorCode code;
    std::string message;
};

/// Stable name of an error code (e.g. "ParseError"), for logs and user-facing output.
[[nodiscard]] std::string_view to_string(ErrorCode code) noexcept;

}  // namespace arise::core

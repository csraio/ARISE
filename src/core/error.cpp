#include "arise/core/error.hpp"

namespace arise::core {

std::string_view to_string(ErrorCode code) noexcept {
    // No default case: -Wswitch -Werror flags any ErrorCode added without a name here.
    switch (code) {
        case ErrorCode::InvalidArgument:
            return "InvalidArgument";
        case ErrorCode::ParseError:
            return "ParseError";
        case ErrorCode::NotFound:
            return "NotFound";
        case ErrorCode::OutOfRange:
            return "OutOfRange";
        case ErrorCode::ChemistryError:
            return "ChemistryError";
        case ErrorCode::Internal:
            return "Internal";
    }
    return "Unknown";  // Unreachable for valid enumerators.
}

}  // namespace arise::core

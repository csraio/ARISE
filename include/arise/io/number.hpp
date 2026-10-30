#pragma once

#include "arise/core/result.hpp"

#include <string_view>

namespace arise::io {

/// Parses a plain decimal number: optional sign, digits with an optional '.', and an
/// optional exponent (`7.26`, `-0.5`, `.5`, `1e-3`). Everything else is a ParseError:
/// surrounding whitespace, hex, `inf`/`nan`, comma decimals, and values that overflow.
///
/// The grammar is checked here, so the result never depends on the C locale; the
/// conversion itself is correctly rounded, matching the compiler's reading of the same
/// literal on every supported platform (golden rule 4).
[[nodiscard]] core::Result<double> parse_decimal(std::string_view text);

}  // namespace arise::io

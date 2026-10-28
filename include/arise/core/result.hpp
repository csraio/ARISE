#pragma once

#include "arise/core/error.hpp"

#include <cassert>
#include <optional>
#include <utility>
#include <variant>

namespace arise::core {

/// Either a value of type T or an Error. ARISE's exception-free error channel.
///
/// Accessing `value()` on an error (or `error()` on a value) is a contract violation,
/// checked by assert in debug builds. Check `ok()` first.
template <typename T>
class [[nodiscard]] Result {
public:
    // Implicit by design so functions can `return value;` or `return Error{...};`.
    Result(T value) : storage_{std::in_place_index<0>, std::move(value)} {}  // NOLINT
    Result(Error error) : storage_{std::in_place_index<1>, std::move(error)} {}  // NOLINT

    [[nodiscard]] bool ok() const noexcept { return storage_.index() == 0; }
    [[nodiscard]] explicit operator bool() const noexcept { return ok(); }

    [[nodiscard]] T& value() & noexcept {
        assert(ok());
        return *std::get_if<0>(&storage_);
    }
    [[nodiscard]] const T& value() const& noexcept {
        assert(ok());
        return *std::get_if<0>(&storage_);
    }
    [[nodiscard]] T&& value() && noexcept {
        assert(ok());
        return std::move(*std::get_if<0>(&storage_));
    }

    [[nodiscard]] const Error& error() const& noexcept {
        assert(!ok());
        return *std::get_if<1>(&storage_);
    }
    [[nodiscard]] Error&& error() && noexcept {
        assert(!ok());
        return std::move(*std::get_if<1>(&storage_));
    }

    template <typename U>
    [[nodiscard]] T value_or(U&& fallback) const& {
        return ok() ? value() : static_cast<T>(std::forward<U>(fallback));
    }
    template <typename U>
    [[nodiscard]] T value_or(U&& fallback) && {
        return ok() ? std::move(*this).value() : static_cast<T>(std::forward<U>(fallback));
    }

private:
    std::variant<T, Error> storage_;
};

/// Success-or-Error for operations that produce no value. Default-constructed means success.
template <>
class [[nodiscard]] Result<void> {
public:
    Result() noexcept = default;
    Result(Error error) : error_{std::move(error)} {}  // NOLINT

    [[nodiscard]] bool ok() const noexcept { return !error_.has_value(); }
    [[nodiscard]] explicit operator bool() const noexcept { return ok(); }

    [[nodiscard]] const Error& error() const& noexcept {
        assert(!ok());
        return *error_;
    }
    [[nodiscard]] Error&& error() && noexcept {
        assert(!ok());
        return std::move(*error_);
    }

private:
    std::optional<Error> error_;
};

}  // namespace arise::core

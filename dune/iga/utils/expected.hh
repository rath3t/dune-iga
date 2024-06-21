// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace Dune::Std {

// The following is a naive implementation of std::expected, parts of the interface are implemented herein
// std::variant is used for the handling of the expected type
// It is currently restricted to default constructable types

// can be replaced with
template <typename E>
requires(std::is_trivial_v<E>)
struct unexpected
{
  unexpected() = delete;
  constexpr explicit unexpected(const E& e)
      : val_(e) {}
  constexpr explicit unexpected(E&& e)
      : val_(std::move(e)) {}

  constexpr const E& error() const& {
    return val_;
  }
  constexpr E& error() & {
    return val_;
  }
  constexpr const E&& error() const&& {
    return std::move(val_);
  }
  constexpr E&& error() && {
    return std::move(val_);
  }

private:
  E val_;
};

// CTAD
template <typename E>
unexpected(E) -> unexpected<E>;

/**
 * \brief Class template Std::expected is a way to store either a expected value or an unexpected value. The unexpected
 * value is most of the type a class enum containing error code, but is not limited to that. This is a mimicry version
 * of C++23 std::expected, it fulfills part of the proposed interface from
 * https://en.cppreference.com/w/cpp/utility/expected. For usage refer to the example therein.
 *
 * @tparam T Type of the expected Element
 * @tparam E Type of the unexpected Element (error code)
 */
template <typename T, typename E>
requires(std::is_default_constructible_v<T> and std::is_default_constructible_v<E>)
struct expected
{
  using value_type      = T;
  using error_type      = E;
  using unexpected_type = unexpected<E>;

  constexpr expected()
      : val_(T{}),
        has_value_(true) {}

  constexpr expected(const expected& other)
      : val_(other.val_),
        has_value_(other.has_value_) {}

  constexpr expected(const unexpected<E>& e)
      : val_(e),
        has_value_(false) {}
  constexpr expected(unexpected<E>&& e)
      : val_(std::move(e)),
        has_value_(false) {}

  constexpr expected(const T& t)
      : val_(t),
        has_value_(true) {}
  constexpr expected(T&& t)
      : val_(std::move(t)),
        has_value_(true) {}

  constexpr expected& operator=(const expected& other) {
    val_       = other.val_;
    has_value_ = other.val_;
    return *this;
  }

  template <class U = T>
  constexpr expected& operator=(U&& v) {
    val_       = v;
    has_value_ = true;
    return *this;
  }

  template <class G>
  constexpr expected& operator=(const unexpected<G>& other) {
    val_       = other;
    has_value_ = false;
    return *this;
  }

  template <class G>
  constexpr expected& operator=(unexpected<G>&& other) {
    val_       = other;
    has_value_ = false;
    return *this;
  }

  constexpr const T* operator->() const {
    assert(has_value_);
    return &std::get<T>(val_);
  }
  constexpr T* operator->() {
    assert(has_value_);
    return &std::get<T>(val_);
  }

  constexpr const T& operator*() const& {
    assert(has_value_);
    return std::get<T>(val_);
  }
  constexpr T& operator*() & {
    assert(has_value_);
    return std::get<T>(val_);
  }

  constexpr explicit operator bool() const {
    return has_value_;
  }
  constexpr bool has_value() const {
    return has_value_;
  }

  constexpr T& value() & {
    assert(has_value_);
    return std::get<T>(val_);
  }
  constexpr const T& value() const& {
    assert(has_value_);
    return std::get<T>(val_);
  }
  constexpr T&& value() && {
    assert(has_value_);
    return std::move(std::get<T>(val_));
  }

  constexpr const T&& value() const&& {
    assert(has_value_);
    return std::move(std::get<T>(val_));
  }

  constexpr const E& error() const& {
    assert(not has_value_);
    return std::get<unexpected<E>>(val_);
  }
  constexpr E& error() & {
    assert(not has_value_);
    return std::get<unexpected<E>>(val_).error();
  }

  constexpr E&& error() && {
    assert(not has_value_);
    return std::move(std::get<unexpected<E>>(val_).error());
  }

  template <class U>
  requires(std::experimental::is_convertible_v<U, T>)
  constexpr T value_or(U&& default_value) const& {
    if (has_value_)
      return std::get<T>(val_);
    return default_value;
  }

private:
  std::variant<T, unexpected<E>> val_;
  bool has_value_;
};

} // namespace Dune::Std

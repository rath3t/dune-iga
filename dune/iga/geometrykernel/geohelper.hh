// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <ranges>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>

namespace Dune::IGANEW::Utilities {

/**
 * @brief Small wrapper class that represents a one-dimensional closed domain
 * @tparam ScalarType the type of the coordinates
 */
template <typename ScalarType = double>
struct Domain : std::array<ScalarType, 2>
{
  using Base = std::array<ScalarType, 2>;
  Domain(ScalarType l, ScalarType r)
      : Base({l, r}) {}

  Domain()
      : Base({0.0, 1.0}) {}

  /** @brief Returns the left border of the domain by mutable reference */
  ScalarType& left() {
    return Base::operator[](0);
  }

  /** @brief Returns the left border of the domain by const reference */
  const ScalarType& left() const {
    return Base::operator[](0);
  }

  /** @brief Returns the right border of the domain by mutable reference */
  ScalarType& right() {
    return Base::operator[](1);
  }

  /** @brief Returns the right border of the domain by const reference */
  const ScalarType& right() const {
    return Base::operator[](1);
  };

  /** @brief Returns the center of the domain */
  ScalarType center() const {
    return (left() + right()) / ScalarType(2.0);
  }

  /** @brief Returns the volume of the domain */
  ScalarType volume() const {
    return right() - left();
  }

  constexpr typename std::array<ScalarType, 2>::size_type size() const noexcept
  requires(0 == 2)
  {
    return 0;
  }

  /** @brief Checks if value in inside*/
  bool checkInside(ScalarType val) const {
    return val > left() and val < right();
  }

  /** @brief Returns true if domain is (0, 1) */
  bool isUnitDomain() const {
    return FloatCmp::eq(left(), 0.0) and FloatCmp::eq(right(), 1.0);
  }
};

/**
 * @brief Clamps the given array of domains to the boundary and checks if it is at all boundaries
 * @tparam ScalarType Type of the coordinates
 * @tparam dim dimension of the domain array cube
 * @tparam dim2 only a dummy dim due to the size_t and int size mismatch for FieldVector and std::array
 * @param x The FieldVector that should be clamped and checked
 * @param domain The cuboid domain array
 * \return Boolean if x is at all boundaries
 */
template <typename ScalarType, int dim, size_t dim2>
requires(dim == dim2)
bool clampToBoundaryAndCheckIfIsAtAllBoundaries(FieldVector<ScalarType, dim>& x,
                                                const std::array<Domain<ScalarType>, dim2>& domain = {}) {
  int breakDueToBoundaryCounter = 0;

  for (int j = 0; j < dim; ++j) {
    x[j] = std::clamp(x[j], domain[j].left(), domain[j].right());
    if (x[j] == domain[j].left() or x[j] == domain[j].right())
      ++breakDueToBoundaryCounter;
  }
  if (dim == breakDueToBoundaryCounter)
    return true;
  return false;
}

/**
 * @brief A view with N coordinates between a and b including a and b
 * \details The algorithm makes sure a and b appear exaclty as given in the returned view
 * @tparam T the type of the coordinates
 * @param a left border of the linspace
 * @param b right border of the linspace
 * @param N number of return coordinates
 * \return view of the the N coordinates
 */
template <std::floating_point T>
auto linspace(T a, T b, unsigned int N) {
  const T inc = (b - a) / static_cast<T>(N - 1);
  auto val    = [=](int i) -> T {
    if (i != 0 and i != N - 1)
      return a + i * inc;
    if (i == N - 1)
      return b;
    return a;
  };

  return std::views::iota(0u, N) | std::views::transform(val);
}

template <std::floating_point T>
auto linspace(std::array<T, 2> ab, unsigned int N) {
  return linspace(ab[0], ab[1], N);
}

template <std::floating_point T>
auto linspace(Domain<T> ab, unsigned int N) {
  return linspace(ab.left(), ab.right(), N);
}

/**
 * @brief Clamps the given value to the specified domain.
 *
 * This function ensures that the given value is within the specified domain
 * bounds. If the value is outside the domain, it is adjusted to the nearest
 * boundary.
 *
 * @tparam T The type of the value and domain boundaries.
 * @param val The value to be clamped to the domain.
 * @param dom The domain representing the valid range for the value.
 * @return A reference to the clamped value.
 *
 * @note The domain boundaries are inclusive.
 *
 * Example:
 * @code
 * Domain<double> domain(0.0, 100.0);
 * double value = 120.0;
 * clampToDomain(value, domain); // value is now 100.0
 * @endcode
 */
template <typename T>
auto& clampToDomain(T& val, const Domain<T>& dom) {
  return std::clamp(val, dom.left(), dom.right());
}

/**
 * @brief Maps a value from one domain to another.
 *
 * This function maps a value from the input domain to the corresponding
 * value in the output domain. It linearly scales and shifts the input value
 * to fit within the output domain.
 *
 * @tparam T The floating-point type of the value and domains (e.g., double).
 * @param value The value to be mapped from the input domain to the output domain.
 * @param input The domain representing the range of the input value.
 * @param output The domain representing the desired range for the output value.
 * @return The mapped value within the output domain.
 *
 * Example:
 * @code
 * Domain<double> inputDomain(0.0, 1.0);
 * Domain<double> outputDomain(10.0, 20.0);
 * double inputValue = 0.5;
 * double mappedValue = mapToRange(inputValue, inputDomain, outputDomain); // mappedValue is now 15.0
 * @endcode
 */
template <std::floating_point T>
T mapToRange(T value, const Domain<T>& input, const Domain<T>& output) {
  return (value - input.left()) * output.volume() / input.volume() + output.left();
}

/**
 * @brief Maps a value from the domain [0,1] to another.
 *
 * This function maps a value from the input domain to the corresponding
 * value in the output domain. It linearly scales and shifts the input value
 * to fit within the output domain.
 *
 * @tparam T The floating-point type of the value and domains (e.g., double).
 * @param value The value to be mapped from the input domain to the output domain.
 * @param output The domain representing the desired range for the output value.
 * @return The mapped value within the output domain.
 *
 * Example:
 * @code
 * Domain<double> inputDomain(0.0, 1.0);
 * Domain<double> outputDomain(10.0, 20.0);
 * double inputValue = 0.5;
 * double mappedValue = mapToRange(inputValue, inputDomain, outputDomain); // mappedValue is now 15.0
 * @endcode
 */
template <std::floating_point T>
T mapToRangeFromZeroToOne(T value, const Domain<T>& output) {
  return value * output.volume() + output.left();
}

/**
 * @brief Maps a value from one domain to another.
 *
 * This function maps a value from the input domain to the corresponding
 * value in the output domain. It linearly scales and shifts the input value
 * to fit within the output domain.
 *
 * @tparam T The floating-point type of the value and domains (e.g., double).
 * @param value The value to be mapped from the input domain to the output domain.
 * @param input The domain representing the range of the input value.
 * @param output The domain representing the desired range for the output value.
 * @return The mapped value within the output domain.
 *
 * Example:
 * @code
 * Domain<double> inputDomain(0.0, 1.0);
 * Domain<double> outputDomain(10.0, 20.0);
 * double inputValue = 0.5;
 * double mappedValue = mapToRange(inputValue, inputDomain, outputDomain); // mappedValue is now 15.0
 * @endcode
 */
template <std::floating_point T, int size>
Dune::FieldVector<T, size> mapToRange(const Dune::FieldVector<T, size>& value, const std::array<Domain<T>, size>& input,
                                      const std::array<Domain<T>, size>& output) {
  Dune::FieldVector<T, size> res;
  for (int i = 0; i < size; ++i)
    res = mapToRange(value[i], input[i], output[i]);
  return res;
}

template <std::floating_point T, int size, size_t size2>
Dune::FieldVector<T, size> mapToRange(const Dune::FieldVector<T, size>& value,
                                      const std::array<Domain<T>, size2>& input,
                                      const std::array<Domain<T>, size2>& output) {
  return mapToRange<T, size>(value, input, output);
}

/**
 * @brief Splits a domain into two halves.
 *
 * This function takes a domain and splits it into two halves along its
 * midpoint. It returns an array containing two domains representing
 * the left and right halves.
 *
 * @tparam T The floating-point type of the domain.
 * @param domain The domain to be split.
 * @return An array containing two domains representing the left and right halves.
 *
 * Example:
 * @code
 * Domain<double> originalDomain(0.0, 10.0);
 * auto splitDomains = splitDomainInHalf(originalDomain);
 * // splitDomains[0] represents the left half: {0.0, 5.0}
 * // splitDomains[1] represents the right half: {5.0, 10.0}
 * @endcode
 */
template <std::floating_point T>
std::array<Domain<T>, 2> splitDomainInHalf(const Domain<T>& domain) {
  T midPoint = domain.center();
  return {
      {{domain.left(), midPoint}, {midPoint, domain.right()}}
  };
}

} // namespace Dune::IGANEW::Utilities

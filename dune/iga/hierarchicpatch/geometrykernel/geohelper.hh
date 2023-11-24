// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <ranges>

#include <dune/common/fvector.hh>

namespace Dune::IGANEW::Utilities {

  /**
   * \brief Small wrapper class that represents a one-dimensional closed domain
   * \tparam ScalarType the type of the coordinates
   */
  template <typename ScalarType = double>
  struct Domain : std::array<ScalarType, 2> {
    using Base = std::array<ScalarType, 2>;
    Domain(ScalarType l, ScalarType r) : Base({l, r}) {}

    Domain() : Base({0.0, 1.0}) {}

    /** \brief Returns the left border of the domain by mutable reference */
    ScalarType& left() { return Base::operator[](0); }

    /** \brief Returns the left border of the domain by const reference */
    const ScalarType& left() const { return Base::operator[](0); }

    /** \brief Returns the right border of the domain by mutable reference */
    ScalarType& right() { return Base::operator[](1); }

    /** \brief Returns the right border of the domain by const reference */
    const ScalarType& right() const { return Base::operator[](1); };

    /** \brief Returns the center of the domain */
    ScalarType center() const { return (left() + right()) / ScalarType(2.0); }

    /** \brief Returns the size of the domain */
    ScalarType size() const { return right() - left(); }
  };

  /**
   * \brief Clamps the given array of domains to the boundary and checks if it is at all boundaries
   * \tparam ScalarType Type of the coordinates
   * \tparam dim dimension of the domain array cube
   * \param x The FieldVector that should be clamped and checked
   * \param domain The cuboid domain array
   * \return Boolean if x is at all boundaries
   */
  template <typename ScalarType, int dim>
  bool clampToBoundaryAndCheckIfIsAtAllBoundaries(FieldVector<ScalarType, dim>& x,
                                                  const std::array<Domain<ScalarType>, dim>& domain = {}) {
    int breakDueToBoundaryCounter = 0;

    for (int j = 0; j < dim; ++j) {
      x[j] = std::clamp(x[j], domain[j].left(), domain[j].right());
      if (x[j] == domain[j].left() or x[j] == domain[j].right()) ++breakDueToBoundaryCounter;
    }
    if (dim == breakDueToBoundaryCounter) return true;
    return false;
  }

#ifndef DOXYGEN
  template <typename ScalarType, int dim, size_t dim2>
  requires(dim
           == dim2) auto clampToBoundaryAndCheckIfIsAtAllBoundaries(FieldVector<ScalarType, dim>& x,
                                                                    const std::array<Domain<ScalarType>, dim2>& domain)
      -> clampToBoundaryAndCheckIfIsAtAllBoundaries<ScalarType, dim, dim2>;
#endif

  /**
   * \brief A view with N coordinates between a and b including a and b
   * \details The algorithm makes sure a and b appear exaclty as given in the returned view
   * \tparam T the type of the coordinatesfi
   * \param a left border of the linspace
   * \param b right border of the linspace
   * \param N number of return coordinates
   * \return view of the the N coordinates
   */
  template <std::floating_point T>
  auto linspace(T a, T b, unsigned int N) {
    const T inc = (b - a) / static_cast<T>(N - 1);
    auto val    = [=](int i) -> T {
      if (i != 0 and i != N - 1) return a + i * inc;
      if (i == N - 1) return b;
      return a;
    };

    return std::views::iota(0u, N) | std::views::transform(val);
  }

  template <typename T>
  auto& clampToDomain(T& val, const Domain<T>& dom) {
    return std::clamp(val, dom.left(), dom.right());
  }

  template <std::floating_point T>
  auto linspace(std::array<T, 2> ab, unsigned int N) {
    return linspace(ab[0], ab[1], N);
  }

  template <std::floating_point T>
  auto linspace(Domain<T> ab, unsigned int N) {
    return linspace(ab.left(), ab.right(), N);
  }

  /// \brief Maps a value of type T that is defined in one domain into another domain
  template <std::floating_point T>
  T mapToRange(T value, const Domain<T>& input, const Domain<T>& output) {
    return (value - input.left()) * output.size() / input.size() + output.left();
  }

  template <std::floating_point T>
  std::array<Domain<T>, 2> splitDomainInHalf(const Domain<T>& domain) {
    T midPoint = domain.center();
    return {{{domain.left(), midPoint}, {midPoint, domain.right()}}};
  }
}  // namespace Dune::IGANEW::Utilities

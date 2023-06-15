// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <ranges>

#include <dune/common/fvector.hh>

namespace Dune::IGA::Utilities {

  template <typename ScalarType = double>
  struct Domain : public std::array<ScalarType, 2> {
    using Base = std::array<ScalarType, 2>;
    Domain(ScalarType l, ScalarType r) : Base({l, r}) {}

    Domain() : Base({0.0, 1.0}) {}

    ScalarType& left() { return Base::operator[](0); };
    const ScalarType& left() const { return Base::operator[](0); };
    ScalarType& right() { return Base::operator[](1); };
    const ScalarType& right() const { return Base::operator[](1); };

    ScalarType center() const { return (left() + right()) / ScalarType(2.0); }
    ScalarType size() const { return right() - left(); }
  };

  template <typename ScalarType, int dim, size_t dim2>
  requires(dim == dim2) auto clampToBoundaryAndCheckIfIsAtAllBoundaries(
      Dune::FieldVector<ScalarType, dim>& x, const std::array<Domain<ScalarType>, dim2>& domain = {}) {
    int breakDueToBoundaryCounter = 0;

    for (int j = 0; j < dim; ++j) {
      x[j] = std::clamp(x[j], domain[j].left(), domain[j].right());
      if (x[j] == domain[j].left() or x[j] == domain[j].right()) ++breakDueToBoundaryCounter;
    }
    if (dim == breakDueToBoundaryCounter)
      return true;
    else
      return false;
  }

  template <std::floating_point T>
  auto linspace(T a, T b, unsigned int N) {
    T inc    = (b - a) / static_cast<T>(N - 1);
    auto val = [=](int i) -> T {
      if (i != 0 and i != N - 1)
        return a + i * inc;
      else if (i == N - 1)
        return b;
      else
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
}  // namespace Dune::IGA::Utilities

//
// Created by Henri on 12.03.2023.
//
#pragma once

namespace Dune::IGA::Utilities {

//  template <std::floating_point T>
//  std::vector<T> linspace(T a, T b, unsigned int N) {
//    T inc    = (b - a) / static_cast<T>(N - 1);
//    auto val = [a, inc](int i) -> T { return a + i * inc; };
//
//    std::vector<T> xs(N);
//    for (auto i : std::views::iota(1u, N - 1))
//      xs[i] = val(i);
//
//    // Make sure the beginning and end are interpolatory
//    xs.front() = a;
//    xs.back()  = b;
//
//    return xs;
//  }

  template <std::floating_point T>
  auto linspace(T a, T b, unsigned int N) {
    T inc    = (b - a) / static_cast<T>(N - 1);
    auto val = [=](int i) -> T {
      if(i!=0 and i!=N-1)
        return a + i * inc;
      else if(i==N-1)
        return b;
      else
      return a;
    };

    return std::views::iota(0u, N) | std::views::transform(val);
  }

  template <std::floating_point T>
  auto linspace(std::array<T, 2> ab, unsigned int N) {
    return linspace(ab[0], ab[1], N);
  }

  /// \brief Maps a value of type T that is defined in one domain into another domain
  template <std::floating_point T>
  T mapToRange(T value, T inputMin, T inputMax, T outputMin, T outputMax) {
    return (value - inputMin) * (outputMax - outputMin) / (inputMax - inputMin) + outputMin;
  }

  template <std::floating_point T>
  std::array<std::array<T, 2>, 2> splitDomains(std::array<T, 2> domain) {
    T midPoint = (domain[0] + domain[1]) / 2;
    return {std::array<T, 2>{domain[0], midPoint}, std::array<T, 2>{midPoint, domain[1]}};
  }
}  // namespace Dune::IGA::Utilities

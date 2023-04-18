//
// Created by Henri on 12.03.2023.
//
#pragma once

namespace Dune::IGA::Utilities {

  struct Parameters {
    int preSample{};
    bool preSampleOnlyCurvedEdges{};
    int preGlobalRefine{};
    int edgeRefinements{};
  };

  Parameters& getParameters() {
    static Parameters instance;
    return instance;
  }

  void setStandardParameters() {
    Parameters parameters{};

    parameters.preSample                = 3;
    parameters.preSampleOnlyCurvedEdges = true;
    parameters.preGlobalRefine          = 0;
    parameters.edgeRefinements          = 0;

    getParameters() = parameters;
  }


  template <std::floating_point T>
  std::vector<T> linspace(T a, T b, unsigned int N) {
    T inc    = (b - a) / static_cast<T>(N - 1);
    auto val = [a, inc](int i) -> T { return a + i * inc; };

    std::vector<T> xs(N);
    for (auto i : std::views::iota(1u, N - 1))
      xs[i] = val(i);

    // Make sure the beginning and end are interpolatory
    xs.front() = a;
    xs.back()  = b;

    return xs;
  }

  template <std::floating_point T>
  std::vector<T> linspace(std::array<T, 2> ab, unsigned int N) {
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


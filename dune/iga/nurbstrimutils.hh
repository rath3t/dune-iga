//
// Created by Henri on 12.03.2023.
//
#pragma once

#ifndef DUNE_IGA_NURBSTRIMUTILS_HH
#  define DUNE_IGA_NURBSTRIMUTILS_HH

#  include <nlohmann/json.hpp>

namespace Dune::IGA::Utilities {
  using json = nlohmann::json;

  struct Parameters {
    std::string inputFile;
    bool plot{};
    int GLOBAL_gridGlobalRefine{};
    std::string gridImplementation;
    int preSample{};
    bool preSampleOnlyCurvedEdges{};
    int preGlobalRefine{};
    int edgeRefinements{};
  };

  Parameters& getParameters() {
    static Parameters instance;
    return instance;
  }

  void from_json(const json& j, Parameters& parameters) {
    j.at("inputFile").get_to(parameters.inputFile);
    j.at("plot").get_to(parameters.plot);
    j.at("GLOBAL_gridGlobalRefine").get_to(parameters.GLOBAL_gridGlobalRefine);
    j.at("gridImplementation").get_to(parameters.gridImplementation);
    j.at("preSample").get_to(parameters.preSample);
    j.at("preSampleOnlyCurvedEdges").get_to(parameters.preSampleOnlyCurvedEdges);
    j.at("preGlobalRefine").get_to(parameters.preGlobalRefine);
    j.at("edgeRefinements").get_to(parameters.edgeRefinements);
  }
  void setStandardParameters() {
    Utilities::Parameters parameters{};

    // Fill in standard Werte (ohne File)
    parameters.preSample                = 1;
    parameters.preSampleOnlyCurvedEdges = true;
    parameters.preGlobalRefine          = 0;
    parameters.edgeRefinements          = 2;

    getParameters() = parameters;
  }

  void readParameters() {
    // Parameters
    Utilities::Parameters parameters{};

    // Load Preferences
    std::ifstream input_Parameters;
    input_Parameters.open("auxiliaryFiles/parameters.json");

    json parameterJSON;
    try {
      parameterJSON = json::parse(input_Parameters);
      parameters    = parameterJSON.get<Utilities::Parameters>();
    } catch (json::parse_error& ex) {
      std::cerr << "parse error at byte " << ex.byte << std::endl;
    }
    getParameters() = parameters;
  }

  template <typename T>
  std::vector<T> linspace(T a, T b, unsigned int N) {
    T inc    = (b - a) / static_cast<T>(N - 1);
    auto val = [a, inc](int i) -> T { return a + i * inc; };

    std::vector<T> xs(N);
    for (auto i : std::ranges::iota_view{1u, N - 1})
      xs[i] = val(i);

    // Make sure the beginning and end are 100 interpolatory (no floating point arithmetic)
    xs.front() = a;
    xs.back()  = b;

    return xs;
  }

  template <std::floating_point T>
  double map(T value, T inputMin, T inputMax, T outputMin, T outputMax) {
    static_assert(!(std::is_integral<T>::value), "map function invoked with Integer-like datatype");
    return (value - inputMin) * (outputMax - outputMin) / (inputMax - inputMin) + outputMin;
  }

  template <std::floating_point T>
  std::array<std::array<T, 2>, 2> splitDomains(std::array<T, 2> domain) {
    T halfPoint = map<T>(0.5, 0, 1, domain[0], domain[1]);
    std::array<std::array<T, 2>, 2> result;
    result[0] = {domain[0], halfPoint};
    result[1] = {halfPoint, domain[1]};

    return result;
  }

  void printVec(auto vec, std::string&& name = "Vec") {
    std::cout << name << ": ";
    for (auto& item : vec)
      std::cout << item << " ";
    std::cout << std::endl;
  }

}  // namespace Dune::IGA::Utilities

#endif  // DUNE_IGA_NURBSTRIMUTILS_HH

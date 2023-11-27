
#pragma once
namespace Dune::IGANEW::Concept {

  template <typename T>
  concept TrimDataContainer = requires(const T& t, size_t i) {
    t[i];
  };

  template <typename T>
  concept Trimmer = TrimDataContainer<typename T::TrimDataContainer> and requires {
    typename T::ParameterSpaceGrid;
    typename T::TrimDataContainer;
    T::isAlwaysTrivial;
    // referenceElement()
  };

}  // namespace Dune::IGANEW::Concept


#pragma once
#include <concepts>
#include <optional>
namespace Dune::IGANEW::Concept {

  template <typename T>
  concept TrimDataContainer = requires(const T& t, size_t i) {
    t[i];
  };




template <typename T>
concept Trimmer = requires(T trimmer,
                                  const typename T::ParameterType& param,
                                  const typename T::ElementTrimData& elementTrimData,
                                  const typename T::PatchTrimData& patchTrimData,
                                  const typename T::ParameterSpaceGrid& paramSpaceGrid) {

    { T::mydimension } -> std::convertible_to<int>;
     typename T::ctype ;
     typename T::ParameterSpaceGrid ;
     typename T::ParameterType ;
     typename T::ElementTrimData ;
     typename T::PatchTrimData ;

    { T() } -> std::convertible_to<T>;  // Default constructor

    { T::template isLocalGeometryLinear<0> } -> std::convertible_to<bool>;
    { T::isAlwaysTrivial } -> std::convertible_to<bool>;

    // { trimmer.referenceElement(std::declval<const EntityType>()) } -> std::convertible_to<typename T::ReferenceElementType>;

    { typename T::ElementTrimDataContainer() } -> std::convertible_to<typename T::ElementTrimDataContainer>;

    // { trimmer.trimData(std::declval<const EntityType>(), std::declval<const GlobalIdSet>()) } -> std::convertible_to<std::optional<std::reference_wrapper<const typename T::ElementTrimData>>>;

    { trimmer.globalRefine(0) } -> std::convertible_to<void>;

    { trimmer.parameterSpaceGrid() } -> std::convertible_to<const typename T::ParameterSpaceGrid&>;
    { trimmer.parameterSpaceGrid() } -> std::convertible_to<typename T::ParameterSpaceGrid&>;


    { trimmer.setup(std::declval<const typename T::ParameterType>()) } -> std::convertible_to<void>;

    // { trimmer.unTrimmedParameterSpaceGrid() } -> std::convertible_to<const typename T::UntrimmedParameterSpaceGrid&>;
    // { trimmer.unTrimmedParameterSpaceGrid() } -> std::convertible_to<typename T::UntrimmedParameterSpaceGrid&>;

    // { trimmer.Grid } -> std::convertible_to<typename T::Grid>;
} and Dune::Concept::Grid<typename T::ParameterSpaceGrid>;


}  // namespace Dune::IGANEW::Concept

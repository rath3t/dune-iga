
#pragma once
#include <concepts>

#include <dune/grid/concepts.hh>

namespace Dune::IGA::Concept {

template <typename T>
concept TrimDataContainer = requires(const T& t, size_t i) { t[i]; };

/**
 * @concept ParameterSpace
 * @ingroup ParameterSpace
 * @brief Concept representing a ParameterSpace class.
 * @tparam T Type to be checked for conformity to the ParameterSpace concept.
 */
template <typename T>
concept ParameterSpace =
    T::isValid and
    requires(T trimmer, const typename T::ParameterType& param, const typename T::ElementTrimData& elementTrimData,
             const typename T::PatchTrimData& patchTrimData, const typename T::ParameterSpaceGrid& paramSpaceGrid) {
      { T::mydimension } -> std::convertible_to<int>;
      typename T::ctype;
      typename T::ParameterSpaceGrid;
      typename T::ParameterType;
      typename T::ElementTrimData;
      typename T::PatchTrimData;

      T(); ///< Check for the presence of a default constructor

      { T::template isLocalGeometryLinear<0> } -> std::convertible_to<bool>;
      { T::isAlwaysTrivial } -> std::convertible_to<bool>;

      // { parameterspace.referenceElement(std::declval<const EntityType>()) } -> std::convertible_to<typename
      // T::ReferenceElementType>; Uncomment the line above when the referenceElement method is added to the
      // ParameterSpace concept.

      // { typename T::ElementTrimDataContainer() } -> std::convertible_to<typename T::ElementTrimDataContainer>;

      // { parameterspace.trimData(std::declval<const EntityType>(), std::declval<const GlobalIdSet>()) } ->
      // std::convertible_to<std::optional<std::reference_wrapper<const typename T::ElementTrimData>>>;
      // Uncomment the line above when the trimData method is added to the ParameterSpace concept.

      { trimmer.globalRefine(0) } -> std::convertible_to<void>;

      { trimmer.parameterSpaceGrid() } -> std::convertible_to<const typename T::ParameterSpaceGrid&>;
      { trimmer.parameterSpaceGrid() } -> std::convertible_to<typename T::ParameterSpaceGrid&>;

      { trimmer.setParameters(std::declval<const typename T::ParameterType>()) } -> std::convertible_to<void>;

      // { parameterspace.unTrimmedParameterSpaceGrid() } -> std::convertible_to<const typename
      // T::UntrimmedParameterSpaceGrid&>; { parameterspace.unTrimmedParameterSpaceGrid() } ->
      // std::convertible_to<typename T::UntrimmedParameterSpaceGrid&>;
    } and
    Dune::Concept::Grid<typename T::ParameterSpaceGrid>;

} // namespace Dune::IGA::Concept

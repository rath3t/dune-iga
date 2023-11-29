
#pragma once
#include <concepts>
#include <optional>
namespace Dune::IGANEW::Concept {

  template <typename T>
  concept TrimDataContainer = requires(const T& t, size_t i) {
    t[i];
  };




/**
 * @concept Trimmer
 * @brief Concept representing a Trimmer class.
 * @tparam T Type to be checked for conformity to the Trimmer concept.
 */
template <typename T>
concept Trimmer = requires(T trimmer,
                            const typename T::ParameterType& param,
                            const typename T::ElementTrimData& elementTrimData,
                            const typename T::PatchTrimData& patchTrimData,
                            const typename T::ParameterSpaceGrid& paramSpaceGrid) {

    { T::mydimension } -> std::convertible_to<int>; ///< Check for the presence of mydimension with int convertible type.
    typename T::ctype; ///< Check for the presence of ctype type.
    typename T::ParameterSpaceGrid; ///< Check for the presence of ParameterSpaceGrid type.
    typename T::ParameterType; ///< Check for the presence of ParameterType type.
    typename T::ElementTrimData; ///< Check for the presence of ElementTrimData type.
    typename T::PatchTrimData; ///< Check for the presence of PatchTrimData type.

    { T() } -> std::convertible_to<T>; ///< Check for the presence of a default constructor returning T.

    { T::template isLocalGeometryLinear<0> } -> std::convertible_to<bool>; ///< Check for the presence of a static constexpr method isLocalGeometryLinear<0> returning bool.
    { T::isAlwaysTrivial } -> std::convertible_to<bool>; ///< Check for the presence of a static constexpr member isAlwaysTrivial returning bool.

    // { trimmer.referenceElement(std::declval<const EntityType>()) } -> std::convertible_to<typename T::ReferenceElementType>;
    // Uncomment the line above when the referenceElement method is added to the Trimmer concept.

    { typename T::ElementTrimDataContainer() } -> std::convertible_to<typename T::ElementTrimDataContainer>; ///< Check for the presence of ElementTrimDataContainer type.

    // { trimmer.trimData(std::declval<const EntityType>(), std::declval<const GlobalIdSet>()) } -> std::convertible_to<std::optional<std::reference_wrapper<const typename T::ElementTrimData>>>;
    // Uncomment the line above when the trimData method is added to the Trimmer concept.

    { trimmer.globalRefine(0) } -> std::convertible_to<void>; ///< Check for the presence of a globalRefine method returning void.

    { trimmer.parameterSpaceGrid() } -> std::convertible_to<const typename T::ParameterSpaceGrid&>; ///< Check for the presence of a const parameterSpaceGrid method returning const ParameterSpaceGrid&.
    { trimmer.parameterSpaceGrid() } -> std::convertible_to<typename T::ParameterSpaceGrid&>; ///< Check for the presence of a parameterSpaceGrid method returning ParameterSpaceGrid&.

    { trimmer.setup(std::declval<const typename T::ParameterType>()) } -> std::convertible_to<void>; ///< Check for the presence of a setup method taking const ParameterType& and returning void.

    // { trimmer.unTrimmedParameterSpaceGrid() } -> std::convertible_to<const typename T::UntrimmedParameterSpaceGrid&>;
    // { trimmer.unTrimmedParameterSpaceGrid() } -> std::convertible_to<typename T::UntrimmedParameterSpaceGrid&>;
    // Uncomment the lines above when unTrimmedParameterSpaceGrid method is added to the Trimmer concept.

    // { trimmer.Grid } -> std::convertible_to<typename T::Grid>;
    // Uncomment the line above when Grid member is added to the Trimmer concept.

} and Dune::Concept::Grid<typename T::ParameterSpaceGrid>;



}  // namespace Dune::IGANEW::Concept

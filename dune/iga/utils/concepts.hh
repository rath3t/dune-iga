// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/grid/concepts/geometry.hh>
#include <dune/iga/geometry/geohelper.hh>

namespace Dune::IGA::Concept {

  template <typename VectorType>
  concept Vector = requires(VectorType v, double a, int index) {
    typename VectorType::value_type;
    { v[index] } -> std::same_as<typename VectorType::value_type&>;
    { v += v } -> std::same_as<VectorType&>;
    { v -= v } -> std::same_as<VectorType&>;
    { v *= a } -> std::same_as<VectorType&>;
    { v /= a } -> std::same_as<VectorType&>;
  };

  /// \concept ControlPoint
  /// \tparam ControlPointType
  template <typename ControlPointType>
  concept ControlPoint = Vector<typename ControlPointType::VectorType> && requires(ControlPointType cp) {
    typename ControlPointType::VectorType;
    typename ControlPointType::VectorType::value_type;
    { cp.p } -> std::same_as<typename ControlPointType::VectorType&>;
    { cp.w } -> std::same_as<typename ControlPointType::VectorType::value_type&>;
  };

  template <typename L, typename R>
  concept MultiplyAble = requires(L x, R y) {
    x* y;
  };

  template <typename L, typename R>
  concept AddAble = requires(L x, R y) {
    x + y;
  };

  template <typename L, typename R>
  concept SubstractAble = requires(L x, R y) {
    x - y;
  };

  template <typename L, typename R>
  concept MultiplyAssignAble = requires(L x, R y) {
    x *= y;
  };

  template <typename L, typename R>
  concept DivideAssignAble = requires(L x, R y) {
    x /= y;
  };

  template <typename L, typename R>
  concept DivideAble = requires(L x, R y) {
    x / y;
  };

  /** \concept This concepts extends the geometry of DUNE grid entities
   *
   * @tparam G
   */
  template <typename G>
  concept NurbsGeometry = Dune::Concept::Geometry<G> && requires(const G g, typename G::GlobalCoordinate global,
                                                                 typename G::LocalCoordinate local) {
    {
      g.zeroFirstAndSecondDerivativeOfPosition(local)
      } -> std::convertible_to<
          std::tuple<typename G::GlobalCoordinate, typename G::JacobianTransposed, typename G::Hessian>>;
    { g.domainMidPoint() } -> std::convertible_to<typename G::LocalCoordinate>;
    { g.domain() } -> std::convertible_to<std::array<Utilities::Domain<double>, G::mydimension>>;
  };

}  // namespace Dune::IGA::Concept

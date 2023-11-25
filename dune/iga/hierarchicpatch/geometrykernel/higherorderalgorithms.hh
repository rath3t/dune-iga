// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <dune/iga/hierarchicpatch/geometrykernel/closestpointprojection.hh>
#include <dune/iga/hierarchicpatch/geometrykernel/geohelper.hh>

namespace Dune::IGANEW::GeometryKernel {

  /**
   * @brief Computes the parameter space coordinate corresponding to a given global coordinate on a geometry.
   *
   * This function calculates the parameter space coordinate for a given global coordinate on a geometry.
   * It uses different strategies based on the dimensionality of the geometry and whether it is a parametric
   * or non-parametric case.
   *
   * @tparam Geometry The geometry type representing the shape.
   * @param geo The geometry object.
   * @param global The global coordinate for which the parameter space coordinate is computed.
   * @return The parameter space coordinate corresponding to the global coordinate.
   * If the jacobian is not invertable std::numeric_limits<ctype>::max() is returned as local coordinates.
   *
   * @note The geometry type must provide appropriate member functions and traits:
   *  - `Geometry::mydimension`: The dimensionality of the parameter space.
   *  - `Geometry::worlddimension`: The dimensionality of the world space.
   *  - `Geometry::ctype`: The type for numerical computations.
   *  - `Geometry::LocalCoordinate`: The type representing a point in the parameter space.
   *  - `Geometry::GlobalCoordinate`: The type representing a point in the world space.
   *  - `Geometry::MatrixHelper`: A helper class providing matrix operations.
   *
   * @throws Dune::MathError if the trust region method fails to find an energy-decreasing direction.
   *
   */
  template <typename Geometry>
  typename Geometry::LocalCoordinate findClosestParameterSpaceCoordinate(
      const Geometry& geo, const typename Geometry::GlobalCoordinate& global) {
    static constexpr int mydimension    = Geometry::mydimension;
    static constexpr int worlddimension = Geometry::worlddimension;
    using ctype                         = typename Geometry::ctype;
    using LocalCoordinate               = typename Geometry::LocalCoordinate;
    using GlobalCoordinate              = typename Geometry::GlobalCoordinate;
    using MatrixHelper                  = typename Geometry::MatrixHelper;
    if constexpr (mydimension == 0)
      return {};
    else if constexpr (mydimension != worlddimension) {
      auto [u, Ru, fu, gap] = closestPointProjectionByTrustRegion(geo, global);
      return u;
    } else {
      const ctype tolerance = ctype(16) * std::numeric_limits<ctype>::epsilon();
      LocalCoordinate x     = geo.domainMidPoint();

      LocalCoordinate dx{};
      do {  // from multilinearGeometry
        const GlobalCoordinate dglobal = geo.global(x) - global;
        MatrixHelper::template xTRightInvA<mydimension, worlddimension>(geo.jacobianTransposed(x), dglobal, dx);
        const bool invertible
            = MatrixHelper::template xTRightInvA<mydimension, worlddimension>(geo.jacobianTransposed(x), dglobal, dx);

        if (!invertible) return LocalCoordinate(std::numeric_limits<ctype>::max());
        x -= dx;
        // if local is outside of maximum knot vector span bound, thus we clamp it to it and return
        // clamp result into boundaries
        if (Dune::IGANEW::Utilities::clampToBoundaryAndCheckIfIsAtAllBoundaries(x, geo.domain())) {
          break;
        }

      } while (dx.two_norm2() > tolerance);
      return x;
    }
  }

  /**
   * @brief Computes the Hessian matrix of a NURBS geometry at a specified local coordinate.
   *
   * The Hessian matrix represents the second-order partial derivatives of the basis functions
   * with respect to the local coordinates. The computation is based on the given local coordinate,
   * NURBS local view, and the net of control points' coordinates.
   *
   * @tparam LocalCoordinate The type representing a point in the local coordinate space.
   * @tparam NurbsLocalView The type providing access to NURBS basis functions and derivatives.
   * @tparam ControlPointCoordinateNetType The type representing the coordinates of control points.
   * @param local The local coordinate at which the Hessian is computed.
   * @param nurbsLocalView The NURBS local view providing basis function derivatives.
   * @param localControlPointNet The net of control points' coordinates.
   * @return The Hessian matrix at the specified local coordinate.
   *
   *
   * @see NurbsLocalView::basisFunctionDerivatives
   *
   */
  template <typename LocalCoordinate, typename NurbsLocalView, typename ControlPointCoordinateNetType>
  auto hessian(const LocalCoordinate& local, const NurbsLocalView& nurbsLocalView,
               const ControlPointCoordinateNetType& localControlPointNet) {
    static constexpr int mydimension    = NurbsLocalView::dimension;
    static constexpr int worlddimension = ControlPointCoordinateNetType::value_type::dimension;
    using ctype                         = typename ControlPointCoordinateNetType::value_type::value_type;

    FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, worlddimension> H;

    const auto basisFunctionDerivatives = nurbsLocalView.basisFunctionDerivatives(local, 2);

    for (int dir = 0; dir < mydimension; ++dir) {
      std::array<unsigned int, mydimension> ithVec{};
      ithVec[dir] = 2;  // second derivative in dir direction
      H[dir]      = dot(basisFunctionDerivatives.get(ithVec), localControlPointNet);
    }
    if constexpr (mydimension > 1) {
      std::array<int, mydimension> mixeDerivs;
      std::ranges::fill(mixeDerivs, 0);  // first mixed derivatives
      if constexpr (mydimension == 2)
        for (int dir = 0; dir < mydimension; ++dir) {
          mixeDerivs[dir] = 1;
        }
      else
        std::ranges::fill_n(mixeDerivs.begin() + 1, 2, 1);  // first mixed derivatives
      int mixedDireCounter = mydimension;
      do {
        H[mixedDireCounter++] = dot(basisFunctionDerivatives.get(mixeDerivs), localControlPointNet);
        if constexpr (mydimension == 2) break;
      } while (std::ranges::next_permutation(mixeDerivs, std::less()).found);
    }
    return H;
  }

  /**
   * @brief Computes the transposed Jacobian matrix of a NURBS geometry at a specified local coordinate.
   *
   * The transposed Jacobian matrix represents the first-order partial derivatives of the basis functions
   * with respect to the local coordinates, transposed for use in geometry transformations. The computation
   * is based on the given local coordinate, NURBS local view, and the net of control points' coordinates.
   *
   * @tparam LocalCoordinate The type representing a point in the local coordinate space.
   * @tparam NurbsLocalView The type providing access to NURBS basis functions and derivatives.
   * @tparam ControlPointCoordinateNetType The type representing the coordinates of control points.
   * @param local The local coordinate at which the transposed Jacobian is computed.
   * @param nurbsLocalView The NURBS local view providing basis function derivatives.
   * @param localControlPointNet The net of control points' coordinates.
   * @return The transposed Jacobian matrix at the specified local coordinate.
   */
  template <typename LocalCoordinate, typename NurbsLocalView, typename ControlPointCoordinateNetType>
  [[nodiscard]] auto jacobianTransposed(const LocalCoordinate& local, const NurbsLocalView& nurbsLocalView,
                                        const ControlPointCoordinateNetType& localControlPointNet) {
    static constexpr int mydimension    = NurbsLocalView::dimension;
    static constexpr int worlddimension = ControlPointCoordinateNetType::value_type::dimension;
    using ctype                         = typename ControlPointCoordinateNetType::value_type::value_type;

    FieldMatrix<ctype, mydimension, worlddimension> result;

    const auto basisFunctionDerivatives = nurbsLocalView.basisFunctionDerivatives(local, 1);

    for (int dir = 0; dir < mydimension; ++dir) {
      std::array<unsigned int, mydimension> ithVec{};
      ithVec[dir] = 1;
      result[dir] = dot(basisFunctionDerivatives.get(ithVec), localControlPointNet);
    }
    return result;
  }

  /**
   * @brief Computes the global position of a NURBS geometry at a specified local coordinate.
   *
   * This function calculates the global position of a NURBS geometry at a given local coordinate.
   * The computation is based on the NURBS basis functions evaluated at the local coordinate and
   * the net of control points' coordinates.
   *
   * @tparam LocalCoordinate The type representing a point in the local coordinate space.
   * @tparam NurbsLocalView The type providing access to NURBS basis functions.
   * @tparam ControlPointCoordinateNetType The type representing the coordinates of control points.
   * @param local The local coordinate at which the global position is computed.
   * @param nurbsLocalView The NURBS local view providing basis functions.
   * @param localControlPointNet The net of control points' coordinates.
   * @return The global position at the specified local coordinate.
   */
  template <typename LocalCoordinate, typename NurbsLocalView, typename ControlPointCoordinateNetType>
  [[nodiscard]] auto position(const LocalCoordinate& local, const NurbsLocalView& nurbsLocalView,
                              const ControlPointCoordinateNetType& localControlPointNet) {
    auto basis = nurbsLocalView.basisFunctions(local);
    return dot(basis, localControlPointNet);
  }

}  // namespace Dune::IGANEW::GeometryKernel

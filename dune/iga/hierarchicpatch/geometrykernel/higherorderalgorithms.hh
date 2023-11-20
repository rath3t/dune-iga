// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <dune/iga/hierarchicpatch/geometrykernel/closestpointprojection.hh>
#include <dune/iga/hierarchicpatch/geometrykernel/geohelper.hh>

namespace Dune::IGANEW::GeometryKernel {

  template <typename Geometry>
  typename Geometry::LocalCoordinate computeParameterSpaceCoordinate(
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

  template <typename LocalCoordinate, typename NurbsLocalView, typename ControlPointCoordinateNetType>
  [[nodiscard]] auto position(const LocalCoordinate& local, const NurbsLocalView& nurbsLocalView,
                              const ControlPointCoordinateNetType& localControlPointNet) {
    auto basis = nurbsLocalView.basisFunctions(local);
    return dot(basis, localControlPointNet);
  }

}  // namespace Dune::IGANEW::GeometryKernel

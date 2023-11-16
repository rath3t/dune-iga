// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once


namespace Dune::IGANEW {

  template<typename Geometry>
  typename Geometry::LocalCoordinate computeParameterSpaceCoordinate(const Geometry& geo, const typename Geometry::GlobalCoordinate& global) {
    static constexpr int mydimension = Geometry::mydimension;
    static constexpr int worlddimension = Geometry::worlddimension;
    using ctype = typename Geometry::ctype;
    using LocalCoordinate = typename Geometry::LocalCoordinate;
    using GlobalCoordinate = typename Geometry::GlobalCoordinate;
    using MatrixHelper = typename Geometry::MatrixHelper;
    if constexpr (mydimension == 0)
      return {};
    else if constexpr (mydimension != worlddimension) {
      auto [u, Ru, fu, gap] = Dune::IGA::closestPointProjectionByTrustRegion(geo, global);
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
        if (IGA::Utilities::clampToBoundaryAndCheckIfIsAtAllBoundaries(x, geo.domain())) {
          break;
        }

      } while (dx.two_norm2() > tolerance);
      return x;
    }
  }


}

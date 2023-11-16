// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once


namespace Dune::IGANEW {

  template<typename Geometry>
  typename Geometry::LocalCoordinate computeParameterSpaceCoordinate(const Geometry& geo, const typename Geometry::GlobalCoordinate& global) {
    static constexpr int mydimension = Geometry::mydimension;
    static constexpr int coorddimension = Geometry::coorddimension;
    using ctype = typename Geometry::ctype;
    using LocalCoordinate = typename Geometry::LocalCoordinate;
    using GlobalCoordinate = typename Geometry::GlobalCoordinate;
    using MatrixHelper = typename Geometry::MatrixHelper;
    if constexpr (mydimension == 0)
      return {};
    else if constexpr (mydimension != coorddimension) {
      auto [u, Ru, fu, gap] = Dune::IGA::closestPointProjectionByTrustRegion(geo, global);
      return u;
    } else {
      const ctype tolerance = ctype(16) * std::numeric_limits<ctype>::epsilon();
      LocalCoordinate x     = geo.domainMidPoint();

      LocalCoordinate dx{};
      do {  // from multilinearGeometry
        const GlobalCoordinate dglobal = geo.global(x) - global;
        MatrixHelper::template xTRightInvA<mydimension, coorddimension>(geo.jacobianTransposed(x), dglobal, dx);
        const bool invertible
            = MatrixHelper::template xTRightInvA<mydimension, coorddimension>(geo.jacobianTransposed(x), dglobal, dx);

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


      auto zeroFirstAndSecondDerivativeOfPosition(const LocalCoordinate& uL) const {
      const auto u = localToSpan(uL);
      GlobalCoordinate pos;

      Hessian H;
      std::array<unsigned int, mydimension> subDirs = getDirectionsOfSubEntityInParameterSpace();

      const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(u, 2);

      std::array<unsigned int, griddim> ithVecZero{};
      pos = Dune::IGA::dot(basisFunctionDerivatives.get(ithVecZero), cpCoordinateNet_);
      JacobianTransposed J;
      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 1;
        J[dir]               = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet_);
        J[dir] *= scaling_.at(subDirs[dir]);  // transform back to 0..1 domain
      }

      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 2;  // second derivative in subDirs[dir] direction
        H[dir]               = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet_);

        H[dir] *= Dune::power(scaling_.at(subDirs[dir]), 2);  // transform back to 0..1
      }
      if constexpr (mydimension > 1 and griddim > 1) {
        std::array<int, griddim> mixeDerivs;
        std::ranges::fill(mixeDerivs, 0);  // first mixed derivatives
        if constexpr (mydimension == 2)
          for (int dir = 0; dir < mydimension; ++dir) {
            mixeDerivs[subDirs[dir]] = 1;
          }
        else
          std::ranges::fill_n(mixeDerivs.begin(), 2, 1);  // first mixed derivatives
        int mixedDireCounter = mydimension;
        do {
          H[mixedDireCounter++] = dot(basisFunctionDerivatives.get(mixeDerivs), cpCoordinateNet_);
          for (int dir = 0; dir < mydimension; ++dir) {
            if (mixeDerivs[dir] == 0) continue;
            H[mixedDireCounter - 1] *= scaling_.at(subDirs[dir]);
          }
          if constexpr (mydimension == 2) break;

        } while (std::ranges::next_permutation(mixeDerivs, std::greater()).found);
      }

      return std::make_tuple(pos, J, H);
    }


}

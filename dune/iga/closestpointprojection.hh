//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// SPDX-FileCopyrightText: 2006  John Maddock
//
// SPDX-License-Identifier: LICENSE_1_0.txt

#pragma once

#include <utility>

#include <dune/iga/nurbspatchgeometry.h>

namespace Dune {

  namespace IGA {
    template <template <std::integral auto, std::integral auto, typename> typename Geo, std::integral auto dim,
              std::integral auto dimworld, typename GeoArgs>
    auto simpleClosestPointProjection(
        const Geo<dim, dimworld, GeoArgs>& geo, auto& point,
        const std::optional<Dune::FieldVector<typename Geo<dim, dimworld, GeoArgs>::ctype, static_cast<size_t>(dim)>>&
            start
        = std::nullopt)
        -> std::tuple<Dune::FieldVector<typename Geo<dim, dimworld, GeoArgs>::ctype, dim>,
                      typename Geo<dim, dimworld, GeoArgs>::ctype, typename Geo<dim, dimworld, GeoArgs>::ctype> requires(dim==1 or dim==2){
      using Geometry = Geo<dim, dimworld, GeoArgs>;
      using ctype    = typename Geometry::ctype;
      auto energy    = [&](auto&& uL) {
        auto pointOnCurve = geo.global(uL);
        auto dist         = pointOnCurve - point;
        return 0.5 * dist.two_norm2();
      };

      auto energyGradAndHess = [&](auto&& uL) {
        auto [pointOnCurve, J, H] = geo.zeroFirstAndSecondDerivativeOfPosition(uL);
        if constexpr (dim == 1) {
          const auto& a1                    = J[0];
          const auto& a1d1                  = H[0];
          auto dist                         = pointOnCurve - point;
          ctype dfdxi1                      = a1d1 * dist;
          ctype df2dxi11                    = a1 * a1 + dfdxi1;
          Dune::FieldMatrix<ctype, 1, 1> HV = df2dxi11;
          Dune::FieldVector<ctype, 1> g     = a1 * dist;
          return std::make_tuple(0.5 * dist.two_norm2(), g, HV);
        } else if constexpr (dim == 2) {
          const auto& a1   = J[0];
          const auto& a2   = J[1];
          const auto& a1d1 = H[0];
          const auto& a2d2 = H[1];
          const auto& a1d2 = H[2];
          auto dist        = pointOnCurve - point;
          ctype dfdxi1     = a1d1 * dist;
          ctype dfdxi2     = a2d2 * dist;
          ctype dfdxi12    = a1d2 * dist;
          ctype df2dxi11   = a1 * a1 + dfdxi1;
          ctype df2dxi22   = a2 * a2 + dfdxi2;
          ctype df2dxi12   = a1 * a2 + dfdxi12;
          Dune::FieldMatrix<ctype, 2, 2> HV;
          HV[0][0] = df2dxi11;
          HV[1][1] = df2dxi22;
          HV[0][1] = HV[1][0] = dfdxi12;
          Dune::FieldVector<ctype, 2> g;
          g[0] = a1 * dist;
          g[1] = a2 * dist;
          return std::make_tuple(0.5 * dist.two_norm2(), g, HV);
        }
      };

      auto isPositiveDefiniteInDirection = [](auto&& H, auto&& du) {
        std::remove_cvref_t<decltype(du)> tmp;
        H.mv(du,tmp);
        auto kappa = du * tmp;
        return kappa > 0.0;
      };

      int maxiter = 100;

      // start with mid-point if nothing else is given
      auto uD = start ? start.value() : geo.domainMidPoint();
      Dune::FieldVector<ctype, dim> u{};
      std::ranges::copy(uD, std::begin(u));
      ctype tol = ctype(16) * std::numeric_limits<ctype>::epsilon();
      ctype energyVal;
      ctype Rnorm;
      auto domain = geo.domain();

      double domainFractionFactor = 10;
      std::array<ctype, dim> domainFraction;
      for (int j = 0; j < dim; ++j)
        domainFraction[j] = (domain[j][1] - domain[j][0]) / domainFractionFactor;
      int i;
      for (i = 0; i < maxiter; ++i) {
        auto [oldEnergy, R, H] = energyGradAndHess(u);
        Dune::FieldVector<ctype, dim> du;
        H.solve(du, -R);
        bool isPositiveDef = isPositiveDefiniteInDirection(H, du);
        // if the increment is too large, we tend to overshoot, thus, we limit it to a fraction of the total domain
        for (int j = 0; j < dim; ++j)
          du[j] = (abs(du[j]) > domainFraction[j]) ? du[j] / abs(du[j]) * domainFraction[j] : du[j];
        auto uold = u;
        u += du;
        // clamp result into boundaries
        int breakDueToBoundaryCounter = 0;
        for (int j = 0; j < dim; ++j) {
          // if clamp does not return the old value we are at the boundary at break,
          // this we check by comparing the addresses
          if (&u[j] != &std::clamp(u[j], domain[j][0], domain[j][1])) ++breakDueToBoundaryCounter;
        }

        if (breakDueToBoundaryCounter == dim) break;

        for (int j = 0; j < 10; ++j) {
          energyVal = energy(u);
          if (energyVal > oldEnergy) {
            u = uold;
            // If our new step didn't move in an energy decreasing direction we have to modify
            if (isPositiveDef) {
              // If our model m(du)= oldEnergy+ dot(R,du) + 0.5*(du^T*H*du) has a positive curvature du^T*Hessian*du, we
              // got
              //  as du the minimizer of the parabola. Since we got energy increase, our model is bad in comparison to
              //  the real energy(u). Thus, we decrease our step size and hope for the best, at least in some (small)
              //  neighborhood our model is good
              u += du / (j + 1);
            } else {
              // If our model is not positive definite (a downward parabola), using du is useless, since it is at the
              // maximizer of m(du)
              //  (classical Newton method results, since it only finds the stationary point)
              //  Thus, we go in the opposite direction, since we know that we go downward the model m(du) and again
              //  check for real energy decrease
              u -= du / (j + 1);
            }
          } else
            break;
        }
        Rnorm = R.two_norm();
        if (Rnorm < tol) break;
      }
      //      std::cout<<"i: "<<i;
      // If we end up at the boundary, and the residual is non-zero, we restart at the opposite domain boundary
      // If this does also not help we return the point at the boundary with smaller energy(distance)
      for (int j = 0; j < dim; ++j)
        u[j] = std::clamp(u[j], domain[j][0], domain[j][1]);
      auto uRestart = u;
      if (Rnorm > tol and not start) {
        for (int j = 0; j < dim; ++j)
          if (abs(domain[j][0] - uRestart[j]) < tol)
            uRestart[j] = domain[j][1];
          else if (abs(domain[j][1] - uRestart[j]) < tol)
            uRestart[j] = domain[j][0];
        decltype(uD) uARestart{};
        std::ranges::copy(uRestart, std::begin(uARestart));
        //        std::cout << "lowUI " << domain[0][0] << " upperUI " << domain[0][1] << std::endl;
        //        std::cout<<"uRestart "<<uRestart<<std::endl;
        auto [u2, Rnorm2, energyVal2] = simpleClosestPointProjection(geo, point, uARestart);
        //        std::cout<<" u2 "<<u2[0]<<" Rnorm2 "<<Rnorm2<<" energyVal2 "<<energyVal2<<std::endl;
        //        std::cout<<" u "<<u<<" Rnorm "<<Rnorm<<" energyVal "<<energyVal<<std::endl;

        return energyVal2 > energyVal ? std::make_tuple(u, Rnorm, energyVal) : std::make_tuple(u2, Rnorm2, energyVal2);
      }

      return std::make_tuple(u, Rnorm, energyVal);
    }
  }  // namespace IGA

}  // namespace Dune

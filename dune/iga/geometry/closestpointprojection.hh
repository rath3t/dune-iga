// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/geometry/geohelper.hh"
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/transpose.hh>

namespace Dune::IGA {
  template <template <std::integral auto, std::integral auto, typename> typename Geo, std::integral auto dim,
            std::integral auto dimworld, typename GeoArgs>
  auto closestPointProjectionByTrustRegion(
      const Geo<dim, dimworld, GeoArgs>& geo, auto& point,
      const std::optional<Dune::FieldVector<typename Geo<dim, dimworld, GeoArgs>::ctype, static_cast<size_t>(dim)>>&
          start
      = std::nullopt)
      -> std::tuple<Dune::FieldVector<typename Geo<dim, dimworld, GeoArgs>::ctype, dim>,
                    typename Geo<dim, dimworld, GeoArgs>::ctype, typename Geo<dim, dimworld, GeoArgs>::ctype,
                    typename Geo<dim, dimworld, GeoArgs>::ctype>
  requires(dim == 1 or dim == 2) {
    using Geometry = Geo<dim, dimworld, GeoArgs>;
    using ctype    = typename Geometry::ctype;

    auto energy = [&](const auto& uL) {
      auto pointOnCurve = geo.global(uL);
      auto dist         = pointOnCurve - point;
      return 0.5 * dist.two_norm2();
    };

    auto energyGradAndHess = [&](const auto& uL) {
      auto [pointOnCurve, J, H] = geo.zeroFirstAndSecondDerivativeOfPosition(uL);
      auto dist                 = pointOnCurve - point;
      Dune::FieldVector<ctype, dim> g;
      J.mv(dist, g);

      Dune::FieldMatrix<ctype, dim*(dim + 1) / 2, 1> dfdxideta;
      H.mv(dist, dfdxideta);

      Dune::FieldMatrix<ctype, dim, dim> HV = J * Dune::transpose(J);
      for (int i = 0; i < dim; ++i)
        HV[i][i] += dfdxideta[i];

      if constexpr (dim == 2) {
        HV[0][1] += dfdxideta[2];
        HV[1][0] += dfdxideta[2];
      }

      return std::make_tuple(0.5 * dist.two_norm2(), g, HV);
    };

    auto isPositiveDefiniteInDirection = [](auto&& H, auto&& du) {
      std::remove_cvref_t<decltype(du)> tmp;
      H.mv(du, tmp);
      auto kappa = du * tmp;
      return kappa > 0.0;
    };

    int maxiter = 100;

    // start with mid-point if nothing else is given
    Dune::FieldVector<ctype, dim> u = start ? start.value() : geo.domainMidPoint();

    ctype tol = 1e-10;
    ctype energyVal;
    ctype Rnorm;
    auto domain = geo.domain();

    double domainFractionFactor = 4;
    ctype domainFraction        = 1;
    for (int j = 0; j < dim; ++j)
      domainFraction *= domain[j].size();
    domainFraction /= domainFractionFactor;
    int i;
    auto [oldEnergy, R, H] = energyGradAndHess(u);
    energyVal              = oldEnergy;
    for (i = 0; i < maxiter; ++i) {
      Dune::FieldVector<ctype, dim> du;
      H.solve(du, -R);
      Rnorm = R.two_norm();

      // if the increment is too large, we tend to overshoot, thus, we limit it to a fraction of the total domain
      ctype duNorm = du.two_norm();
      if (duNorm > domainFraction) du *= domainFraction / duNorm;

      const bool isPositiveDef = isPositiveDefiniteInDirection(H, du);
      auto uold                = u;
      oldEnergy                = energy(u);
      u += du;
      if (Utilities::clampToBoundaryAndCheckIfIsAtAllBoundaries(u, domain)) break;

      int testSteps = 20;
      for (int j = 0; j < testSteps; ++j) {
        energyVal = energy(u);
        if (energyVal > oldEnergy and duNorm > 1e-8) {  // we directly accept very small steps
          u                 = uold;
          ctype scaleFactor = Dune::power(2.0, j + 1);
          // If our new step didn't move in an energy-decreasing direction, we have to modify
          if (isPositiveDef) {
            // If our model m(du)= oldEnergy+ dot(R,du) + 0.5*(du^T*H*du) has a positive curvature du^T*Hessian*du, we
            // got
            //  as du the minimizer of the parabola. Since we got energy increase, our model is bad in comparison to
            //  the real energy(u). Thus, we decrease our step size and hope for the best, at least in some (small)
            //  neighborhood our model is good
            u += du / scaleFactor;
          } else {
            // If our model is not positive definite (a downward parabola), using du is useless, since it is at the
            // maximizer of m(du)
            //  (classical Newton method results, since it only finds the stationary point)
            //  Thus, we go in the opposite direction, since we know that we go downward the model m(du) and again
            //  check for real energy decrease
            u -= du / scaleFactor;
          }
        } else
          break;

        if (j == testSteps - 1) DUNE_THROW(Dune::MathError, "We couldn't find an energy decreasing direction");
      }
      std::tie(oldEnergy, R, H) = energyGradAndHess(u);

      Rnorm = R.two_norm();
      if (Rnorm < tol) break;
    }
    // If we end up at the boundary, and the residual is non-zero, we restart at the opposite domain boundary
    // If this does also not help we return the point at the boundary with smaller energy(distance)
    for (int j = 0; j < dim; ++j)
      u[j] = clampToDomain(u[j], domain[j]);
    //    auto uRestart = u;
    //    if (Rnorm > tol and not start) {
    //      for (int j = 0; j < dim; ++j)
    //        if (abs(domain[j].left() - uRestart[j]) < tol)
    //          uRestart[j] = domain[j].right();
    //        else if (abs(domain[j].right() - uRestart[j]) < tol)
    //          uRestart[j] = domain[j].left();
    //
    //      auto [u2, Rnorm2, energyVal2, gap] = closestPointProjectionByTrustRegion(geo, point, uRestart);
    //
    //      return energyVal2 > energyVal ? std::make_tuple(u, Rnorm, energyVal, (geo.global(u) - point).two_norm())
    //                                    : std::make_tuple(u2, Rnorm2, energyVal2, (geo.global(u2) -
    //                                    point).two_norm());
    //    }
    //    std::cout << "u " << u << " Rnorm " << Rnorm << " energyVal " << energyVal << "geo " << geo.global(u) <<
    //    std::endl;
    return std::make_tuple(u, Rnorm, energyVal, (geo.global(u) - point).two_norm());
  }
}  // namespace Dune::IGA

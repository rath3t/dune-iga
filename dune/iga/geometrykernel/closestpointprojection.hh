// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/geometrykernel/geohelper.hh"
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/transpose.hh>

namespace Dune::IGA {
/**
 * @brief Finds the closest point projection on a NURBS curve or surface by trust region method.
 *
 * This function computes the closest point projection of a given point onto a NURBS curve or surface.
 * The trust region method is used to iteratively refine the solution. The function returns a tuple
 * containing the local coordinates of the closest point projection, the residual norm, the energy value,
 * and the distance between the projected point and the original point.
 *
 * @tparam Geo The NURBS geometry type
 * @param geo The NURBS geometry object.
 * @param point The point to be projected onto the NURBS geometry.
 * @param start Optional starting point for the iteration. If not provided, the midpoint of the domain is used.
 * @return A tuple containing the local coordinates, residual norm, energy value, and distance from the point.
 *
 * @throws Dune::MathError if the trust region method fails to find an energy-decreasing direction.
 *
 */
template <typename Geo> // @todo add concept, HJ: maybe add custom return type instead of tuple
auto closestPointProjectionByTrustRegion(const Geo& geo, auto& point,
                                         const std::optional<typename Geo::LocalCoordinate>& start = std::nullopt)
    -> std::tuple<typename Geo::LocalCoordinate, typename Geo::ctype, typename Geo::ctype, typename Geo::ctype>
requires(Geo::mydimension == 1 or Geo::mydimension == 2)
{
  using ctype              = typename Geo::ctype;
  static constexpr int dim = Geo::mydimension;
  auto energy              = [&](const auto& uL) {
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

  ctype tol = ctype(16) * std::numeric_limits<ctype>::epsilon();
  ctype Rnorm;
  auto domain = geo.domain();

  double domainFractionFactor = 4;
  ctype domainFraction        = 1;
  for (int j = 0; j < dim; ++j)
    domainFraction *= domain[j].volume();
  domainFraction /= domainFractionFactor;
  int i;
  auto [oldEnergy, R, H] = energyGradAndHess(u);
  ctype energyVal        = oldEnergy;
  ctype r                = oldEnergy;
  for (i = 0; i < maxiter; ++i) {
    FieldVector<ctype, dim> du;
    H.solve(du, -R);
    // if the increment is too large, we tend to overshoot, thus, we limit it to a fraction of the total domain
    ctype duNorm = du.two_norm();
    if (duNorm > domainFraction)
      du *= domainFraction / duNorm;

    const bool isPositiveDef = isPositiveDefiniteInDirection(H, du);
    auto uold                = u;
    oldEnergy                = energy(u);
    u += du;
    if (Utilities::clampToBoundaryAndCheckIfIsAtAllBoundaries(u, domain))
      break;

    int testSteps = 20;
    for (int j = 0; j < testSteps; ++j) {
      energyVal = energy(u);
      if (energyVal > oldEnergy and duNorm > 1e-8) { // we directly accept very small steps
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

      if (j == testSteps - 1)
        DUNE_THROW(Dune::MathError, "We couldn't find an energy decreasing direction");
    }
    std::tie(oldEnergy, R, H) = energyGradAndHess(u);

    Rnorm = R.two_norm();
    if (Rnorm < tol)
      break;
  }

  // clamp again to boundary
  for (int j = 0; j < dim; ++j)
    u[j] = clampToDomain(u[j], domain[j]);

  return std::make_tuple(u, Rnorm, energyVal, (geo.global(u) - point).two_norm());
}
} // namespace Dune::IGA

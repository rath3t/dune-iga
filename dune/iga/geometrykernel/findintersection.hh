// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#define DUNE_FMatrix_WITH_CHECKING

#include "dune/common/float_cmp.hh"
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune::IGA {

enum class IntersectionCurveAndLine
{
  disjoint,
  intersect,
  parallel
};

template <typename GeoCurve, typename ScalarType, int dim>
requires(dim == 2)
auto findIntersectionLinearCurveAndLine(const GeoCurve& geoCurve, const FieldVector<ScalarType, dim>& pos,
                                        const FieldVector<ScalarType, dim>& dir, FieldVector<ScalarType, 2> tParameter,
                                        ScalarType tol = 1e-10) {
  assert(geoCurve.degree()[0] == 1 && "Curve has to be linear to be calling this function");
  assert(geoCurve.patchData().controlPoints.size() == 2 && "The linear curve should only have two control points");

  const auto& curveP0 = geoCurve.patchData().controlPoints.directGet(0).p;
  const auto& curveP1 = geoCurve.patchData().controlPoints.directGet(1).p;

  const FieldMatrix<ScalarType, dim, 2> linearSystem({curveP1 - curveP0, -dir});
  const FieldVector<ScalarType, dim> b = pos - curveP0;

  if (not FloatCmp::eq(linearSystem.determinant(), 0.0, tol)) {
    FieldVector<ScalarType, dim> sol;
    transpose(linearSystem).solve(sol, b);
    // Map the solution to the knotSpan, thus [0,1] to [domain.left(), domain.right()]
    sol[0] = mapToRangeFromZeroToOne(sol[0], geoCurve.domain()[0]);

    // Check if point is in domain
    if (not geoCurve.domain()[0].checkInside(sol[0]))
      return std::make_tuple(IntersectionCurveAndLine::disjoint, sol, curveP0);
    return std::make_tuple(IntersectionCurveAndLine::intersect, sol, geoCurve.global(sol[0]));
  }
  // If system is not solvable, the lines are paralell and therefore have no intersection
  return std::make_tuple(IntersectionCurveAndLine::parallel, tParameter, curveP0);
}

template <typename GeoCurve, typename ScalarType, int dim>
auto findIntersectionCurveAndLine(const GeoCurve& geoCurve, const FieldVector<ScalarType, dim>& pos,
                                  const FieldVector<ScalarType, dim>& dir, FieldVector<ScalarType, 2> tParameter,
                                  ScalarType tol = 1e-10) {
  if (geoCurve.degree()[0] == 1)
    return findIntersectionLinearCurveAndLine(geoCurve, pos, dir, tParameter, tol);

  struct Line
  {
    FieldVector<ScalarType, dim> operator()(ScalarType t) {
      return pos + t * dir;
    }
    const FieldVector<ScalarType, dim>& jacobian() {
      return dir;
    }

    const FieldVector<ScalarType, dim>& pos;
    const FieldVector<ScalarType, dim>& dir;
  };

  auto line = Line({pos, dir});

  // Get some trivial cases
  if (FloatCmp::eq(geoCurve.global({tParameter[0]}), geoCurve.corner(1), tol)) {
    return std::make_tuple(IntersectionCurveAndLine::intersect, tParameter, geoCurve.corner(1));
  }

  tParameter[0] = std::clamp(tParameter[0], geoCurve.domain()[0].front(), geoCurve.domain()[0].back());

  bool sucess             = false;
  const int maxIterations = 1000;
  FieldVector<ScalarType, 2> curvePoint;
  auto domain = geoCurve.domain();

  const FieldVector<ScalarType, dim>& lineDerivative = line.jacobian();
  for (int iter = 0; iter < maxIterations; ++iter) {
    curvePoint                              = geoCurve.global(tParameter[0]);
    const FieldVector curveDerivative       = geoCurve.jacobianTransposed(tParameter[0])[0];
    const FieldVector curvSecondeDerivative = geoCurve.hessian(tParameter[0])[0];
    const FieldVector dist                  = curvePoint - line(tParameter[1]);

    const ScalarType energy = 0.5 * dist.two_norm2();

    FieldVector<ScalarType, 2> grad({curveDerivative * dist, -lineDerivative * dist});

    FieldMatrix<ScalarType, 2> Hess;
    Hess[0][0] = curveDerivative.two_norm2() + dist * curvSecondeDerivative;
    Hess[1][1] = lineDerivative.two_norm2();
    Hess[0][1] = Hess[1][0] = -curveDerivative * lineDerivative;

    // Newton-Raphson update
    FieldVector<ScalarType, 2> deltaT;
    Hess.solve(deltaT, grad);
    tParameter -= deltaT;

    // Clamp
    tParameter[0] = std::clamp(tParameter[0], domain[0].front(), domain[0].back());

    // Check for convergence
    if (deltaT.two_norm() < tol) {
      sucess     = true;
      curvePoint = geoCurve.global(tParameter[0]);
      break;
    }
  }

  if (not domain[0].checkInside(tParameter[0]))
    sucess = false;
  return std::make_tuple(sucess ? IntersectionCurveAndLine::intersect : IntersectionCurveAndLine::disjoint, tParameter,
                         curvePoint);
}
} // namespace Dune::IGA

#undef DUNE_FMatrix_WITH_CHECKING

// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
namespace Dune::IGANEW {

  enum class FindIntersectionCurveAndLineResult {
    noSucess, sucess, linesParallel
  };

  template <typename GeoCurve, typename ScalarType, int dim>
  auto findIntersectionLinearCurveAndLine(const GeoCurve& geoCurve, const FieldVector<ScalarType, dim>& pos,
                                    const FieldVector<ScalarType, dim>& dir, FieldVector<ScalarType, 2> tParameter,
                                    ScalarType tol = 1e-10) {
    assert(geoCurve.degree()[0] == 1 && "Curve has to be linear to be calling this function");

    const FieldVector curveDerivative       = geoCurve.jacobianTransposed(tParameter[0])[0];
    const FieldVector curvePos = geoCurve.global({geoCurve.domain()[0].front()});

    FieldMatrix<ScalarType, dim, 2> linearSystem;
    FieldVector<ScalarType, dim> b;

    for (int i = 0; i < dim; i++) {
      linearSystem[i][0] = curveDerivative[i];
      linearSystem[i][1] = -dir[i];
      b[i] = pos[i] - curvePos[i];
    }

    // If system is not solvable, the lines are paralell and therefore have no intersection
    if (FloatCmp::eq(linearSystem.determinant(), 0.0, tol))
      return std::make_tuple(FindIntersectionCurveAndLineResult::linesParallel, tParameter, geoCurve.global(tParameter[0]));

    FieldVector<ScalarType, dim> sol;
    linearSystem.solve(sol, b);

    // Check if point is in domain
    if (not geoCurve.domain()[0].checkInside(tParameter[0]))
      return std::make_tuple(FindIntersectionCurveAndLineResult::noSucess, sol, geoCurve.global(sol[0]));
    return std::make_tuple(FindIntersectionCurveAndLineResult::sucess, sol, geoCurve.global(sol[0]));
  }

  template <typename GeoCurve, typename ScalarType, int dim>
  auto findIntersectionCurveAndLine(const GeoCurve& geoCurve, const FieldVector<ScalarType, dim>& pos,
                                    const FieldVector<ScalarType, dim>& dir, FieldVector<ScalarType, 2> tParameter,
                                    ScalarType tol = 1e-10) {
    if (geoCurve.degree()[0] == 1)
      return findIntersectionLinearCurveAndLine(geoCurve, pos, dir, tParameter, tol);

    struct Line {
      FieldVector<ScalarType, dim> operator()(ScalarType t) { return pos + t * dir; }
      FieldVector<ScalarType, dim> jacobian() { return dir; }

      const FieldVector<ScalarType, dim>& pos;
      const FieldVector<ScalarType, dim>& dir;
    };

    auto line = Line({pos, dir});

    bool sucess             = false;
    const int maxIterations = 1000;
    FieldVector<ScalarType, 2> curvePoint;

    const FieldVector lineDerivative = line.jacobian();
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

      // Check for convergence
      if (deltaT.two_norm() < tol) {
        sucess     = true;
        curvePoint = geoCurve.global(tParameter[0]);
        break;
      }
    }
    auto domain = geoCurve.domain();
    if (not domain[0].checkInside(tParameter[0])) sucess = false;
    return std::make_tuple(sucess ? FindIntersectionCurveAndLineResult::sucess : FindIntersectionCurveAndLineResult::noSucess, tParameter, curvePoint);
  }
}  // namespace Dune::IGANEW

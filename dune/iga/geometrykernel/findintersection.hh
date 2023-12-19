// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
namespace Dune::IGANEW {

  template <typename GeoCurve, typename ScalarType, int dim>
  auto findIntersectionCurveAndLine(const GeoCurve& geoCurve, const FieldVector<ScalarType, dim>& pos,
                                    const FieldVector<ScalarType, dim>& dir, FieldVector<ScalarType, dim> tParameter,
                                    ScalarType tol = 1e-10) {
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
    return std::make_tuple(sucess, tParameter, curvePoint);
  }
}  // namespace Dune::IGANEW

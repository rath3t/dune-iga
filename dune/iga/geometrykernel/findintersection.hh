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
    const auto domain = geoCurve.domain();
    FieldVector<ScalarType, 2> curvePoint;

    const FieldVector lineDerivative = line.jacobian();
    const auto energyLambda = [&](auto& tParameter_ ) {
      curvePoint                              = geoCurve.global(tParameter_[0]);
      const FieldVector dist                  = curvePoint - line(tParameter_[1]);
  return  0.5 * dist.two_norm2();
    };

    for (int iter = 0; iter < maxIterations; ++iter) {
      const FieldVector curveDerivative       = geoCurve.jacobianTransposed(tParameter[0])[0];
      const FieldVector curvSecondeDerivative = geoCurve.hessian(tParameter[0])[0];
      const FieldVector dist                  = curvePoint - line(tParameter[1]);

      const ScalarType energy =energyLambda(tParameter);

      FieldVector<ScalarType, 2> grad({curveDerivative * dist, -lineDerivative * dist});

      FieldMatrix<ScalarType, 2> Hess;
      Hess[0][0] = curveDerivative.two_norm2() + dist * curvSecondeDerivative;
      Hess[1][1] = lineDerivative.two_norm2();
      Hess[0][1] = Hess[1][0] = -curveDerivative * lineDerivative;

      // Newton-Raphson update
      FieldVector<ScalarType, 2> deltaT;
      Hess.solve(deltaT, grad);
      tParameter -= deltaT;
      tParameter[0]= std::clamp(tParameter[0],domain[0].left(),domain[0].right());

        if (tParameter[0] == domain[0].left() or tParameter[0] == domain[0].right()) {
          if(energyLambda(tParameter)<tol)
            return std::make_tuple(true, tParameter, geoCurve.global(tParameter[0]));
          else
            return std::make_tuple(false, tParameter, geoCurve.global(tParameter[0]));
        }


      // Check for convergence
      if (deltaT.two_norm() < tol) {
        sucess     = true;
        curvePoint = geoCurve.global(tParameter[0]);
        break;
      }
    }

    return std::make_tuple(sucess, tParameter, curvePoint);
  }
}  // namespace Dune::IGANEW

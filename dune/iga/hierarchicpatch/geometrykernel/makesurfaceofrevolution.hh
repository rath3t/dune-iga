

#pragma once
#include<dune/iga/hierarchicpatch/splines/nurbspatchdata.hh>
#include<dune/iga/hierarchicpatch/geometrykernel/basealgorithms.hh>
namespace Dune::IGANEW
{

  /* \brief Create a surface of revolution from a generating curve, from Piegl and Tiller, Algo 8.1
  * \tparam ScalarType the field type (use float, double, complex, etc)
  * \param generatrix 1D nurbs curve, which should be revolved
  * \param point the origin of revolution
  * \param revolutionAxis the axis of revolution
  * \param revolutionAngle how far the curve should be revolved in degrees
  * \return NURBSPatchData representing the surface
  */
  template <typename ScalarType = double>
  auto makeSurfaceOfRevolution(const NURBSPatchData<1, 3, ScalarType>& generatrix,
                               const Dune::FieldVector<ScalarType, 3>& point,
                               const Dune::FieldVector<ScalarType, 3>& revolutionAxis,
                               const ScalarType& revolutionAngle = 360.0) {
    const auto& genCP          = generatrix.controlPoints;
    using ControlPoint         = typename NURBSPatchData<2, 3, ScalarType>::ControlPointType;
    using GlobalCoordinateType = typename NURBSPatchData<2, 3, ScalarType>::GlobalCoordinateType;
    const auto pi              = std::numbers::pi_v<ScalarType>;

    const auto revolutionAxisN = revolutionAxis / revolutionAxis.two_norm();
    auto newKnotsArray        = std::array<std::vector<double>, 2UL>();
    newKnotsArray[1]          = generatrix.knotSpans[0];
    auto& U                   = newKnotsArray[0];

    const int narcs = std::ceil(revolutionAngle / 90);
    U.resize(2 * (narcs + 2));
    if (revolutionAngle <= 90.0) {
    } else if (revolutionAngle <= 180.0) {
      U[3] = U[4] = 0.5;
    } else if (revolutionAngle <= 270.0) {
      U[3] = U[4] = 1.0 / 3.0;
      U[5] = U[6] = 2.0 / 3.0;
    } else {
      U[3] = U[4] = 0.25;
      U[5] = U[6] = 0.5;
      U[7] = U[8] = 0.75;
    }
    const ScalarType dtheta = revolutionAngle / narcs * pi / 180;
    std::ranges::fill_n(U.begin(), 3, 0.0);
    std::ranges::fill_n(std::ranges::reverse_view(U).begin(), 3, 1.0);

    typename NURBSPatchData<2UL, 3UL, ScalarType>::ControlPointNetType surfaceCP(
        2 * narcs + 1, generatrix.controlPoints.strideSizes()[0]);
    using std::cos;
    using std::sin;
    const ScalarType wm = cos(dtheta / 2.0);
    std::vector<ScalarType> cosines(narcs);
    std::vector<ScalarType> sines(narcs);
    ScalarType angle = 0.0;
    for (int i = 0; i < narcs; i++) {
      angle += dtheta;
      cosines[i] = cos(angle);
      sines[i]   = sin(angle);
    }
    ControlPoint PO = genCP.directGet(0);
    for (int j = 0; j < genCP.strideSizes()[0]; j++) {
      const GlobalCoordinateType Om = projectPointOntoLine(point, revolutionAxisN, genCP.directGet(j).p);
      GlobalCoordinateType X        = genCP.directGet(j).p - Om;
      const ScalarType r            = X.two_norm();
      X /= r;
      const GlobalCoordinateType Y = cross(revolutionAxisN, X);
      surfaceCP(0, j) = PO    = genCP.directGet(j);
      GlobalCoordinateType TO = Y;
      for (int index = 0, i = 0; i < narcs; ++i) {
        const GlobalCoordinateType P2 = Om + r * cosines[i] * X + r * sines[i] * Y;
        surfaceCP(index + 2, j)       = {.p = P2, .w = genCP.directGet(j).w};

        const GlobalCoordinateType T2 = -sines[i] * X + cosines[i] * Y;
        surfaceCP(index + 1, j).p     = intersect3DLines(PO.p, TO, P2, T2);
        surfaceCP(index + 1, j).w     = wm * genCP.directGet(j).w;
        index += 2;
        if (i < narcs - 1) {
          PO.p = P2;
          TO   = T2;
        }
      }
    }
    return NURBSPatchData<2UL, 3UL, ScalarType>(newKnotsArray, surfaceCP, {2, generatrix.degree[0]});
  }
  }

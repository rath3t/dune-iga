// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <dune/iga/hierarchicpatch/geometrykernel/basealgorithms.hh>
#include <dune/iga/hierarchicpatch/splines/nurbspatchdata.hh>

namespace Dune::IGANEW {
  // TODO Alex quote algo from NURBS book
  /**
   * @brief Create a circular arc as a NURBS patch with customizable options.
   *
   * This function generates a circular arc with specified customization options, including radius, start and end
   * angles, origin, and base vectors of the plane where the arc should reside.
   *
   * @tparam ScalarType The field type for coordinates (e.g., float, double, complex).
   * @param radius Radius of the circular arc.
   * @param startAngle Starting angle of the arc in degrees.
   * @param endAngle End angle of the arc in degrees.
   * @param origin Center of the circle.
   * @param X First base vector of the plane where the arc should reside.
   * @param Y Second base vector of the plane where the arc should reside.
   * @return NURBSPatchData representing the circular arc.
   *
   * @note The angles are specified in degrees.
   *
   * Example:
   * @code
   * using ScalarType = double; // or any other numeric type
   * auto circularArc = makeCircularArc<double>(2.0, 45.0, 180.0);
   * @endcode
   */
  template <typename ScalarType = double>
  NURBSPatchData<1, 3, ScalarType> makeCircularArc(const ScalarType radius = 1.0, const ScalarType startAngle = 0.0,
                                                   ScalarType endAngle                            = 360.0,
                                                   const Dune::FieldVector<ScalarType, 3>& origin = {0, 0, 0},
                                                   const Dune::FieldVector<ScalarType, 3>& X      = {1, 0, 0},
                                                   const Dune::FieldVector<ScalarType, 3>& Y      = {0, 1, 0}) {
    using GlobalCoordinateType = typename NURBSPatchData<1, 3, ScalarType>::GlobalCoordinateType;
    const auto pi              = std::numbers::pi_v<ScalarType>;

    if (endAngle < startAngle) endAngle += 360.0;
    const ScalarType theta = endAngle - startAngle;
    const int narcs        = std::ceil(theta / 90);

    typename NURBSPatchData<1, 3, ScalarType>::ControlPointNetType circleCPs(2 * narcs + 1);
    const ScalarType dtheta  = theta / narcs * pi / 180;
    const int n              = 2 * narcs;
    const ScalarType w1      = cos(dtheta / 2.0);
    GlobalCoordinateType PO  = origin + radius * cos(startAngle) * X + radius * sin(startAngle) * Y;
    GlobalCoordinateType TO  = -sin(startAngle) * X + cos(startAngle) * Y;
    circleCPs.directGet(0).p = PO;
    ScalarType angle         = startAngle;
    for (int index = 0, i = 0; i < narcs; ++i) {
      angle += dtheta;
      const GlobalCoordinateType P2    = origin + radius * cos(angle) * X + radius * sin(angle) * Y;
      circleCPs.directGet(index + 2).p = P2;
      const GlobalCoordinateType T2    = -sin(angle) * X + cos(angle) * Y;
      const GlobalCoordinateType P1    = intersect3DLines(PO, TO, P2, T2);
      circleCPs.directGet(index + 1)   = {.p = P1, .w = w1};
      index += 2;
      if (i < narcs - 1) {
        PO = P2;
        TO = T2;
      }
    }

    auto knotVec = std::array<std::vector<double>, 1>{};
    auto& U      = knotVec[0];
    U.resize(2 * (narcs + 2));
    if (narcs == 1) {
    } else if (narcs == 2) {
      U[3] = U[4] = 1.0 / 2.0;
    } else if (narcs == 3) {
      U[3] = U[4] = 1.0 / 3.0;
      U[5] = U[6] = 2.0 / 3.0;
    } else {
      U[3] = U[4] = 1.0 / 4.0;
      U[5] = U[6] = 2.0 / 4.0;
      U[7] = U[8] = 3.0 / 4.0;
    }

    std::ranges::fill_n(U.begin(), 3, 0.0);
    std::ranges::fill_n(std::ranges::reverse_view(U).begin(), 3, 1.0);
    return NURBSPatchData<1, 3, ScalarType>(knotVec, circleCPs, {2});
  }
}  // namespace Dune::IGANEW

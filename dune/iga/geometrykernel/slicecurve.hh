// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/splines/bsplinealgorithms.hh>
#include <dune/iga/splines/nurbsalgorithms.hh>

namespace Dune::IGANEW {

  // Curtesy to https://github.com/pradeep-pyro/tinynurbs/blob/master/include/tinynurbs/core/modify.h
  template <typename GeoCurve, typename ScalarType = typename GeoCurve::ctype>
  auto splitCurve(const GeoCurve& geoCurve, ScalarType u) -> std::pair<GeoCurve, GeoCurve> {
    static_assert(GeoCurve::mydimension == 1);

    auto knotMultiplicity = [](auto& knotVec, auto u_) {
      return std::ranges::count_if(knotVec, [=](double val) { return FloatCmp::eq(u_, val, 1e-8); });
    };

    auto& patchData     = geoCurve.patchData();
    auto knots          = patchData.knotSpans.front();
    auto control_points = patchData.controlPoints.directGetAll();
    auto domain         = geoCurve.domain();
    auto degree         = patchData.degree[0];

    auto span = Splines::findSpan(degree, u, knots);
    int r     = degree - knotMultiplicity(knots, u);

    std::vector<ScalarType> newKnots(r);
    std::ranges::fill(newKnots, u);
    auto tmpPatchData = Splines::knotRefinement(patchData, newKnots, 0);

    std::vector<ScalarType> tmpKnots{tmpPatchData.knotSpans.front()};
    decltype(control_points) tmpCp{tmpPatchData.controlPoints.directGetAll()};

    std::vector<ScalarType> leftKnots;
    const int span_l = Splines::findSpan(degree, u, tmpKnots) + 1;
    for (int i = 0; i < span_l; ++i)
      leftKnots.push_back(tmpKnots[i]);
    leftKnots.push_back(u);

    std::vector<ScalarType> rightKnots;
    for (int i = 0; i < degree + 1; ++i)
      rightKnots.push_back(u);
    for (int i = span_l; i < tmpKnots.size(); ++i)
      rightKnots.push_back(tmpKnots[i]);

    decltype(control_points) leftControlPoints;
    const int ks = span - degree + 1;
    for (int i = 0; i < ks + r; ++i)
      leftControlPoints.push_back(tmpCp[i]);

    decltype(control_points) rightControlPoints;
    for (int i = ks + r - 1; i < tmpCp.size(); ++i)
      rightControlPoints.push_back(tmpCp[i]);

    typename GeoCurve::PatchData leftPatchData{
        {leftKnots}, typename GeoCurve::ControlPointNetType{leftControlPoints}, patchData.degree};
    typename GeoCurve::PatchData rightPatchData{
        {rightKnots}, typename GeoCurve::ControlPointNetType{rightControlPoints}, patchData.degree};

    return std::make_pair(GeoCurve{leftPatchData}, GeoCurve{rightPatchData});
  }

  template <typename GeoCurve, typename ScalarType = typename GeoCurve::ctype>
  auto sliceCurve(const GeoCurve& geoCurve, std::array<ScalarType, 2> t) -> GeoCurve {
    static_assert(GeoCurve::mydimension == 1);

    auto domain = geoCurve.domain();
    if (domain[0].front() == t[0] and domain[0].back() == t[1]) return geoCurve;
    if (domain[0].front() == t[0]) return std::get<0>(splitCurve(geoCurve, t[1]));
    if (domain[0].back() == t[1]) return std::get<1>(splitCurve(geoCurve, t[0]));

    auto tmpCurve = std::get<1>(splitCurve(geoCurve, t[0]));
    return std::get<0>(splitCurve(tmpCurve, t[1]));
  }

}  // namespace Dune::IGANEW
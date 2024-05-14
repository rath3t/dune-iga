// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/splines/bsplinealgorithms.hh>
#include <dune/iga/splines/nurbsalgorithms.hh>

namespace Dune::IGA {

// Inspired by https://github.com/pradeep-pyro/tinynurbs/blob/master/include/tinynurbs/core/modify.h
template <typename GeoCurve, typename ScalarType = typename GeoCurve::ctype>
requires(GeoCurve::mydimension == 1)
auto splitCurve(const GeoCurve& geoCurve, ScalarType u) -> std::pair<GeoCurve, GeoCurve> {
  using ControlPointType    = typename GeoCurve::ControlPointType;
  using PatchData           = typename GeoCurve::PatchData;
  using ControlPointNetType = typename GeoCurve::ControlPointNetType;

  const PatchData& patchData = geoCurve.patchData();
  auto knots                 = patchData.knotSpans.front();
  auto domain                = geoCurve.domain();
  auto degree                = patchData.degree[0];

  auto [multiplicity, span] = Splines::multiplicity(knots, u);
  const int r               = degree - multiplicity;

  std::vector<ScalarType> newKnots(r);
  std::ranges::fill(newKnots, u);
  auto tmpPatchData = Splines::knotRefinement(patchData, newKnots, 0);

  std::vector<ScalarType> tmpKnots{tmpPatchData.knotSpans.front()};
  const auto& tmpCp = tmpPatchData.controlPoints.directGetAll();

  std::vector<ScalarType> leftKnots;
  const int span_l = Splines::findSpan(degree, u, tmpKnots) + 1;
  // insert all old knots left of u
  leftKnots.insert(leftKnots.end(), tmpKnots.begin(), std::next(tmpKnots.begin(), span_l));
  // append u once since we already inserted it above r times, so we need to insert it r+1 times, which provides C^0
  // continuity there
  leftKnots.push_back(u);

  std::vector<ScalarType> rightKnots;
  // insert u degree+1 times at the front to provide C^0 continuity there
  rightKnots.insert(rightKnots.end(), degree + 1, u);
  // append the rest of the old knots as they are
  rightKnots.insert(rightKnots.end(), std::next(tmpKnots.begin(), span_l), tmpKnots.end());

  std::vector<ControlPointType> leftControlPoints;
  const int ks = span - degree + 1;
  // add all controlpoints that live left of u, which is the range [0,ks+r[ of controlpoints
  leftControlPoints.insert(leftControlPoints.end(), tmpCp.begin(), std::next(tmpCp.begin(), ks + r));

  std::vector<ControlPointType> rightControlPoints;
  // add all controlpoints that live right of u, which is the range [ks+r-1,end[ of controlpoints
  rightControlPoints.insert(rightControlPoints.end(), std::next(tmpCp.begin(), ks + r - 1), tmpCp.end());

  assert(std::ranges::is_sorted(leftKnots) && "The left knotvector should be sorted, but is not.");
  assert(std::ranges::is_sorted(rightKnots) && "The right knotvector should be sorted, but is not.");

  assert(std::ranges::count(leftKnots.begin(), leftKnots.begin() + degree + 1, leftKnots.front()) == degree + 1);
  assert(std::ranges::count(leftKnots.end() - degree - 1, leftKnots.end(), leftKnots.back()) == degree + 1);

  assert(std::ranges::count(rightKnots.begin(), rightKnots.begin() + degree + 1, rightKnots.front()) == degree + 1);
  assert(std::ranges::count(rightKnots.end() - degree - 1, rightKnots.end(), rightKnots.back()) == degree + 1);

  PatchData leftPatchData{{leftKnots}, ControlPointNetType{leftControlPoints}, patchData.degree};
  PatchData rightPatchData{{rightKnots}, ControlPointNetType{rightControlPoints}, patchData.degree};

  return std::make_pair(GeoCurve{leftPatchData}, GeoCurve{rightPatchData});
}

template <typename GeoCurve, typename ScalarType = typename GeoCurve::ctype>
requires(GeoCurve::mydimension == 1)
auto sliceCurve(const GeoCurve& geoCurve, std::array<ScalarType, 2> t) -> GeoCurve {
  auto domain = geoCurve.domain();
  if (FloatCmp::eq(domain[0].front(), t[0]) and FloatCmp::eq(domain[0].back(), t[1]))
    return geoCurve;
  if (FloatCmp::eq(domain[0].front(), t[0]))
    return std::get<0>(splitCurve(geoCurve, t[1]));
  if (FloatCmp::eq(domain[0].back(), t[1]))
    return std::get<1>(splitCurve(geoCurve, t[0]));

  auto tmpCurve = std::get<1>(splitCurve(geoCurve, t[0]));
  return std::get<0>(splitCurve(tmpCurve, t[1]));
}

} // namespace Dune::IGA

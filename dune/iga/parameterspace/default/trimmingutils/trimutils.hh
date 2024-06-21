// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <clipper2/clipper.core.h>
#include <ranges>

#include "dune/iga/geometrykernel/findintersection.hh"
#include "dune/iga/geometrykernel/nurbspatchgeometry.hh"
#include "dune/iga/parameterspace/default/trimmingutils/cliputils.hh"
#include <dune/common/float_cmp.hh>
#include <dune/iga/geometrykernel/geohelper.hh>

namespace Dune::IGA::DefaultParameterSpace::Util {
template <typename ScalarType, int dim>
auto approxSamePoint(const Clipper2Lib::PointD& pt1, const FieldVector<ScalarType, dim>& pt2, const double prec)
    -> bool {
  return FloatCmp::eq(pt1.x, pt2[0], prec) and FloatCmp::eq(pt1.y, pt2[1], prec);
}

auto distance(const Clipper2Lib::PointD& p1, const auto& p2) -> double {
  return std::hypot(p1.x - p2[0], p1.y - p2[1]);
}

auto findGoodStartingPoint(const auto& curve, const Clipper2Lib::PointD& pt) -> double {
  int N          = std::max(curve.numberOfControlPoints().front() * 30, 200);
  auto linSpace  = Utilities::linspace(curve.domain().front(), N);
  auto distances = std::ranges::transform_view(linSpace, [&](const auto u) { return distance(pt, curve.global(u)); });
  auto min_idx   = std::ranges::distance(distances.begin(), std::ranges::min_element(distances));
  return linSpace[min_idx];
}

template <typename TrimmingCurve>
auto createTrimmingCurveSlice(const TrimmingCurve& curve, double t1, double t2) -> TrimmingCurve {
  return sliceCurve(curve, {t1, t2});
}

template <typename TrimmingCurve>
auto createHostGeometry(auto& vertex1, auto& vertex2) -> TrimmingCurve {
  const std::array<std::vector<double>, 1> knotSpans = {{{0, 0, 1, 1}}};

  const std::vector<typename TrimmingCurve::ControlPointType> controlPoints = {
      {{.p = vertex1, .w = 1}, {.p = vertex2, .w = 1}}
  };
  auto controlNet = NURBSPatchData<1, 2>::ControlPointNetType(controlPoints);
  typename TrimmingCurve::PatchData patchData(knotSpans, controlNet, {1});

  return GeometryKernel::NURBSPatch(patchData);
}

enum class FindIntersectionError
{
  parallel,
  notSuccessful
};
struct FindIntersectionResult
{
  double tParam{};
  FieldVector<double, 2> curvePoint{};
};

auto callFindIntersection(const auto& curvePatchGeo, int edgeIdx, const auto& ip, const auto& corners)
    -> Std::expected<FindIntersectionResult, FindIntersectionError> {
  auto pos         = corners[edgeIdx];
  auto dir         = edgeDirections[edgeIdx];
  double lineGuess = (edgeIdx == 0 or edgeIdx == 2) ? (ip.x - pos[0]) / dir[0] : (ip.y - pos[1]) / dir[1];

  // Catch a trivial but difficult case (algorithm sometimes is not able to converge to endPoint)
  if (auto endPoint = curvePatchGeo.corner(1);
      approxSamePoint(ip, endPoint, 1e-8) and not FloatCmp::eq(endPoint, curvePatchGeo.corner(0))) {
    return FindIntersectionResult{curvePatchGeo.domain()[0][1], endPoint};
  }

  // Todo use z-Value to get a good starting point, not brute-force
  auto guessTParam = FieldVector<double, 2>{findGoodStartingPoint(curvePatchGeo, ip), lineGuess};

  auto [success, tParam, curvePoint] = findIntersectionCurveAndLine(curvePatchGeo, pos, dir, guessTParam);
  if (success == IntersectionCurveAndLine::intersect)
    return FindIntersectionResult(tParam[0], curvePoint);
  if (success == IntersectionCurveAndLine::parallel) {
    return Std::unexpected(FindIntersectionError::parallel);
  }

  auto domain = curvePatchGeo.domain();
  if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].front()})))
    return FindIntersectionResult(domain[0][0], curvePoint);
  if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].back()})))
    return FindIntersectionResult(domain[0].back(), curvePoint);

  return Std::unexpected(FindIntersectionError::notSuccessful);
};

template <typename Entity>
Entity coarsestFather(const Entity& ele) {
  assert(ele.hasFather());
  Entity father = ele.father();
  while (father.hasFather()) {
    father = father.father();
  }
  assert(father.level() == 0);
  return father;
}

} // namespace Dune::IGA::DefaultParameterSpace::Util

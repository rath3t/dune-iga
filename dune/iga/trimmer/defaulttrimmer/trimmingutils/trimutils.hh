// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include <clipper2/clipper.core.h>
#include <ranges>

#include "dune/iga/geometrykernel/findintersection.hh"
#include "dune/iga/geometrykernel/nurbspatchgeometry.hh"
#include "dune/iga/trimmer/defaulttrimmer/trimmingutils/cliputils.hh"
#include <dune/common/float_cmp.hh>
#include <dune/iga/geometrykernel/geohelper.hh>

namespace Dune::IGANEW::DefaultTrim::Util {
template <typename ScalarType, int dim>
auto approxSamePoint(const Clipper2Lib::PointD& pt1, const FieldVector<ScalarType, dim>& pt2,
                     const double prec) -> bool {
  return FloatCmp::eq(pt1.x, pt2[0], prec) and FloatCmp::eq(pt1.y, pt2[1], prec);
}

auto distance(const Clipper2Lib::PointD& p1, const auto& p2) -> double {
  return std::hypot(p1.x - p2[0], p1.y - p2[1]);
}

auto findGoodStartingPoint(const auto& curve, const Clipper2Lib::PointD& pt, int N = 200) -> double {
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

auto callFindIntersection(const auto& curvePatchGeo, int edgeIdx, const auto& ip,
                          const auto& corners) -> std::pair<double, FieldVector<double, 2>> {
  auto pos = corners[edgeIdx];
  auto dir = edgeDirections[edgeIdx];

  double lineGuess = (edgeIdx == 0 or edgeIdx == 2) ? (ip.x - pos[0]) / dir[0] : (ip.y - pos[1]) / dir[1];
  auto guessTParam = FieldVector<double, 2>{findGoodStartingPoint(curvePatchGeo, ip), lineGuess};

  auto [success, tParam, curvePoint] = findIntersectionCurveAndLine(curvePatchGeo, pos, dir, guessTParam);
  if (success == IntersectionCurveAndLine::intersect)
    return std::make_pair(tParam[0], curvePoint);
  if (success == IntersectionCurveAndLine::parallel) {
    DUNE_THROW(Dune::GridError, "Couldn't find intersection Point, lines are parallel");
  }

  auto domain = curvePatchGeo.domain();
  if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].front()})))
    return std::make_pair(domain[0][0], curvePoint);
  if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].back()})))
    return std::make_pair(domain[0].back(), curvePoint);

  DUNE_THROW(Dune::GridError, "Couldn't find intersection Point");
};

} // namespace Dune::IGANEW::DefaultTrim::Util

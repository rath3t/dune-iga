#pragma once

namespace Dune::IGANEW::DefaultTrim::Util {
  template <typename ScalarType, int dim>
  auto approxSamePoint(const Clipper2Lib::PointD& pt1, const FieldVector<ScalarType, dim>& pt2, const double prec)
      -> bool {
    return FloatCmp::eq(pt1.x, pt2[0], prec) and FloatCmp::eq(pt1.y, pt2[1], prec);
  }

  auto distance(const Clipper2Lib::PointD& p1, const auto& p2) -> double {
    return std::hypot(p1.x - p2[0], p1.y - p2[1]);
  }

  auto findGoodStartingPoint(const auto& curve, const Clipper2Lib::PointD& pt, int N = 100) -> double {
    auto linSpace = Utilities::linspace(curve.domain().front(), N);
    std::vector<double> distances;
    std::ranges::transform(linSpace, std::back_inserter(distances),
                           [&](const auto u) { return distance(pt, curve.global(u)); });
    auto min_idx = std::ranges::distance(distances.begin(), std::ranges::min_element(distances));
    return linSpace[min_idx];
  }

  template <typename TrimmingCurve>
  auto createTrimmingCurveSlice(const TrimmingCurve& curve, double t1, double t2) -> TrimmingCurve {
    return sliceCurve(curve, {t1, t2});
  }

  template <typename TrimmingCurve>
  auto createHostGeometry(auto& vertex1, auto& vertex2) -> TrimmingCurve {
    const std::array<std::vector<double>, 1> knotSpans = {{{0, 0, 1, 1}}};

    const std::vector<typename TrimmingCurve::ControlPointType> controlPoints
        = {{{.p = vertex1, .w = 1}, {.p = vertex2, .w = 1}}};
    auto controlNet = NURBSPatchData<1, 2>::ControlPointNetType(controlPoints);
    typename TrimmingCurve::PatchData patchData(knotSpans, controlNet, {1});

    return GeometryKernel::NURBSPatch(patchData);
  }

  auto callFindIntersection(const auto& curvePatchGeo, const int edgeIdx, const auto& ip, const auto& corners)
      -> std::pair<double, FieldVector<double, 2>> {
    // @todo gerneralize for more than one loop
    auto pos = corners[edgeIdx];
    auto dir = edgeDirections[edgeIdx];

    auto guessTParam = FieldVector<double, 2>{findGoodStartingPoint(curvePatchGeo, ip), 0.5};

    auto [success, tParam, curvePoint] = findIntersectionCurveAndLine(curvePatchGeo, pos, dir, guessTParam);
    if (success == FindIntersectionCurveAndLineResult::sucess) return std::make_pair(tParam[0], curvePoint);
    if (success == FindIntersectionCurveAndLineResult::linesParallel) {
      DUNE_THROW(Dune::GridError, "Couldn't find intersection Point, lines are parallel");
    }

    assert(Util::approxSamePoint(ip, curvePoint, 1e-4) && "Intersection Point and Curve Point are not the same");

    auto domain = curvePatchGeo.domain();
    if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].front()})))
      return std::make_pair(domain[0][0], curvePoint);
    if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].back()})))
      return std::make_pair(domain[0].back(), curvePoint);

    DUNE_THROW(Dune::GridError, "Couldn't find intersection Point");
  };
}  // namespace Dune::IGANEW::DefaultTrim::Util
//
// Created by Henri on 14.03.2023.
//

#pragma once

namespace Dune::IGA::Impl::Trim {

  using ctype        = double;
  using ClipperPoint = Clipper2Lib::Point<ctype>;
  using Point        = FieldVector<ctype, 2>;
  using PointVector  = std::vector<Point>;
  
  using CurveGeometry = NURBSPatchGeometry<1, 2>;

  // To keep track of which point belongs to which edge or node we will use a pointMap
  struct IntersectionPoint {
    Point point;
    bool alreadyVisited = false;

    explicit IntersectionPoint(ClipperPoint& _clipperPoint) : point({_clipperPoint.x, _clipperPoint.y}){};
  };

  using IntersectionPointMap    = std::map<int, std::vector<IntersectionPoint>>;
  using IntersectionPointMapPtr = std::shared_ptr<std::map<int, std::vector<IntersectionPoint>>>;

  PointVector getElementCorners(const auto& element) {
    auto geo = element.geometry();

    std::vector<Point> corners(4);
    corners[0] = geo.corner(0);
    corners[1] = geo.corner(1);
    corners[2] = geo.corner(3);  // see dune book page 127 Figure 5.12
    corners[3] = geo.corner(2);

    return corners;
  }

  Clipper2Lib::PathD getElementEdges(const auto& element) {
    auto corners = getElementCorners(element);
    Clipper2Lib::PathD edges;
    for (int i = 0; i < 4; ++i)
      edges.emplace_back(corners[i][0], corners[i][1]);
    return edges;
  }

  bool areEdgesParallel(CurveGeometry& curve1, CurveGeometry& curve2) {
    if (curve1.degree()[0] > 1 || curve2.degree()[0] > 1)
      return false;

    const double tolerance = double(16) * std::numeric_limits<double>::epsilon();

    auto domain1 = curve1.domain()[0];
    auto domain2 = curve2.domain()[0];

    auto p1 = curve1(domain1[0]);
    auto p2 = curve1(domain1[1]);

    auto p3 = curve2(domain2[0]);
    auto p4 = curve2(domain2[1]);

    // Compute direction vectors of both curves
    auto d1 = p2 - p1;
    auto d2 = p4 - p3;

    auto cos_angle = Dune::dot(d1, d2) / (d1.two_norm() * d2.two_norm());
    auto angle     = std::acos(std::min(std::max(cos_angle, -1.0), 1.0));

    // Check if the angle is small enough to consider the curves on the same line
    return (Dune::FloatCmp::eq(angle, 0.0, tolerance) ||
            Dune::FloatCmp::eq(angle, std::numbers::pi, tolerance));

  }

  struct IntersectionResult {
    double localResult{};
    Point globalResult;
    bool found = false;
  };

  IntersectionResult newtonRaphson(CurveGeometry& curve, Point& intersectionPoint) {
    const double pointEqualityTolerance = 1e-7;
    const double objectiveFctTolerance = 1e-4;

    // There are some easy cases where the inital Guess is at the back or at the front of the domain, check these cases
    // first
    auto domain = curve.domain()[0];

    if ((curve(domain.front()) - intersectionPoint).two_norm() < pointEqualityTolerance)
      return {domain.front(), curve(domain.front()), true};
    if ((curve(domain.back()) - intersectionPoint).two_norm() < pointEqualityTolerance)
      return {domain.back(), curve(domain.back()), true};

    // At first, we assume that the intersectionPoint is in some definition exact and the corresponding u has to be
    // found
    Dune::FieldVector<double, 1> u{curve.domainMidPoint()[0]};
    Dune::FieldVector<double, 1> du;

    // Set up the algorithm
    Point dGlobal;
    int max_iterations = 250;
    int i              = 0;
    do {
      dGlobal = curve(u[0]) - intersectionPoint;

      using MatrixHelper = typename Dune::MultiLinearGeometryTraits<double>::MatrixHelper;
      MatrixHelper::template xTRightInvA<1, 2>(curve.jacobianTransposed(u[0]), dGlobal, du);

      // Update the guess
      u -= du;

      ++i;
      if (i > max_iterations) return {u[0], intersectionPoint, false};

      // Clamp result, this might be already an indicator that there is no solution
      if (Dune::FloatCmp::gt(u[0], curve.patchData_->knotSpans[0].back())) u[0] = curve.patchData_->knotSpans[0].back();
      if (Dune::FloatCmp::lt(u[0], curve.patchData_->knotSpans[0].front())) u[0] = curve.patchData_->knotSpans[0].front();

    } while (dGlobal.two_norm() > objectiveFctTolerance);

    double localResult = u[0];
    auto globalResult  = curve(localResult);

    return {localResult, globalResult, true};
  }

  bool approxSamePoint(Point a, Point b, double tol = 1e-2) {
    return Dune::FloatCmp::eq(a, b, tol);
  }

  int getEdgeOrientation(int edge) {
    return  (edge == 0 || edge == 2) ? 0 : 1;
  }

  int iPlusX(int i, int x){
    return (i + x) % 4;
  }

  Boundary boundaryForEdge(PointVector corners, int edge) {
    PointVector v{corners[iPlusX(edge, 0)], corners[iPlusX(edge, 1)]};
    return Boundary{v};
  }

  void sortIntersectionPoints(const IntersectionPointMapPtr& pointMapPtr, int edge) {
    auto orientation = getEdgeOrientation(edge);
    auto cmp = [orientation, edge](const IntersectionPoint& a, const IntersectionPoint& b) {
      return (edge == 0 || edge == 1) ? a.point[orientation] < b.point[orientation]
                                      : a.point[orientation] > b.point[orientation];
    };
    std::ranges::sort((*pointMapPtr)[edge], cmp);
  }

  std::pair<Point, bool> getNextIntersectionPointOnEdge(const IntersectionPointMapPtr& pointMapPtr, int edge) {
    if ((*pointMapPtr)[edge].empty()) return {Point(), false};

    for (const IntersectionPoint& ip : (*pointMapPtr)[edge]) {
      if (!ip.alreadyVisited) return {ip.point, true};
    }
    return {Point(), false};
  }

  void markIntersectionPointAsVisited(const IntersectionPointMapPtr& pointMapPtr, int edge, Point& point) {
    auto it = std::find_if((*pointMapPtr)[edge].begin(), (*pointMapPtr)[edge].end(), [&point](IntersectionPoint x) {
      return (x.point[0] == point[0] && x.point[1] == point[1]);
    });
    assert(it != (*pointMapPtr)[edge].end());
    (*it).alreadyVisited = true;
  }

  std::tuple<int, IntersectionResult, bool> findBoundaryThatHasIntersectionPoint(
      const IntersectionPointMapPtr& pointMapPtr, int edge, PointVector& corners,
      std::vector<Boundary>& globalBoundaries) {
    auto edgeCurve = boundaryForEdge(corners, edge);

    Point intersectionPointGuess = getNextIntersectionPointOnEdge(pointMapPtr, edge).first;

    // Get the exact edge Intersection Point
    int boundaryIndexThatHasIntersectionPoint = -1;
    IntersectionResult edgeIntersectResult;
    for (int counter = 0; auto boundary : globalBoundaries) {

      if (areEdgesParallel(boundary.nurbsGeometry, edgeCurve.nurbsGeometry)) {
        counter++;
        continue;
      }

      edgeIntersectResult = newtonRaphson(boundary.nurbsGeometry, intersectionPointGuess);

      if (edgeIntersectResult.found) {
        boundaryIndexThatHasIntersectionPoint = counter;
        break;
      }
      counter++;
    }

    if (boundaryIndexThatHasIntersectionPoint >= 0) {
      // Mark this intersectionPoint as visited
      markIntersectionPointAsVisited(pointMapPtr, edge, intersectionPointGuess);

      return {boundaryIndexThatHasIntersectionPoint, edgeIntersectResult, true};
    } else {
      return {-1, IntersectionResult(), false};
    }
  }

  std::tuple<double, Point, int> traceCurve(const IntersectionPointMapPtr& pointMapPtr, PointVector& corners,
                                            CurveGeometry& curveToTrace, const double startU) {
    // Todo: evtl. hört bei einem Loch eine Kurven Parametrisierung auf und es beginnt eine neue, die am gleichen Punkt
    // weitergeht

    // We need to find the next guess for an intersection Point on this curve, then find the intersection
    // Point exactly If there is no exact Intersection Point, then there is no intersection
    bool otherDirectionChecked   = false;
    bool checkedConsecutiveCurve = false;

    bool found    = false;
    auto currentI = 0;
    Point foundInterSectionPoint;
    int foundOnEdgeNr = -1;
    double intersectU;

    // Gehe durch die Punkte der Kurve durch, schaue nach möglichen Kandidaten für IntertsectionPoints, kann relativ
    // ungenau sein das onAnyEdge muss weg
    auto maxU   = curveToTrace.domain()[0].back();
    auto beginU = startU;

  startOver:
    auto linSpace = Utilities::linspace(beginU, maxU, 350);
    for (int k = currentI; k < linSpace.size(); ++k) {
      // Evaluate curve
      auto currentPointOnCurveToTrace = curveToTrace(linSpace[k]);

      for (int current_edge = 0; current_edge < 4; ++current_edge) {
        std::pair<Point, bool> nextIntersectionPointOnEdgeResult
            = getNextIntersectionPointOnEdge(pointMapPtr, current_edge);
        if (!nextIntersectionPointOnEdgeResult.second) continue;

        // Teste Punkt currentPointOnCurveToTrace gegen den Möglichen IntersectionPoint
        if (approxSamePoint(currentPointOnCurveToTrace, nextIntersectionPointOnEdgeResult.first)) {
          // Now we can calculate the exact U for the intersectionPoint

          Boundary edgeBoundary = boundaryForEdge(corners, current_edge);

          int dir                 = getEdgeOrientation(current_edge);
          IntersectionResult res2 = newtonRaphson(curveToTrace, nextIntersectionPointOnEdgeResult.first);

          if (res2.found) {
            foundInterSectionPoint = res2.globalResult;
            intersectU             = res2.localResult;
            foundOnEdgeNr          = current_edge;

            // Mark as visited
            markIntersectionPointAsVisited(pointMapPtr, current_edge, nextIntersectionPointOnEdgeResult.first);

            break;
          }
        }
      }

      if (foundOnEdgeNr >= 0) {
        found = true;
        break;
      }
    }
    if (!found && !otherDirectionChecked) {
      beginU                = curveToTrace.domain()[0].front();
      maxU                  = startU;
      otherDirectionChecked = true;
      goto startOver;
    }
    if (!found) throw std::runtime_error("No Intersection Point found while tracing curve.");


    return {intersectU, foundInterSectionPoint, foundOnEdgeNr};
  }

  bool pointOnAnyEdge(ClipperPoint point, auto path) {
    auto distance = [](ClipperPoint p1, ClipperPoint p2) {
      return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
    };

    auto pointOnLine = [distance](ClipperPoint p, ClipperPoint a, ClipperPoint b) {
      // https://stackoverflow.com/a/17693146
      const double tolerance = 1e-8;
      return Dune::FloatCmp::eq(distance(a, p) + distance(b, p), distance(a, b), tolerance);
    };

    for (int i = 0; i < 4; ++i)
      if (pointOnLine(point, path[i], path[iPlusX(i, 1)])) return true;

    return false;
  }

}  // namespace Dune::IGA::Impl::Trim
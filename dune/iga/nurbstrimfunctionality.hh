//
// Created by Henri on 14.03.2023.
//

#pragma once

#include <algorithm>
#include <bits/ranges_algo.h>
#include <bits/ranges_base.h>
#include <clipper2/clipper.core.h>
#include <clipper2/clipper.h>
#include <iterator>
#include <memory>
#include <optional>

#include "dune/common/float_cmp.hh"
#include "dune/common/fvector.hh"
#include "dune/iga/nurbstrimutils.hh"
#include <dune/iga/nurbspatchgeometry.h>
#include <dune/iga/nurbstrimboundary.hh>

namespace Dune::IGA {
  enum class ElementTrimFlag { full, empty, trimmed };
}

namespace Dune::IGA::Impl::Trim {

  // The following numbering is used for edges and nodes:
  /*
   * 3 -- 2 -- 2
   * |         |
   * 3         1
   * |         |
   * 0 -- 0 -- 1
   */

  // An int64 is able to hold numbers up to 9e18, so scaling by 1e12 should be fine
  using ctype = int64_t;

  constexpr int scale      = 12;
  const double scaleFactor = std::pow(10, scale);

  constexpr int pathSamples  = 800;
  constexpr double tolerance = 1e-8;

  /* Newton Raphson Control Variables */
  constexpr double objectiveFctTolerance       = 1e-4;
  constexpr int probesForStartingPoint         = 5;
  constexpr int fallbackProbesForStartingPoint = 200;
  constexpr int maxIterations                  = 50;

  // constexpr double epsilonPrecision = double(16) * std::numeric_limits<double>::epsilon();

  using ClipperPoint = Clipper2Lib::Point<ctype>;
  using ClipperPath  = Clipper2Lib::Path<ctype>;
  using ClipperPaths = Clipper2Lib::Paths<ctype>;

  using Point    = FieldVector<double, 2>;
  using IntPoint = FieldVector<ctype, 2>;

  using CurveGeometry = NURBSPatchGeometry<1, 2>;

  ctype toIntDomain(double x) { return x * scaleFactor; }
  FieldVector<ctype, 2> toIntDomain(FieldVector<double, 2> x) { return x * scaleFactor; }

  double toFloatDomain(ctype x) { return x / scaleFactor; }
  FieldVector<double, 2> toFloatDomain(FieldVector<ctype, 2> x) { return x / scaleFactor; }

  // To keep track of which point belongs to which edge or node we will use a pointMap IntersectionPoints
  // IntersectionPoint lives in the FloatDomain
  struct IntersectionPoint {
    Point point;
    bool alreadyVisited = false;

    explicit IntersectionPoint(const ClipperPoint& _clipperPoint)
        : point({toFloatDomain(_clipperPoint.x), toFloatDomain(_clipperPoint.y)}){};
  };

  using IntersectionPointMap = std::array<std::vector<IntersectionPoint>, 4>;

  // Helper functions
  double distance(const ClipperPoint p1, const ClipperPoint p2) { return std::hypot(p1.x - p2.x, p1.y - p2.y); }

  double distance(const Point p1, const Point p2) { return std::hypot(p1[0] - p2[0], p1[1] - p2[1]); }

  bool samePoint(const IntPoint& a, const ClipperPoint& b) { return a == IntPoint{b.x, b.y}; }

  inline int getEdgeOrientation(const int edge) { return (edge == 0 || edge == 2) ? 0 : 1; }

  inline int nextEntityIdx(const int i, const int x) { return (i + x) % 4; }

  std::array<Point, 2> curveStartEndPoint(const CurveGeometry& curve) {
    std::array<int, 2> indices{0, static_cast<int>(curve.patchData_->controlPoints.size()[0]) - 1};
    return {curve.patchData_->controlPoints.get({indices[0]}).p, curve.patchData_->controlPoints.get({indices[1]}).p};
  }

  template <typename T = ctype, int dim = 2>
  std::vector<FieldVector<T, dim>> getElementCorners(const auto& element) {
    auto geo = element.geometry();

    // The corners are sorted counter-clockwise to get a closed loop, also see top of this file for numbering
    return {geo.corner(0), geo.corner(1), geo.corner(3), geo.corner(2)};
  }

  Clipper2Lib::Path<ctype> getElementEdgesFromElementCorners(const std::vector<IntPoint>& corners) {
    Clipper2Lib::Path<ctype> edges;
    for (int i = 0; i < 4; ++i)
      edges.emplace_back(corners[i][0], corners[i][1]);
    return edges;
  }

  Clipper2Lib::Path<ctype> getElementEdges(const auto& element) {
    auto corners = getElementCorners(element);
    return getElementEdgesFromElementCorners(corners);
  }

  Clipper2Lib::Paths<ctype> getClip(TrimData* trimData) {
    Clipper2Lib::Paths<ctype> clipPaths;
    Clipper2Lib::Path<ctype> tempPath;
    for (auto& loop : trimData->boundaryLoops) {
      tempPath.clear();
      for (auto& boundary : loop.boundaries)
        for (auto& point : boundary.path<ctype>(pathSamples, scale))
          tempPath.push_back(point);

      // Sanitize Path, so we don't get duplicate intersection Points on corners
      tempPath = Clipper2Lib::TrimCollinear(tempPath);
      clipPaths.push_back(tempPath);
    }

    return clipPaths;
  }

  /// \brief Checks if a point p is on the line spanned by the points a and b within a given tolerance
  bool pointOnLine(const ClipperPoint& p, const ClipperPoint& a, const ClipperPoint& b, const double tol = tolerance) {
    // https://stackoverflow.com/a/17693146
    return distance(a, p) + distance(b, p) == distance(a, b);
  }

  std::pair<bool, int> pointOnAnyEdge(const ClipperPoint& point, const auto& elementEdges) {
    for (int i = 0; i < 4; ++i)
      if (pointOnLine(point, elementEdges[i], elementEdges[nextEntityIdx(i, 1)])) return std::make_pair(true, i);

    return std::make_pair(false, -1);
  }

  std::pair<bool, int> pointOnAnyVertex(const ClipperPoint& point, const std::vector<IntPoint>& corners) {
    for (size_t i = 0; i < corners.size(); ++i)
      if (samePoint(corners[i], point)) return std::make_pair(true, static_cast<int>(i));

    return std::make_pair(false, -1);
  }

  // Sorts the Intersection Points counter-clockwise around the element (crucial step!)
  void sortIntersectionPoints(IntersectionPointMap& pointMap) {
    for (int edge = 0; edge < 2; ++edge)
      std::ranges::sort(pointMap[edge], [&](const IntersectionPoint& a, const IntersectionPoint& b) {
        auto orientation = getEdgeOrientation(edge);
        return a.point[orientation] < b.point[orientation];
      });
    for (int edge = 2; edge < 4; ++edge)
      std::ranges::sort(pointMap[edge], [&](const IntersectionPoint& a, const IntersectionPoint& b) {
        auto orientation = getEdgeOrientation(edge);
        return a.point[orientation] > b.point[orientation];
      });
  }

  bool pointInElement(const Point& point, std::vector<IntPoint>& corners) {
    auto lowerLeft  = toFloatDomain(corners[0]);
    auto upperRight = toFloatDomain(corners[2]);

    return FloatCmp::gt(point[0], lowerLeft[0]) and FloatCmp::lt(point[0], upperRight[0])
           and FloatCmp::gt(point[1], lowerLeft[1]) and FloatCmp::lt(point[1], upperRight[1]);
  }
  bool pointInElementOrOnEdge(const Point& point, std::vector<IntPoint>& corners) {
    auto lowerLeft  = toFloatDomain(corners[0]);
    auto upperRight = toFloatDomain(corners[2]);

    return FloatCmp::ge(point[0], lowerLeft[0], tolerance) and FloatCmp::le(point[0], upperRight[0], tolerance)
           and FloatCmp::ge(point[1], lowerLeft[1], tolerance) and FloatCmp::le(point[1], upperRight[1], tolerance);
  }

  bool isFullElement(ClipperPath& clippedEdges, std::vector<IntPoint>& corners) {
    if (clippedEdges.size() != 4) return false;

    return std::ranges::all_of(clippedEdges, [&corners](const ClipperPoint& point) {
      auto it = std::ranges::find_if(corners, [&point](const IntPoint& corner) { return samePoint(corner, point); });
      return it != corners.end();
    });
  }

  struct ClippingResult {
    std::unique_ptr<IntersectionPointMap> pointMapPtr;
    std::array<int, 4> edgeCounter;
    std::array<int, 4> nodeCounter;
  };

  std::pair<ElementTrimFlag, std::optional<ClippingResult>> clipElement(const auto& element, ClipperPaths& clip) {
    // Prepare Result
    ElementTrimFlag trimFlag;

    auto corners      = getElementCorners(element);
    auto elementEdges = getElementEdges(element);

    ClipperPaths clippedEdges;
    if (clip.size() == 1) {
      Clipper2Lib::Rect<ctype> elementRect{corners[0][0], corners[1][1], corners[1][0], corners[3][1]};
      clippedEdges = Clipper2Lib::RectClip(elementRect, clip);
    } else
      clippedEdges
          = Clipper2Lib::Intersect(ClipperPaths{getElementEdges(element)}, clip, Clipper2Lib::FillRule::EvenOdd);

    // At the moment there is no hole, if there are more than 2 Paths, then there is a hole
    if (clippedEdges.size() > 1) {
      std::cerr << "Hole detected in element, hole gets ignored" << std::endl;
      // clippedEdges.erase(clippedEdges.begin() + 1, clippedEdges.end());
      return std::make_pair(ElementTrimFlag::full, std::nullopt);
    }

    // If the clippedEdges are empty, this means the element is outside the clip -> empty
    // And if the clip has only 4 corners and these corners are the same as the element -> full
    if (clippedEdges.empty())
      trimFlag = ElementTrimFlag::empty;
    else {
      if (isFullElement(clippedEdges.front(), corners)) {
        trimFlag = ElementTrimFlag::full;
      } else
        trimFlag = ElementTrimFlag::trimmed;
    }

    if (trimFlag != ElementTrimFlag::trimmed) return std::make_pair(trimFlag, std::nullopt);

    // If the code is here → Elements are trimmed

    // Fill pointMap to store the intersection points in regard to their corresponding edges
    auto pointMap = std::make_unique<IntersectionPointMap>();

    // This keeps track of how many intersection Points there are at a given edge
    std::array<int, 4> edgeCounter{0, 0, 0, 0};
    std::array<int, 4> nodeCounter{0, 0, 0, 0};

    for (auto& point : clippedEdges.front()) {
      auto [isOnAnyEdge, edgeThePointIsOn] = pointOnAnyEdge(point, elementEdges);
      auto [isOnAnyNode, nodeThePointIsOn] = pointOnAnyVertex(point, corners);

      // If is onAnyNode is true, then isOnAnyEdge automatically is also true, so test onAnyNode first
      if (isOnAnyNode)
        nodeCounter[nodeThePointIsOn] = 1;
      else if (isOnAnyEdge) {
        (*pointMap)[edgeThePointIsOn].emplace_back(point);
        ++edgeCounter[edgeThePointIsOn];
      }
    }

    assert(not std::ranges::all_of(nodeCounter, [](auto x) { return x == 1; }));

    // Sort IntersectionPoints counter-clockwise
    sortIntersectionPoints(*pointMap);

    return std::make_pair(trimFlag,
                          std::make_optional<ClippingResult>({std::move(pointMap), edgeCounter, nodeCounter}));
  }

  bool areLinesParallel(const CurveGeometry& line1, const CurveGeometry& line2) {
    if (line1.degree()[0] > 1 || line2.degree()[0] > 1) return false;

    auto [p1, p2] = curveStartEndPoint(line1);
    auto [p3, p4] = curveStartEndPoint(line2);

    // Compute direction vectors of both curves
    auto d1 = p2 - p1;
    auto d2 = p4 - p3;

    auto cos_angle = static_cast<double>(Dune::dot(d1, d2) / (d1.two_norm() * d2.two_norm()));
    auto angle     = std::acos(std::min(std::max(cos_angle, -1.0), 1.0));

    // Check if the angle is small enough to consider the curves on the same line (or parallel)
    return (Dune::FloatCmp::eq(angle, 0.0, tolerance) || Dune::FloatCmp::eq(angle, std::numbers::pi, tolerance));
  }

  struct IntersectionResult {
    double localResult{};
    Point globalResult;
    bool found = false;
  };

  double findGoodStartingPoint(CurveGeometry& curve, Point& intersectionPoint, int N) {
    auto linSpace = Utilities::linspace(curve.domain()[0], N);
    std::vector<double> distances;
    std::ranges::transform(linSpace, std::back_inserter(distances),
                           [&](const auto& u) { return distance(curve(u), intersectionPoint); });
    auto min_it  = std::ranges::min_element(distances);
    auto min_idx = std::ranges::distance(distances.begin(), min_it);
    return linSpace[min_idx];
  }

  IntersectionResult newtonRaphson(CurveGeometry& curve, Point& intersectionPoint, int probes) {
    // There are some easy cases where the IP is at the back or at the front of the domain, check these cases first
    auto [start, end] = curveStartEndPoint(curve);

    auto [lowU, upperU] = curve.domain()[0];

    if (Dune::FloatCmp::eq(start, intersectionPoint, tolerance)) return {lowU, start, true};
    if (Dune::FloatCmp::eq(end, intersectionPoint, tolerance)) return {upperU, end, true};

    // We assume that the intersectionPoint is in some definition exact and the corresponding u has to be found
    Dune::FieldVector<double, 1> u{findGoodStartingPoint(curve, intersectionPoint, probes)};
    Dune::FieldVector<double, 1> du;

    // Set up the algorithm
    Point dGlobal;
    int i = 0;
    do {
      dGlobal = curve(u[0]) - intersectionPoint;

      using MatrixHelper = typename Dune::MultiLinearGeometryTraits<double>::MatrixHelper;
      MatrixHelper::template xTRightInvA<1, 2>(curve.jacobianTransposed(u[0]), dGlobal, du);

      // Update the guess
      u -= du;

      if (++i > maxIterations) return {u[0], intersectionPoint, false};

      // Clamp result, this might be already an indicator that there is no solution
      if (Dune::FloatCmp::gt(u[0], curve.patchData_->knotSpans[0].back())) u[0] = upperU;
      if (Dune::FloatCmp::lt(u[0], curve.patchData_->knotSpans[0].front())) u[0] = lowU;

    } while (dGlobal.two_norm() > objectiveFctTolerance);

    double localResult = u[0];
    auto globalResult  = curve(localResult);

    return {localResult, globalResult, true};
  }

  /// This scales back to the floatDomain
  Boundary boundaryForEdge(const std::vector<IntPoint>& corners, int edge) {
    return Boundary(toFloatDomain(corners[nextEntityIdx(edge, 0)]), toFloatDomain(corners[nextEntityIdx(edge, 1)]));
  }

  std::pair<Point, bool> getNextIntersectionPointOnEdge(IntersectionPointMap* pointMapPtr, int edge) {
    if ((*pointMapPtr)[edge].empty()) return {Point(), false};

    for (const IntersectionPoint& ip : (*pointMapPtr)[edge]) {
      if (!ip.alreadyVisited) return {ip.point, true};
    }
    return {Point(), false};
  }

  void markIntersectionPointAsVisited(IntersectionPointMap* pointMapPtr, int edge, const Point& point) {
    auto it = std::find_if((*pointMapPtr)[edge].begin(), (*pointMapPtr)[edge].end(),
                           [&point](IntersectionPoint x) { return (Dune::FloatCmp::eq(x.point, point)); });
    assert(it != (*pointMapPtr)[edge].end());
    (*it).alreadyVisited = true;
  }

  struct FindNextBoundaryResult {
    int edge{};
    IntersectionResult result;
  };


  std::optional<FindNextBoundaryResult> findBoundaryThatHasIntersectionPoint(
      IntersectionPointMap* pointMapPtr, int edge, std::vector<IntPoint>& corners,
      std::vector<Boundary>& globalBoundaries) {

    auto edgeCurve               = boundaryForEdge(corners, edge);
    Point intersectionPointGuess = getNextIntersectionPointOnEdge(pointMapPtr, edge).first;

    for (auto probes : {probesForStartingPoint, fallbackProbesForStartingPoint}) {
      for (int i = 0; auto& boundary : globalBoundaries) {
        if (areLinesParallel(boundary.nurbsGeometry, edgeCurve.nurbsGeometry)) {
          ++i;
          continue;
        }

        auto edgeIntersectResult = newtonRaphson(boundary.nurbsGeometry, intersectionPointGuess, probes);

        if (edgeIntersectResult.found) {
          markIntersectionPointAsVisited(pointMapPtr, edge, intersectionPointGuess);
          assert(Dune::FloatCmp::eq(edgeIntersectResult.globalResult, intersectionPointGuess, 1e-3));
          return std::make_optional(FindNextBoundaryResult(i, edgeIntersectResult));
        }
        ++i;
      }
    }

    return std::nullopt;
  }

  struct TraceCurveResult {
    double intersectU;
    Point intersectionPoint;
    int foundOnEdge;
  };

  std::optional<TraceCurveResult> traceCurveImpl(IntersectionPointMap* pointMapPtr, CurveGeometry& curveToTrace,
                                                 const int startEdge) {
    // Look through edges and determine the next Intersection Point (the IPs are sorted counter-clockwise)

    for (auto probes : {probesForStartingPoint, fallbackProbesForStartingPoint}) {
      for (int i = 1; i < 5; ++i) {
        auto currentEdge                               = nextEntityIdx(startEdge, i);
        auto [intersectionPoint, hasIntersectionPoint] = getNextIntersectionPointOnEdge(pointMapPtr, currentEdge);
        if (hasIntersectionPoint) {
          // Check if the curveToTrace goes through the intersection point
          auto [localResult, globalResult, found] = newtonRaphson(curveToTrace, intersectionPoint, probes);
          if (found) {
            markIntersectionPointAsVisited(pointMapPtr, currentEdge, intersectionPoint);
            return std::make_optional<TraceCurveResult>({localResult, globalResult, currentEdge});
          }
        }
      }
    }
    return std::nullopt;
  }

  std::vector<Boundary> extractBoundaries(const std::shared_ptr<TrimData>& data) {
    std::vector<Boundary> boundaries;
    for (auto& loop : data->boundaryLoops) {
      std::ranges::copy(loop.boundaries, std::back_inserter(boundaries));
    }
    return boundaries;
  }

  struct FindNextBoundaryLoopState {
    // Variables that need to be set
    const int nodeToBegin;
    std::vector<IntPoint> corners;
    std::shared_ptr<TrimData> trimData;
    std::vector<Boundary> boundaries;
    std::shared_ptr<IntersectionPointMap> pointMapPtr;
    const std::array<int, 4> edgeCounter;
    const std::array<int, 4> nodeCounter;

    // Loop States
    bool isCurrentlyOnNode            = true;
    bool moreIntersectionPointsOnEdge = false;
    int edgeTheLoopIsOn               = -1;  // This variable is only filed when !loopIsCurrentlyOnNode

    // Current positions (are initially set through the constructor)
    int currentNode;
    Point currentNodePoint;

    // Helper Functions
    bool hasIntersectionPointOnEdgeNr(int nr) { return edgeCounter[nr] > 0; };
    bool hasIntersectionPointOnNodeNr(int nr) { return nodeCounter[nr] > 0; };

    FindNextBoundaryLoopState(int _nodeToBegin, std::vector<IntPoint>& _corners, ClippingResult& _result,
                              const std::shared_ptr<TrimData>& _trimData)
        : nodeToBegin(_nodeToBegin),
          corners(_corners),
          pointMapPtr(std::move(_result.pointMapPtr)),
          edgeCounter(_result.edgeCounter),
          nodeCounter(_result.nodeCounter),
          trimData(_trimData),
          boundaries(extractBoundaries(_trimData)) {
      currentNode      = nodeToBegin;
      currentNodePoint = corners[currentNode];

    };

    FindNextBoundaryLoopState() = delete;
    FindNextBoundaryLoopState(const FindNextBoundaryLoopState& ) = delete;
  };

  struct TraceCurveInput {
    int boundaryIdx;
    double startU;
    int startEdge;
  };
  struct TraceCurveOutput {
    int boundaryIdx;
    double beginU;
    double intersectU;
    Point intersectionPoint;
    int foundOnEdge;
  };

  std::optional<std::vector<TraceCurveOutput>> traceCurve(FindNextBoundaryLoopState* state,
                                                          TraceCurveInput traceCurveInput) {
    CurveGeometry curveToTrace  = state->boundaries[traceCurveInput.boundaryIdx].nurbsGeometry;
    auto [startPoint, endPoint] = curveStartEndPoint(curveToTrace);

    // There is an edge case where we have a closed curve (e.g. a circle) that has endPoint and startPoint at the same
    // point and that point is in the element
    // TODO: Get rid of initializer lists
    if (Dune::FloatCmp::eq(startPoint, endPoint) && pointInElement(startPoint, state->corners)) {
      auto traceCurveResult = traceCurveImpl(state->pointMapPtr.get(), curveToTrace, traceCurveInput.startEdge);
      if (traceCurveResult.has_value())
        return std::make_optional<std::vector<TraceCurveOutput>>(
            {{traceCurveInput.boundaryIdx, traceCurveInput.startU, curveToTrace.domain()[0].back(),
              traceCurveResult->intersectionPoint, -1},
             {traceCurveInput.boundaryIdx, curveToTrace.domain()[0].front(), traceCurveResult->intersectU,
              traceCurveResult->intersectionPoint, traceCurveResult->foundOnEdge}});
      else
        return std::nullopt;
    }
    // ... or the point is on an edge, even nastier
    if (Dune::FloatCmp::eq(startPoint, endPoint) && pointInElementOrOnEdge(startPoint, state->corners)) {
      if (Dune::FloatCmp::ne(curveToTrace(traceCurveInput.startU), startPoint)) {
        auto traceCurveResult = traceCurveImpl(state->pointMapPtr.get(), curveToTrace, traceCurveInput.startEdge);
        if (traceCurveResult.has_value())
          return std::make_optional<std::vector<TraceCurveOutput>>(
              {{traceCurveInput.boundaryIdx, traceCurveInput.startU, curveToTrace.domain()[0].back(),
                traceCurveResult->intersectionPoint, traceCurveResult->foundOnEdge}});
      }
    }

    // Check the boundary that was given in traceInput
    auto traceCurveResult = traceCurveImpl(state->pointMapPtr.get(), curveToTrace, traceCurveInput.startEdge);

    if (traceCurveResult.has_value())
      return std::make_optional<std::vector<TraceCurveOutput>>(
          {{traceCurveInput.boundaryIdx, traceCurveInput.startU, traceCurveResult->intersectU,
            traceCurveResult->intersectionPoint, traceCurveResult->foundOnEdge}});

    // Now check if the startPoint or endPoint of the curve is in the element
    int status = -1;  // -1 means not found, 0 means startPoint, 1 means endPoint

    if (pointInElement(endPoint, state->corners))
      status = 1;
    else if (pointInElement(startPoint, state->corners))
      status = 0;

    // If neither start nor endpoint is in the element then traceCurve most likely made a mistake
    if (status == -1) return std::nullopt;

    // TODO take into account orientation of curve -> Simplification possible test with kragarm_trim_circle3 with 0 refinement

    // Now we check if there is another boundary that begins at this point and trace this potential curve
    Point checkPoint      = (status == 1) ? endPoint : startPoint;
    double endUOfCurve1   = (status == 1) ? curveToTrace.domain()[0].back() : curveToTrace.domain()[0].front();
    double startUOfCurve2 = std::numeric_limits<double>::quiet_NaN();

    // Find the next BoundarySegment
    auto it = std::ranges::find_if(state->boundaries, [&](auto& _boundary) {
      auto [startPointOfBoundaryToCheck, endPointOfBoundaryToCheck] = curveStartEndPoint(_boundary.nurbsGeometry);
      // Check if it's the same boundary Geometry
      if (startPointOfBoundaryToCheck == startPoint && endPointOfBoundaryToCheck == endPoint) return false;

      if (Dune::FloatCmp::eq(checkPoint, startPointOfBoundaryToCheck, 1e-8)) {
        startUOfCurve2 = _boundary.domain[0];
        return true;
      } else if (Dune::FloatCmp::eq(checkPoint, endPointOfBoundaryToCheck, 1e-8)) {
        startUOfCurve2 = _boundary.domain[1];
        return true;
      }
      return false;
    });

    if (it != state->boundaries.end()) {
      int newBoundaryIdx     = static_cast<int>(std::ranges::distance(state->boundaries.begin(), it));
      auto curveToTrace2     = state->boundaries[newBoundaryIdx].nurbsGeometry;
      auto traceCurveResult2 = traceCurveImpl(state->pointMapPtr.get(), curveToTrace2, traceCurveInput.startEdge);

      assert(Dune::FloatCmp::eq(curveToTrace(endUOfCurve1), curveToTrace2(startUOfCurve2), 1e-8));

      if (traceCurveResult2.has_value())
        return std::make_optional<std::vector<TraceCurveOutput>>(
            {{traceCurveInput.boundaryIdx, traceCurveInput.startU, endUOfCurve1, checkPoint, -1},
             {newBoundaryIdx, startUOfCurve2, traceCurveResult2->intersectU, traceCurveResult2->intersectionPoint,
              traceCurveResult2->foundOnEdge}});
    }

    return std::nullopt;
  }

  bool findNextBoundary(FindNextBoundaryLoopState* state, std::vector<Boundary>& elementBoundaries) {
    int node = state->currentNode;
    int edge = node;

    // Check if edge has no intersections and next node has one
    if (state->isCurrentlyOnNode && !state->hasIntersectionPointOnEdgeNr(edge)
        && state->hasIntersectionPointOnNodeNr(nextEntityIdx(node, 1))) {
      elementBoundaries.push_back(boundaryForEdge(state->corners, edge));

      // Update Node and Point
      node                    = nextEntityIdx(node, 1);
      state->currentNode      = node;
      state->currentNodePoint = toFloatDomain(state->corners[node]);

      return (node == state->nodeToBegin);
    }
    if (state->isCurrentlyOnNode) {
      // Determine the next edge with an IP
      for (int i = 0; i < 4; ++i) {
        if (state->hasIntersectionPointOnEdgeNr(nextEntityIdx(node, i))) {
          edge = nextEntityIdx(node, i);
          break;
        }
      }
      if (!(state->hasIntersectionPointOnEdgeNr(edge)))
        throw std::runtime_error(
            "No next edge with an intersection Point found. Maybe the IP is on an edge, which is not yet implemented");

      // Set node for the next step
      state->currentNode      = node;
      state->currentNodePoint = toFloatDomain(state->corners[node]);

      /// Find the intersection Point on the current Edge

      auto findNextBoundaryResult
          = findBoundaryThatHasIntersectionPoint(state->pointMapPtr.get(), edge, state->corners, state->boundaries);

      if (not findNextBoundaryResult.has_value()) throw std::runtime_error("findBoundaryThatHasIntersectionPoint not successful");

      auto [boundaryIndexThatHasIntersectionPoint, edgeIntersectResult] = findNextBoundaryResult.value();
      elementBoundaries.emplace_back(state->currentNodePoint, edgeIntersectResult.globalResult);

      /// Curve Tracing

      auto traceResult
          = traceCurve(state, {boundaryIndexThatHasIntersectionPoint, edgeIntersectResult.localResult, edge});

      if (not traceResult.has_value()) throw std::runtime_error("tracCurve not successful");

      for (const auto& result : traceResult.value()) {
        auto curveToTrace    = state->boundaries[result.boundaryIdx].nurbsGeometry;
        Boundary newBoundary = Boundary(curveToTrace, std::array<double, 2>{result.beginU, result.intersectU});
        elementBoundaries.push_back(newBoundary);
      }

      // Now we assume we are on an edge (has to be atm)
      state->isCurrentlyOnNode = false;
      state->currentNodePoint  = traceResult.value().back().intersectionPoint;
      state->edgeTheLoopIsOn   = traceResult.value().back().foundOnEdge;

      state->moreIntersectionPointsOnEdge
          = getNextIntersectionPointOnEdge(state->pointMapPtr.get(), state->edgeTheLoopIsOn).second;

      return false;
    }
    // There are two possibilities when the loop is not on a node, but on an edge
    // 1. we have more than one intersectionPoint on this edge (here) or not (later)
    if (!(state->isCurrentlyOnNode) && state->moreIntersectionPointsOnEdge) {
      edge = state->edgeTheLoopIsOn;

      /// Find the intersection Point on the current Edge

      auto findNextBoundaryResult
          = findBoundaryThatHasIntersectionPoint(state->pointMapPtr.get(), edge, state->corners, state->boundaries);

      if (not findNextBoundaryResult.has_value()) throw std::runtime_error("findBoundaryThatHasIntersectionPoint not successful");

      auto [boundaryIndexThatHasIntersectionPoint, edgeIntersectResult] = findNextBoundaryResult.value();
      elementBoundaries.emplace_back(state->currentNodePoint, edgeIntersectResult.globalResult);

      /// Curve Tracing

      auto traceResult
          = traceCurve(state, {boundaryIndexThatHasIntersectionPoint, edgeIntersectResult.localResult, edge});

      if (not traceResult.has_value()) throw std::runtime_error("tracCurve not successful");

      for (const auto& result : traceResult.value()) {
        auto curveToTrace    = state->boundaries[result.boundaryIdx].nurbsGeometry;
        Boundary newBoundary = Boundary(curveToTrace, std::array<double, 2>{result.beginU, result.intersectU});
        elementBoundaries.push_back(newBoundary);
      }

      // Now we assume we are on an edge (has to be atm)
      state->isCurrentlyOnNode = false;
      state->currentNodePoint  = traceResult.value().back().intersectionPoint;
      state->edgeTheLoopIsOn   = traceResult.value().back().foundOnEdge;

      state->moreIntersectionPointsOnEdge
          = getNextIntersectionPointOnEdge(state->pointMapPtr.get(), state->edgeTheLoopIsOn).second;

      return false;
    }

    if (!(state->isCurrentlyOnNode) && !(state->moreIntersectionPointsOnEdge)) {
      elementBoundaries.emplace_back(state->currentNodePoint,
                                     toFloatDomain(state->corners[nextEntityIdx(state->edgeTheLoopIsOn, 1)]));

      state->isCurrentlyOnNode = true;
      node                     = nextEntityIdx(state->edgeTheLoopIsOn, 1);

      // Update Node and Point
      state->currentNode      = node;
      state->currentNodePoint = toFloatDomain(state->corners[node]);
      state->edgeTheLoopIsOn  = -1;

      return (node == state->nodeToBegin);
    }
    __builtin_unreachable();
  }

  std::optional<std::vector<Boundary>> constructElementBoundaries(ClippingResult& result,
                                                                  std::vector<IntPoint>& corners,
                                                                  std::weak_ptr<TrimData> data) {
    try {
      std::vector<Boundary> elementBoundaries;

      // Determine the node where we begin to trace
      const auto nodeCounter = result.nodeCounter;
      auto it_nodes          = std::ranges::find_if(nodeCounter, [](auto x) { return x > 0; });

      if (it_nodes == nodeCounter.end())
        throw std::runtime_error(
            "No Curve tracing scheme implemented for elements that are only trimmed on single element edges.");

      auto nodeToBegin = std::distance(nodeCounter.begin(), it_nodes);
      auto state = std::make_unique<FindNextBoundaryLoopState>(nodeToBegin, corners, result, std::shared_ptr(data));

      // Iteratively find all the elementBoundaries
      bool found;
      do {
        found = findNextBoundary(state.get(), elementBoundaries);
      } while (!found);

      return std::make_optional<std::vector<Boundary>>(elementBoundaries);
    } catch (const std::runtime_error& e) {
      std::cout << "Trimming encountered problem: " << e.what() << "\n Element is ignored" << std::endl;
      return std::nullopt;
    }
  }

}  // namespace Dune::IGA::Impl::Trim
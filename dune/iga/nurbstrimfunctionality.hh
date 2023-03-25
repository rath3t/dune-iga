//
// Created by Henri on 14.03.2023.
//

#pragma once

namespace Dune::IGA::Impl::Trim {

  // The following numbering is used for edges and nodes:
  /*
   * 3 -- 2 -- 2
   * |         |
   * 3         1
   * |         |
   * 0 -- 0 -- 1
   */

  constexpr int clipperPrecision = 8;
  constexpr int pathSamples      = 200;
  constexpr double tolerance     = 1e-8;
  constexpr double fullElementTolerance = 1e-5;
  // constexpr double epsilonPrecision = double(16) * std::numeric_limits<double>::epsilon();

  using ctype        = double;
  using ClipperPoint = Clipper2Lib::Point<ctype>;
  using Point        = FieldVector<ctype, 2>;
  using PointVector  = std::vector<Point>;

  using CurveGeometry = NURBSPatchGeometry<1, 2>;

  // To keep track of which point belongs to which edge or node we will use a pointMap
  struct IntersectionPoint {
    Point point;
    bool alreadyVisited = false;

    explicit IntersectionPoint(const ClipperPoint& _clipperPoint) : point({_clipperPoint.x, _clipperPoint.y}){};
  };

  using IntersectionPointMap    = std::map<int, std::vector<IntersectionPoint>>;

  // Helper functions
  double distance(const ClipperPoint p1, const ClipperPoint p2) {
    return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
  }

  bool approxSamePoint(const Point& a, const ClipperPoint& b, const double tol = tolerance) {
    return Dune::FloatCmp::eq(a, Point{b.x, b.y}, tol);
  }

  bool approxSamePoint(Point a, Point b, const double tol = tolerance) { return Dune::FloatCmp::eq(a, b, tol); }

  bool approxSamePoint(ClipperPoint a, ClipperPoint b, const double tol = tolerance) {
    return Dune::FloatCmp::eq(distance(a, b), 0.0, tol);
  }

  int getEdgeOrientation(const int edge) { return (edge == 0 || edge == 2) ? 0 : 1; }

  int iPlusX(int i, int x) { return (i + x) % 4; }

  std::array<Point, 2> curveStartEndPoint(const CurveGeometry& curve) {
    std::array<int, 2> indices{0, static_cast<int>(curve.patchData_->controlPoints.size()[0]) - 1};
    return {curve.patchData_->controlPoints.get({indices[0]}).p, curve.patchData_->controlPoints.get({indices[1]}).p};
  }

  PointVector getElementCorners(const auto& element) {
    auto geo = element.geometry();

    // The corners are sorted counter-clockwise to get a closed loop, also see top of this file for numbering
    return {geo.corner(0), geo.corner(1), geo.corner(3), geo.corner(2)};
  }

  Clipper2Lib::PathD getElementEdgesFromElementCorners(const PointVector& corners) {
    Clipper2Lib::PathD edges;
    for (int i = 0; i < 4; ++i)
      edges.emplace_back(corners[i][0], corners[i][1]);
    return edges;
  }

  Clipper2Lib::PathD getElementEdges(const auto& element) {
    auto corners = getElementCorners(element);
    return getElementEdgesFromElementCorners(corners);
  }

  Clipper2Lib::PathsD getClip(std::vector<Boundary>& boundaries) {
    Clipper2Lib::PathD pathD;
    for (auto& boundary : boundaries)
      for (auto& point : boundary.path(pathSamples))
        pathD.push_back(point);

    // Trim Patch, so we don't get duplicate intersection Points on corners
    pathD = Clipper2Lib::TrimCollinear(pathD, clipperPrecision);

    return {pathD};
  }

  /// Checks if a point p is on the line spanned by the points a and b within a given tolerance
  bool pointOnLine(const ClipperPoint& p, const ClipperPoint& a, const ClipperPoint& b, const double tol = tolerance) {
    // https://stackoverflow.com/a/17693146
    return Dune::FloatCmp::eq(distance(a, p) + distance(b, p), distance(a, b), tol);
  }

  std::pair<bool, int> pointOnAnyEdge(const ClipperPoint& point, const auto& elementEdges) {
    for (int i = 0; i < 4; ++i)
      if (pointOnLine(point, elementEdges[i], elementEdges[iPlusX(i, 1)])) return std::make_pair(true, i);

    return std::make_pair(false, -1);
  }

  std::pair<bool, int> pointOnAnyVertex(const ClipperPoint& point, const PointVector& corners) {
    for (size_t i = 0; i < corners.size(); ++i)
      if (approxSamePoint(corners[i], point)) return std::make_pair(true, static_cast<int>(i));

    return std::make_pair(false, -1);
  }

  // Sorts the Intersection Points counter-clockwise around the element (crucial step!)
  void sortIntersectionPoints(const std::shared_ptr<IntersectionPointMap>& pointMapPtr, const int edge) {
    auto orientation = getEdgeOrientation(edge);
    auto cmp         = [orientation, edge](const IntersectionPoint& a, const IntersectionPoint& b) {
      return (edge == 0 || edge == 1) ? a.point[orientation] < b.point[orientation]
                                              : a.point[orientation] > b.point[orientation];
    };
    std::ranges::sort((*pointMapPtr)[edge], cmp);
  }

  bool pointInElement(const Point& point, const Clipper2Lib::PathD& elementEdges) {
    auto res = Clipper2Lib::PointInPolygon({point[0], point[1]}, elementEdges);
    return (res == Clipper2Lib::PointInPolygonResult::IsInside);
  }

  struct ClippingResult {
    std::shared_ptr<IntersectionPointMap> pointMapPtr;
    std::array<int, 4> edgeCounter;
    std::array<int, 4> nodeCounter;
  };

  std::pair<ElementTrimFlag, std::optional<ClippingResult>> clipElement(const auto& element,
                                                                        Clipper2Lib::PathsD& clip) {
    // Prepare Result
    ElementTrimFlag trimFlag;

    auto corners      = getElementCorners(element);
    auto elementEdges = getElementEdges(element);

    Clipper2Lib::RectD elementRect{corners[0][0], corners[1][1], corners[1][0], corners[3][1]};

    // Clip with RectClip (more efficient than normal Intersect Clip, works because our ParameterSpace Grid
    // has always rectangular elements
    Clipper2Lib::PathsD clippedEdges = Clipper2Lib::RectClip(elementRect, clip, clipperPrecision);

    // If the clippedEdges are empty, this means the element is outside of the clip -> empty
    // And if the area of the clip is equal to the area of the element -> full
    if (clippedEdges.empty())
      trimFlag = ElementTrimFlag::empty;
    else {
      if (FloatCmp::eq(Clipper2Lib::Area(clippedEdges), element.geometry().volume(), fullElementTolerance))
        trimFlag = ElementTrimFlag::full;
      else
        trimFlag = ElementTrimFlag::trimmed;
    }

    if (trimFlag != ElementTrimFlag::trimmed) return std::make_pair(trimFlag, std::nullopt);

    // If the code is here → Elements are trimmed
    assert(clippedEdges.size() == 1);
    clippedEdges.front() = Clipper2Lib::TrimCollinear(clippedEdges.front(), clipperPrecision);

    // Create a local pointMap to store the intersection points in regard to their corresponding edges or nodes
    // Node Intersection Points are stored with indices 10 - 13 & edges with 0 - 3
    IntersectionPointMap pointMap;

    for (auto& point : clippedEdges[0]) {
      // Check if point is on any Edge
      auto [isOnAnyEdge, edgeThePointIsOn] = pointOnAnyEdge(point, elementEdges);
      auto [isOnAnyNode, nodeThePointIsOn] = pointOnAnyVertex(point, corners);

      // If is onAnyNode is true, then isOnAnyEdge has to be also true
      if (isOnAnyNode)
        pointMap[nodeThePointIsOn + 10].emplace_back(point);
      else if (isOnAnyEdge)
        pointMap[edgeThePointIsOn].emplace_back(point);
    }

    // This keeps track of how many intersection Points there are at a given edge
    std::array<int, 4> edgeCounter{0, 0, 0, 0};
    std::array<int, 4> nodeCounter{0, 0, 0, 0};

    for (const auto& [key, value] : pointMap) {
      if (key < 4)
        edgeCounter[key] = static_cast<int>(value.size());
      else
        nodeCounter[key - 10] = static_cast<int>(value.size());
    }

    // Make a shared Pointer out of pointMap
    std::shared_ptr<IntersectionPointMap> pointMapPtr = std::make_shared<IntersectionPointMap>(pointMap);

    // Sort IntersectionPoints counter-clockwise
    for (int i = 0; i < 4; ++i)
      sortIntersectionPoints(pointMapPtr, i);

    return std::make_pair(trimFlag,
                          std::make_optional<ClippingResult>({std::move(pointMapPtr), edgeCounter, nodeCounter}));
  }

  bool areLinesParallel(CurveGeometry& line1, CurveGeometry& line2) {
    if (line1.degree()[0] > 1 || line2.degree()[0] > 1) return false;

    auto [p1, p2] = curveStartEndPoint(line1);
    auto [p3, p4] = curveStartEndPoint(line2);

    // Compute direction vectors of both curves
    auto d1 = p2 - p1;
    auto d2 = p4 - p3;

    auto cos_angle = Dune::dot(d1, d2) / (d1.two_norm() * d2.two_norm());
    auto angle     = std::acos(std::min(std::max(cos_angle, -1.0), 1.0));

    // Check if the angle is small enough to consider the curves on the same line (or parallel)
    return (Dune::FloatCmp::eq(angle, 0.0, tolerance) || Dune::FloatCmp::eq(angle, std::numbers::pi, tolerance));
  }

  struct IntersectionResult {
    double localResult{};
    Point globalResult;
    bool found = false;
  };

  IntersectionResult newtonRaphson(CurveGeometry& curve, Point& intersectionPoint) {
    const double objectiveFctTolerance  = 1e-4;
    const int max_iterations            = 250;

    // There are some easy cases where the IP is at the back or at the front of the domain, check these cases first
    auto [start, end] = curveStartEndPoint(curve);

    if (Dune::FloatCmp::eq(start, intersectionPoint, tolerance))
      return {curve.domain()[0].front(), start, true};
    if (Dune::FloatCmp::eq(end, intersectionPoint, tolerance))
      return {curve.domain()[0].back(), end, true};

    // At first, we assume that the intersectionPoint is in some definition exact and the corresponding u has to be
    // found
    Dune::FieldVector<double, 1> u{curve.domainMidPoint()[0]};
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

      ++i;
      if (i > max_iterations) return {u[0], intersectionPoint, false};

      // Clamp result, this might be already an indicator that there is no solution
      if (Dune::FloatCmp::gt(u[0], curve.patchData_->knotSpans[0].back())) u[0] = curve.patchData_->knotSpans[0].back();
      if (Dune::FloatCmp::lt(u[0], curve.patchData_->knotSpans[0].front()))
        u[0] = curve.patchData_->knotSpans[0].front();

    } while (dGlobal.two_norm() > objectiveFctTolerance);

    double localResult = u[0];
    auto globalResult  = curve(localResult);

    return {localResult, globalResult, true};
  }

  Boundary boundaryForEdge(const PointVector& corners, int edge) {
    return Boundary{corners[iPlusX(edge, 0)], corners[iPlusX(edge, 1)]};
  }

  std::pair<Point, bool> getNextIntersectionPointOnEdge(const std::shared_ptr<IntersectionPointMap>& pointMapPtr, int edge) {
    if ((*pointMapPtr)[edge].empty()) return {Point(), false};

    for (const IntersectionPoint& ip : (*pointMapPtr)[edge]) {
      if (!ip.alreadyVisited) return {ip.point, true};
    }
    return {Point(), false};
  }

  void markIntersectionPointAsVisited(const std::shared_ptr<IntersectionPointMap>& pointMapPtr, int edge, Point& point) {
    auto it = std::find_if((*pointMapPtr)[edge].begin(), (*pointMapPtr)[edge].end(), [&point](IntersectionPoint x) {
      return (Dune::FloatCmp::eq(x.point, point));
    });
    assert(it != (*pointMapPtr)[edge].end());
    (*it).alreadyVisited = true;
  }

  std::tuple<int, IntersectionResult, bool> findBoundaryThatHasIntersectionPoint(
      const std::shared_ptr<IntersectionPointMap>& pointMapPtr, int edge, PointVector& corners,
      std::vector<Boundary>& globalBoundaries) {
    auto edgeCurve               = boundaryForEdge(corners, edge);
    Point intersectionPointGuess = getNextIntersectionPointOnEdge(pointMapPtr, edge).first;

    for (int i = 0; auto& boundary : globalBoundaries) {
      if (areLinesParallel(boundary.nurbsGeometry, edgeCurve.nurbsGeometry)) {
        ++i;
        continue;
      }

      auto edgeIntersectResult = newtonRaphson(boundary.nurbsGeometry, intersectionPointGuess);

      if (edgeIntersectResult.found) {
        markIntersectionPointAsVisited(pointMapPtr, edge, intersectionPointGuess);
        return {i, edgeIntersectResult, true};
      }
      ++i;
    }
    return {-1, IntersectionResult(), false};
  }

  struct TraceCurveResult {
    double intersectU;
    Point intersectionPoint;
    int foundOnEdge;
  };

  std::optional<TraceCurveResult> traceCurveImpl(const std::shared_ptr<IntersectionPointMap>& pointMapPtr,
                                                  CurveGeometry& curveToTrace, int startEdge) {

    // Look through edges and determine the next Intersection Point (the IPs are sorted counter-clockwise)

    for (int i = startEdge; i < startEdge + 4; ++i) {
      auto currentEdge = iPlusX(startEdge, i);
      auto [intersectionPoint, hasIntersectionPoint] = getNextIntersectionPointOnEdge(pointMapPtr, currentEdge);
      if (hasIntersectionPoint) {
        // Check if the curveToTrace goes through the intersectinPoint
        auto [localResult, globalResult, found] = newtonRaphson(curveToTrace, intersectionPoint);
        if (found) {
          markIntersectionPointAsVisited(pointMapPtr, currentEdge, intersectionPoint);
          return std::make_optional<TraceCurveResult>({localResult, globalResult, currentEdge});
        }
      }
    }
    return std::nullopt;
  }

  struct FindNextBoundaryLoopState {
    // Variables that need to be set
    const int nodeToBegin;
    PointVector corners;
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
    constexpr bool hasIntersectionPointOnEdgeNr(int nr) { return edgeCounter[nr] > 0; };
    constexpr bool hasIntersectionPointOnNodeNr(int nr) { return nodeCounter[nr] > 0; };

    FindNextBoundaryLoopState(int _nodeToBegin, PointVector& _corners, ClippingResult& _result,
                              std::vector<Boundary>& _boundaries)
        : nodeToBegin(_nodeToBegin),
          corners(_corners),
          pointMapPtr(std::move(_result.pointMapPtr)),
          edgeCounter(_result.edgeCounter),
          nodeCounter(_result.nodeCounter),
          boundaries(_boundaries) {
      currentNode      = nodeToBegin;
      currentNodePoint = corners[currentNode];
    };

    FindNextBoundaryLoopState() = delete;
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

  std::optional<std::vector<TraceCurveOutput>> traceCurve(const std::unique_ptr<FindNextBoundaryLoopState>& state,
                                                          TraceCurveInput traceCurveInput) {
    // Check the boundary that was given in traceInput
    CurveGeometry curveToTrace = state->boundaries[traceCurveInput.boundaryIdx].nurbsGeometry;
    auto traceCurveResult = traceCurveImpl(state->pointMapPtr, curveToTrace, traceCurveInput.startEdge);

    if (traceCurveResult.has_value())
      return std::make_optional<std::vector<TraceCurveOutput>>(
          {{traceCurveInput.boundaryIdx, traceCurveInput.startU, traceCurveResult->intersectU, traceCurveResult->intersectionPoint, traceCurveResult->foundOnEdge}});

    // Now check if the start or endPoint of the curve is in the element
    auto elementEdges = getElementEdgesFromElementCorners(state->corners);
    int status        = -1;  // -1 means not found, 0 means start, 1 means end
    auto [start, end] = curveStartEndPoint(curveToTrace);

    // TODO Check that start and End are not known IntersectionPoints or are on an edge (this is important for elements on the boundary)

    if (pointInElement(end, elementEdges))
      status = 1;
    else if (pointInElement(start, elementEdges))
      status = 0;

    if (status == -1) return std::nullopt;

    // Now we check if there is another boundary that begins at this point and trace this potential curve
    Point checkPoint = (status == 1) ? end : start;
    double endU = (status == 1) ? curveToTrace.domain()[0].back() : curveToTrace.domain()[0].front();

    double startU;
    auto it = std::ranges::find_if(state->boundaries, [&checkPoint, &startU](auto& _boundary) {
      auto [startBi, endBi] = curveStartEndPoint(_boundary.nurbsGeometry);
      if (Dune::FloatCmp::eq(checkPoint, startBi, 1e-8)) {
        startU = _boundary.domain[0];
        return true;
      } else if (Dune::FloatCmp::eq(checkPoint, endBi, 1e-8)) {
        startU = _boundary.domain[1];
        return true;
      }
      return false;
    });
    if (it != state->boundaries.end()) {
      int newBoundaryIdx = static_cast<int>(std::ranges::distance(state->boundaries.begin(), it));
      curveToTrace       = state->boundaries[newBoundaryIdx].nurbsGeometry;

     auto traceCurveResult2 = traceCurveImpl(state->pointMapPtr, curveToTrace, traceCurveInput.startEdge);

      if (traceCurveResult2.has_value())
        return std::make_optional<std::vector<TraceCurveOutput>>(
            {{traceCurveInput.boundaryIdx, traceCurveInput.startU, endU, checkPoint, -1},
             {newBoundaryIdx, startU, traceCurveResult2->intersectU, traceCurveResult2->intersectionPoint, traceCurveResult2->foundOnEdge}});
    }

    return std::nullopt;
  }

  bool findNextBoundary(const std::unique_ptr<FindNextBoundaryLoopState>& state,
                        std::vector<Boundary>& elementBoundaries) {
    int node = state->currentNode;
    int edge = node;

    // Check if edge has no intersections and next node has one
    if (state->isCurrentlyOnNode && !state->hasIntersectionPointOnEdgeNr(edge)
        && state->hasIntersectionPointOnNodeNr(iPlusX(node, 1))) {
      elementBoundaries.push_back(boundaryForEdge(state->corners, edge));

      // Update Node and Point
      node                    = iPlusX(node, 1);
      state->currentNode      = node;
      state->currentNodePoint = state->corners[node];

      return (node == state->nodeToBegin);
    }
    if (state->isCurrentlyOnNode) {
      // Determine the next edge with an IP
      for (int i = 0; i < 4; ++i) {
        if (state->hasIntersectionPointOnEdgeNr(iPlusX(node, i))) {
          edge = iPlusX(node, i);
          break;
        }
      }
      if (!(state->hasIntersectionPointOnEdgeNr(edge)))
        throw std::runtime_error(
            "No next edge with an intersection Point found. Maybe the IP is on an edge, which is not yet implemented");

      // Set node for the next step
      state->currentNode      = node;
      state->currentNodePoint = state->corners[node];

      ////////
      // Find the intersection Point on the current Edge
      ////////

      auto [boundaryIndexThatHasIntersectionPoint, edgeIntersectResult, found]
          = findBoundaryThatHasIntersectionPoint(state->pointMapPtr, edge, state->corners, state->boundaries);

      // There has to be an intersectionPoint
      if (!found) throw std::runtime_error("No Intersection Point found on a edge where there has to be one.");

      // Make a parametrisation of the section before the intersectionPoint
      elementBoundaries.emplace_back(state->currentNodePoint, edgeIntersectResult.globalResult);

      ////////
      // Curve Tracing
      ////////

      auto elementEdges = getElementEdgesFromElementCorners(state->corners);
      auto traceResult  = traceCurve(state, {boundaryIndexThatHasIntersectionPoint, edgeIntersectResult.localResult, edge});

      if (!(traceResult.has_value()))
        throw std::runtime_error("No Intersection Point found on a edge where there has to be one.");

      for (const auto& result : *traceResult) {
        // After we found the next IntersectionPoint on this curve, add the so found part of the elementBoundary
        auto curveToTrace    = state->boundaries[result.boundaryIdx].nurbsGeometry;
        Boundary newBoundary = Boundary(curveToTrace, std::array<double, 2>{result.beginU, result.intersectU});
        elementBoundaries.push_back(newBoundary);
      }

      // Now we assume we are on an edge (has to be atm)
      state->isCurrentlyOnNode = false;
      state->currentNodePoint  = traceResult.value().back().intersectionPoint;
      state->edgeTheLoopIsOn   = traceResult.value().back().foundOnEdge;

      state->moreIntersectionPointsOnEdge
          = getNextIntersectionPointOnEdge(state->pointMapPtr, state->edgeTheLoopIsOn).second;

      return false;
    }
    // There are two possibilities when the loop is not on a node, but on an edge
    // 1. we have more than one intersectionPoint on this edge (here) or not (later)
    if (!(state->isCurrentlyOnNode) && state->moreIntersectionPointsOnEdge) {
      edge = state->edgeTheLoopIsOn;

      ////////
      // Find the intersection Point on the current Edge
      ////////

      auto [boundaryIndexThatHasIntersectionPoint, edgeIntersectResult, found]
          = findBoundaryThatHasIntersectionPoint(state->pointMapPtr, edge, state->corners, state->boundaries);

      // There has to be an intersectionPoint otherwise we wouldn't be here
      if (!found) throw std::runtime_error("No Intersection Point found on a edge where there has to be one.");

      // Make a parametrisation of the section before the intersectionPoint
      elementBoundaries.emplace_back(state->currentNodePoint, edgeIntersectResult.globalResult);

      ////////
      // Curve Tracing
      ////////

      auto traceResult = traceCurve(state, {boundaryIndexThatHasIntersectionPoint, edgeIntersectResult.localResult, edge});

      if (!(traceResult.has_value()))
        throw std::runtime_error("No Intersection Point found on a edge where there has to be one.");

      for (const auto& result : *traceResult) {
        // After we found the next IntersectionPoint on this curve, add the so found part of the elementBoundary
        auto curveToTrace    = state->boundaries[result.boundaryIdx].nurbsGeometry;
        Boundary newBoundary = Boundary(curveToTrace, std::array<double, 2>{result.beginU, result.intersectU});
        elementBoundaries.push_back(newBoundary);
      }

      // Now we assume we are on an edge (has to be atm)
      state->isCurrentlyOnNode = false;
      state->currentNodePoint  = traceResult.value().back().intersectionPoint;
      state->edgeTheLoopIsOn   = traceResult.value().back().foundOnEdge;

      state->moreIntersectionPointsOnEdge
          = getNextIntersectionPointOnEdge(state->pointMapPtr, state->edgeTheLoopIsOn).second;

      return false;
    }

    if (!(state->isCurrentlyOnNode) && !(state->moreIntersectionPointsOnEdge)) {
      elementBoundaries.emplace_back(state->currentNodePoint, state->corners[iPlusX(state->edgeTheLoopIsOn, 1)]);

      state->isCurrentlyOnNode = true;
      node                     = iPlusX(state->edgeTheLoopIsOn, 1);

      // Update Node and Point
      state->currentNode      = node;
      state->currentNodePoint = state->corners[node];
      state->edgeTheLoopIsOn  = -1;

      return (node == state->nodeToBegin);
    }
    __builtin_unreachable();
  }

  std::optional<std::vector<Boundary>> constructElementBoundaries(ClippingResult& result, PointVector& corners,
                                                                  std::vector<Boundary>& boundaries) {
    std::vector<Boundary> elementBoundaries;

    // Determine the node where we begin to trace
    const std::array<int, 4> nodeCounter = result.nodeCounter;
    auto it_nodes                        = std::ranges::find_if(nodeCounter, [](auto x) { return x > 0; });

    if (it_nodes == nodeCounter.end())
      throw std::runtime_error("No Curve tracing scheme implemented for elements that are only trimmed on single element edges.");


    auto nodeToBegin = std::distance(nodeCounter.begin(), it_nodes);
    auto state       = std::make_unique<FindNextBoundaryLoopState>(nodeToBegin, corners, result, boundaries);

    try {
      // Iteratively find all the elementBoundaries
      bool found;
      do {
        found = findNextBoundary(state, elementBoundaries);
      } while (!found);
      return std::make_optional<std::vector<Boundary>>(elementBoundaries);
    } catch (const std::runtime_error& e) {
      std::cout << "Trimming encountered problem: " << e.what() << "\n Element is ignored" << std::endl;
      return std::nullopt;
    }
  }

}  // namespace Dune::IGA::Impl::Trim
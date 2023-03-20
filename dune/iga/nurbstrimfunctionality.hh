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

  // Helper functions

  bool approxSamePoint(Point a, Point b, const double tol = 1e-2) { return Dune::FloatCmp::eq(a, b, tol); }

  double distance(ClipperPoint p1, ClipperPoint p2) {
    return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
  }

  bool approxSamePoint(const Point& a, const ClipperPoint& b, const double tol = 1e-8) {
    return Dune::FloatCmp::eq(a, Point{b.x, b.y}, tol);
  }

  bool approxSamePoint(ClipperPoint a, ClipperPoint b, const double tol = 1e-8) {
    return Dune::FloatCmp::eq(distance(a, b), 0.0, tol);
  }

  int getEdgeOrientation(int edge) { return (edge == 0 || edge == 2) ? 0 : 1; }

  int iPlusX(int i, int x) { return (i + x) % 4; }

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

  Clipper2Lib::PathsD getClip(std::vector<Boundary>& boundaries, const int pathSamples = 200,
                              const int clipperPrecision = 8) {
    Clipper2Lib::PathD pathD;
    for (auto& boundary : boundaries)
      for (auto& point : boundary.path(pathSamples))
        pathD.push_back(point);

    // Trim Patch, so we don't get duplicate intersection Points on corners
    pathD = Clipper2Lib::TrimCollinear(pathD, clipperPrecision);

    return {pathD};
  }

  /// Checks if a point p is on the line spanned by the points a and b within a given tolerance
  bool pointOnLine(const ClipperPoint& p, const ClipperPoint& a, const ClipperPoint& b, const double tolerance = 1e-8) {
    // https://stackoverflow.com/a/17693146
    return Dune::FloatCmp::eq(distance(a, p) + distance(b, p), distance(a, b), tolerance);
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

  void sortIntersectionPoints(const IntersectionPointMapPtr& pointMapPtr, int edge) {
    auto orientation = getEdgeOrientation(edge);
    auto cmp         = [orientation, edge](const IntersectionPoint& a, const IntersectionPoint& b) {
      return (edge == 0 || edge == 1) ? a.point[orientation] < b.point[orientation]
                                              : a.point[orientation] > b.point[orientation];
    };
    std::ranges::sort((*pointMapPtr)[edge], cmp);
  }

  struct ClippingResult {
    IntersectionPointMapPtr pointMapPtr;
    std::array<int, 4> edgeCounter;
    std::array<int, 4> nodeCounter;
  };

  std::pair<ElementTrimFlag, std::optional<ClippingResult>> clipElement(const auto& element,
                                                                        Clipper2Lib::PathsD& clip) {
    constexpr int clipperPrecision = 8;
    constexpr double tolerance     = 1e-8;

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
      if (FloatCmp::eq(Clipper2Lib::Area(clippedEdges), element.geometry().volume(), tolerance))
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
      else if (key > 9)
        nodeCounter[key - 10] = static_cast<int>(value.size());
    }

    // Make a shared Pointer out of pointMap
    IntersectionPointMapPtr pointMapPtr = std::make_shared<IntersectionPointMap>(pointMap);

    // Sort IntersectionPoints counter-clockwise
    for (int i = 0; i < 4; ++i)
      sortIntersectionPoints(pointMapPtr, i);

    return std::make_pair(trimFlag,
                          std::make_optional<ClippingResult>({std::move(pointMapPtr), edgeCounter, nodeCounter}));
  }

  bool areEdgesParallel(CurveGeometry& curve1, CurveGeometry& curve2) {
    if (curve1.degree()[0] > 1 || curve2.degree()[0] > 1) return false;

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
    return (Dune::FloatCmp::eq(angle, 0.0, tolerance) || Dune::FloatCmp::eq(angle, std::numbers::pi, tolerance));
  }

  struct IntersectionResult {
    double localResult{};
    Point globalResult;
    bool found = false;
  };

  IntersectionResult newtonRaphson(CurveGeometry& curve, Point& intersectionPoint) {
    const double pointEqualityTolerance = 1e-7;
    const double objectiveFctTolerance  = 1e-4;

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
      if (Dune::FloatCmp::lt(u[0], curve.patchData_->knotSpans[0].front()))
        u[0] = curve.patchData_->knotSpans[0].front();

    } while (dGlobal.two_norm() > objectiveFctTolerance);

    double localResult = u[0];
    auto globalResult  = curve(localResult);

    return {localResult, globalResult, true};
  }

  Boundary boundaryForEdge(const PointVector& corners, int edge) {
    PointVector v{corners[iPlusX(edge, 0)], corners[iPlusX(edge, 1)]};
    return Boundary{v};
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

  struct FindNextBoundaryLoopState {
    // Variables that need to be set
    const int nodeToBegin;
    PointVector corners;
    std::vector<Boundary> boundaries;
    IntersectionPointMapPtr pointMapPtr;
    const std::array<int, 4> edgeCounter;
    const std::array<int, 4> nodeCounter;

    // Loop States
    bool isCurrentlyOnNode            = true;
    bool moreIntersectionPointsOnEdge = false;
    int edgeTheLoopIsOn               = -1;  // This variable is only filed when !loopIsCurrentlyOnNode

    // Current positions (are set through the constructor
    int currentNode;
    Point currentNodePoint;

    // Lambdas
    bool hasIntersectionPointOnEdgeNr(int nr) { return edgeCounter[nr] > 0; };
    bool hasIntersectionPointOnNodeNr(int nr) { return nodeCounter[nr] > 0; };

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

  bool findNextBoundary(const std::unique_ptr<FindNextBoundaryLoopState>& state,
                        std::vector<Boundary>& elementBoundaries) {
    int node = state->currentNode;
    int edge = node;

    // Check if edge has no intersections and next node has one
    if (state->isCurrentlyOnNode && !state->hasIntersectionPointOnEdgeNr(edge)
        && state->hasIntersectionPointOnNodeNr(iPlusX(node, 1))) {
      elementBoundaries.push_back(boundaryForEdge(state->corners, edge));

      // Update Node and Point
      node = iPlusX(node, 1);
      state->currentNode      = node;
      state->currentNodePoint = state->corners[node];

      return (node == state->nodeToBegin);
    }

    if (state->isCurrentlyOnNode && !state->hasIntersectionPointOnEdgeNr(edge)
        && !(state->hasIntersectionPointOnNodeNr(iPlusX(node, 1)))) {
      std::cerr << "Not yet implemented" << std::endl;

      // Todo: So lange durch die Nodes loopen, bis ein Intersection Point gefunden wurde, aufpassen dass dieser
      // auch der StartNode sein kann
      return false;
    }
    if (state->isCurrentlyOnNode) {
      // If we are here, then there is an intersectionPoint on the current edge, we have to find the curve
      // that is intersecting there and then trace it_nodes to find the next intersection Point on the next edge
      // (or node)
      assert(state->hasIntersectionPointOnEdgeNr(edge));

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
      PointVector newEdgeParametrisation{state->currentNodePoint, edgeIntersectResult.globalResult};
      elementBoundaries.emplace_back(newEdgeParametrisation);

      ////////
      // Curve Tracing
      ////////

      // Then we need the parametrisation for the curve to come
      Boundary boundaryToTrace   = state->boundaries[boundaryIndexThatHasIntersectionPoint];
      CurveGeometry curveToTrace = boundaryToTrace.nurbsGeometry;

      // Call traceCurve
      // Todo: as for now we are only checking against edge Intersection Points, next intersection could be on a
      // node (highly unlikely but with bad luck due to precision in Intersection Point, an intersection could be
      // registered on node, instead of an edge)
      double beginU = edgeIntersectResult.localResult;
      auto [intersectU, foundInterSectionPoint, foundOnEdgeNr]
          = traceCurve(state->pointMapPtr, state->corners, curveToTrace, beginU);

      // After we found the next IntersectionPoint on this curve, add the so found part of the elementBoundary
      Boundary newBoundary = Boundary(curveToTrace, std::array<double, 2>{beginU, intersectU});
      elementBoundaries.push_back(newBoundary);

      // Now we assume we are on an edge (has to be atm)
      state->isCurrentlyOnNode = false;
      state->currentNodePoint  = foundInterSectionPoint;
      state->edgeTheLoopIsOn   = foundOnEdgeNr;

      state->moreIntersectionPointsOnEdge
          = getNextIntersectionPointOnEdge(state->pointMapPtr, state->edgeTheLoopIsOn).second;

      return false;
    }
    // There are two possibilities when the loop is not on a node, but on a edge+
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
      PointVector newEdgeParametrisation{state->currentNodePoint, edgeIntersectResult.globalResult};
      elementBoundaries.emplace_back(newEdgeParametrisation);

      ////////
      // Curve Tracing
      ////////

      // Then we need the parametrisation for the curve to come
      Boundary boundaryToTrace   = state->boundaries[boundaryIndexThatHasIntersectionPoint];
      CurveGeometry curveToTrace = boundaryToTrace.nurbsGeometry;

      // Call traceCurve
      double beginU = edgeIntersectResult.localResult;
      auto [intersectU, foundInterSectionPoint, foundOnEdgeNr]
          = traceCurve(state->pointMapPtr, state->corners, curveToTrace, beginU);

      // After we found the next IntersectionPoint on this curve, add the so found part of the elementBoundary
      Boundary newBoundary = Boundary(curveToTrace, std::array<double, 2>{beginU, intersectU});
      elementBoundaries.push_back(newBoundary);

      // Now for now we assume we are still on an edge (Todo: again check against node Intersection Points)
      state->isCurrentlyOnNode = false;
      state->currentNodePoint  = foundInterSectionPoint;
      state->edgeTheLoopIsOn   = foundOnEdgeNr;

      state->moreIntersectionPointsOnEdge
          = getNextIntersectionPointOnEdge(state->pointMapPtr, state->edgeTheLoopIsOn).second;

      return false;
    }

    if (!(state->isCurrentlyOnNode) && !(state->moreIntersectionPointsOnEdge)) {
      PointVector v1{state->currentNodePoint, state->corners[iPlusX(state->edgeTheLoopIsOn, 1)]};
      elementBoundaries.emplace_back(v1);

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

  std::optional<std::vector<Boundary>> constructElementBoundaries(ClippingResult& result, PointVector& corners, std::vector<Boundary>& boundaries) {
    std::vector<Boundary> elementBoundaries;

    // Determine the node where we begin to trace
    const std::array<int, 4> nodeCounter = result.nodeCounter;
    auto it_nodes = std::ranges::find_if(nodeCounter, [](auto x) { return x > 0; });

    if (it_nodes == nodeCounter.end()) {
      std::cerr << "No Curve tracing scheme implemented for elements that are only trimmed on single element edges. Element is ignored" << std::endl;
      return std::nullopt;
    }

    auto nodeToBegin = std::distance(nodeCounter.begin(), it_nodes);
    auto state = std::make_unique<FindNextBoundaryLoopState>(nodeToBegin, corners, result, boundaries);

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
// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "nurbstrimboundary.hh"

#include <clipper2/clipper.core.h>
#include <clipper2/clipper.h>

// #include "dune/iga/geometry/closestpointprojection.hh"
#include "dune/iga/nurbspatchgeometry.hh"
#include <dune/grid/yaspgrid.hh>

namespace Dune::IGA {
enum class ElementTrimFlag
{
  full,
  empty,
  trimmed
};
}

namespace Dune::IGA::Trim {

struct ElementBoundaries
{
  std::vector<Boundary> outerBoundaries;
  std::optional<std::vector<std::vector<Boundary>>> innerBoundaries;
};

// The following numbering is used for edges and nodes:
/*
 * 3 -- 2 -- 2
 * |         |
 * 3         1
 * |         |
 * 0 -- 0 -- 1
 */

template <typename intType_ = int64_t, int sc_ = 12>
class NURBSPatchTrimmer
{
public:
  static constexpr auto dim  = 2;
  static constexpr int scale = sc_;
  inline static double scaleFactor() {
    return std::pow(10, scale);
  };

  /// Parameters
  static constexpr int pathSamples  = 800;
  static constexpr double tolerance = 1e-8;

  /// How far the curve to point distance is allowed to be, to be considered "on" the curve
  static constexpr double gapTolerance = 1e-3;

  using intType      = intType_;
  using ClipperPoint = Clipper2Lib::Point<intType>;
  using ClipperPath  = Clipper2Lib::Path<intType>;
  using ClipperPaths = Clipper2Lib::Paths<intType>;

  using Point    = FieldVector<double, dim>;
  using IntPoint = FieldVector<intType, dim>;

  using CurveGeometry             = NURBSPatchGeometry<1, dim>;
  using TensorGrid                = Dune::YaspGrid<dim, Dune::TensorProductCoordinates<intType, dim>>;
  using TensorGridView            = TensorGrid::LeafGridView;
  using TensorGridElementGeometry = Dune::Entity<0, dim, TensorGrid, Dune::YaspEntity>::Geometry;

  NURBSPatchTrimmer(const TensorGridView& gridView, TrimData* trimData)
      : tensorGridView_(gridView),
        indexSet_(gridView.indexSet()),
        trimData_(trimData),
        globalBoundaries_(extractBoundaries(trimData)),
        clip_(getClip(trimData)) {
    setUp();
  }

  NURBSPatchTrimmer()                         = delete;
  NURBSPatchTrimmer(const NURBSPatchTrimmer&) = delete;

  /////////////////
  /// Public Interface
  /////////////////

  std::tuple<ElementTrimFlag, std::optional<ElementBoundaries>, bool> trimElement(unsigned int directIndex) {
    if (trimFlags_[directIndex] != ElementTrimFlag::trimmed)
      return std::make_tuple(trimFlags_[directIndex], std::nullopt, false);
    else {
      auto elementBoundaries = constructElementBoundaries(directIndex);
      if (elementBoundaries)
        return std::make_tuple(ElementTrimFlag::trimmed, elementBoundaries.value(), false);
      else
        return std::make_tuple(ElementTrimFlag::empty, std::nullopt, true);
    }
  }

private:
  // Forward declare some structs and data types
  struct ClippingResult;
  struct InnerLoops;
  struct FindNextBoundaryLoopState;
  struct IntersectionPoint;
  struct FindNextBoundaryResult;
  struct IntersectionResult;
  struct TraceCurveResult;
  struct TraceCurveOutput;
  struct TraceCurveInput;

  using IntersectionPointMap = std::array<std::vector<IntersectionPoint>, 4>;

  // Private data members
  const TensorGridView& tensorGridView_;
  const TensorGridView::IndexSet& indexSet_;
  TrimData* trimData_;
  std::vector<Boundary> globalBoundaries_;
  ClipperPaths clip_;
  std::vector<ElementTrimFlag> trimFlags_;
  std::vector<std::optional<ClippingResult>> clipResults_;
  std::vector<TensorGridElementGeometry> elementGeometries_;

  /////////////////
  /// Private trim methods
  /////////////////

  // starting point for trim functionality
  void setUp() {
    // Reserve space for vectors
    auto n = tensorGridView_.size(0);
    trimFlags_.reserve(n);
    elementGeometries_.reserve(n);

    clipResults_.resize(n);
    for (auto&& element : elements(tensorGridView_)) {
      elementGeometries_.push_back(element.geometry());
      clipElement(element);
    }
  }

  void findInnerLoops(std::optional<InnerLoops>& innerLoops, const ClipperPaths& clippedEdges) const {
    assert(clip_.size() > 1);
    assert(clippedEdges.size() > 1);
    std::vector<int> foundLoops;

    auto checkEquality = [](const ClipperPath& p1, const ClipperPath& p2) -> bool {
      auto accumulator = [](intType rhs, const ClipperPoint& point) -> intType { return rhs + point.x + point.y; };
      auto sumP1       = std::accumulate(p1.begin(), p1.end(), intType(0), accumulator);
      auto sumP2       = std::accumulate(p2.begin(), p2.end(), intType(0), accumulator);

      return intCmpEq(sumP1, sumP2, 100l);
    };

    for (auto it1 = clippedEdges.cbegin() + 1; it1 < clippedEdges.cend(); ++it1)
      for (auto it2 = clip_.cbegin() + 1; it2 < clip_.cend(); ++it2)
        if (checkEquality(*it1, *it2))
          foundLoops.push_back(std::distance(clip_.begin(), it2));

    if (not foundLoops.empty())
      innerLoops.emplace(foundLoops);
  }

  void clipElement(const auto& element) {
    auto index        = indexSet_.index(element);
    auto corners      = getElementCorners(index);
    auto elementEdges = getElementEdgesFromElementCorners(corners);

    ClipperPaths clippedEdges;
    // Rule NonZero to ensure that inner boundaries are defined clockwise and outer boundaries counter-clockwise
    /*
    if (clip_.size() == 1) {
      Clipper2Lib::Rect<intType> elementRect{corners[0][0], corners[1][1], corners[1][0], corners[3][1]};
      clippedEdges = Clipper2Lib::RectClip(elementRect, clip_);
    } else
     */
    clippedEdges = Clipper2Lib::Intersect(ClipperPaths{elementEdges}, clip_, Clipper2Lib::FillRule::NonZero);

    // At the moment there is no hole, if there are more than 2 Paths, then there is a hole
    std::optional<InnerLoops> innerLoops = std::nullopt;
    if (clippedEdges.size() > 1) {
      findInnerLoops(innerLoops, clippedEdges);
    }

    // If the clippedEdges are empty, this means the element is outside the clip -> empty
    // And if the clip has only 4 corners and these corners are the same as the element -> full
    if (clippedEdges.empty())
      trimFlags_.push_back(ElementTrimFlag::empty);
    else {
      if (isFullElement(clippedEdges.front(), corners) and not innerLoops.has_value())
        trimFlags_.push_back(ElementTrimFlag::full);
      else
        trimFlags_.push_back(ElementTrimFlag::trimmed);
    }

    if (trimFlags_.back() != ElementTrimFlag::trimmed) {
      clipResults_.push_back(std::nullopt);
      return;
    }

    // Fill pointMap to store the intersection points in regard to their corresponding edges
    auto pointMap = std::make_unique<IntersectionPointMap>();

    // Look if there are two distinct regions // Not yet implemented
    int regions = clippedEdges.size() - (innerLoops.has_value() ? innerLoops.value().loopIndices.size() : 0);
    if (regions > 1)
      std::cerr << "More than one region, undefined behaviour, inspect output " << std::endl;

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

    // Sort IntersectionPoints counter-clockwise
    sortIntersectionPoints(*pointMap);

    clipResults_[index].emplace(std::move(pointMap), edgeCounter, nodeCounter, std::move(innerLoops));
  }

  std::optional<ElementBoundaries> constructElementBoundaries(int directIndex) {
    std::optional<ElementBoundaries> result{};
    try {
      std::vector<Boundary> elementBoundaries;
      auto& clipResult = clipResults_[directIndex].value();
      auto corners     = getElementCorners(directIndex);

      // Determine the node where we begin to trace
      const auto nodeCounter = clipResult.nodeCounter;
      auto it_nodes          = std::ranges::find_if(nodeCounter, [](auto x) { return x > 0; });

      if (it_nodes == nodeCounter.end())
        throw std::runtime_error(
            "No Curve tracing scheme implemented for elements that are only trimmed on single element edges.");

      auto nodeToBegin = std::distance(nodeCounter.begin(), it_nodes);
      auto state       = FindNextBoundaryLoopState(nodeToBegin, corners, clipResult);

      // Iteratively find all the elementBoundaries
      bool found;
      do {
        found = findNextBoundary(&state, elementBoundaries);
      } while (!found);

      // Success
      result.emplace(elementBoundaries, extractInnerBoundaries(clipResult));
    } catch (const std::runtime_error& e) {
      std::cout << "Trimming encountered problem: " << e.what() << "\n Element is ignored" << std::endl;
    }
    return result;
  }

  std::optional<std::vector<std::vector<Boundary>>> extractInnerBoundaries(const ClippingResult& clipResult) const {
    if (not clipResult.innerLoops.has_value())
      return std::nullopt;

    std::vector<std::vector<Boundary>> innerBoundaries;
    innerBoundaries.resize(clipResult.innerLoops.value().loopIndices.size());

    for (auto c = 0u; auto& loopIndex : clipResult.innerLoops.value().loopIndices) {
      assert(loopIndex < trimData_->boundaryLoops.size());
      innerBoundaries[c].reserve(trimData_->boundaryLoops[loopIndex].size());
      std::ranges::copy(trimData_->boundaryLoops[loopIndex], std::back_inserter(innerBoundaries[c]));
      ++c;
    }

    return std::make_optional(innerBoundaries);
  }

  // @todo Get rid of weird and redundant naming of outer-loop variables
  bool findNextBoundary(FindNextBoundaryLoopState* state, std::vector<Boundary>& elementBoundaries) {
    int node       = state->currentNode;
    int nextEntity = node;

    // Check if nextEntity has no intersections and next node has one
    if (state->isCurrentlyOnNode && !state->hasIntersectionPointOnEdgeNr(nextEntity) &&
        state->hasIntersectionPointOnNodeNr(nextEntityIdx(node, 1))) {
      elementBoundaries.push_back(boundaryForEdge(state->corners, nextEntity));

      // Update Node and Point
      node                    = nextEntityIdx(node, 1);
      state->currentNode      = node;
      state->currentNodePoint = toFloatDomain(state->corners[node]);

      return (node == state->nodeToBegin);
    }
    if (state->isCurrentlyOnNode) {
      // Determine the next nextEntity with an IP
      for (int i = 0; i < 4; ++i) {
        auto nextEdge = nextEntityIdx(node, i);
        auto nextNode = nextEntityIdx(node, i + 1);

        if (state->hasIntersectionPointOnEdgeNr(nextEdge)) {
          nextEntity = nextEdge;
          break;
        }
        if (state->hasIntersectionPointOnNodeNr(nextNode)) {
          nextEntity = nextNode;
          break;
        }
      }
      // @todo Assert if nextEntity was found or not

      // Set node for the next step
      state->currentNode      = node;
      state->currentNodePoint = toFloatDomain(state->corners[node]);

      int boundaryIndexThatHasIntersectionPoint;
      IntersectionResult edgeIntersectResult;

      if (node == nextEntity) {
        /// Find the intersection Point on the current Edge
        auto findNextBoundaryResult =
            findBoundaryThatHasIntersectionPoint(state->pointMapPtr.get(), nextEntity, state->corners);

        if (not findNextBoundaryResult.has_value())
          throw std::runtime_error("findBoundaryThatHasIntersectionPoint not successful");

        boundaryIndexThatHasIntersectionPoint = findNextBoundaryResult.value().edge;
        edgeIntersectResult                   = findNextBoundaryResult.value().result;

        elementBoundaries.emplace_back(state->currentNodePoint, edgeIntersectResult.globalResult);
      } else {
        auto findNextBoundaryResult = findNextBoundaryThatHasIntersectionPointOnNode(state, node);

        if (not findNextBoundaryResult.has_value())
          throw std::runtime_error("findNextBoundaryThatHasIntersectionPointOnNode not successful");

        boundaryIndexThatHasIntersectionPoint = findNextBoundaryResult.value().edge;
        edgeIntersectResult                   = findNextBoundaryResult.value().result;
      }

      /// Curve Tracing
      auto traceResult =
          traceCurve(state, {boundaryIndexThatHasIntersectionPoint, edgeIntersectResult.localResult, nextEntity});

      if (not traceResult.has_value())
        throw std::runtime_error("traceCurve not successful");

      for (const auto& result : traceResult.value()) {
        auto curveToTrace    = globalBoundaries_[result.boundaryIdx].nurbsGeometry;
        Boundary newBoundary = Boundary(curveToTrace, Utilities::Domain<double>{result.beginU, result.intersectU});
        elementBoundaries.push_back(newBoundary);
      }

      auto lastResult          = traceResult.value().back();
      state->isCurrentlyOnNode = lastResult.onNode;
      state->currentNodePoint  = lastResult.intersectionPoint;

      if (lastResult.onNode)
        state->currentNode = lastResult.edge;
      else {
        state->edgeTheLoopIsOn = lastResult.edge;
        state->moreIntersectionPointsOnEdge =
            getNextIntersectionPointOnEdge(state->pointMapPtr.get(), state->edgeTheLoopIsOn).second;
      }
      return lastResult.onNode && state->currentNode == state->nodeToBegin;
    }

    if (!(state->isCurrentlyOnNode) && state->moreIntersectionPointsOnEdge) {
      nextEntity = state->edgeTheLoopIsOn;

      /// Find the intersection Point on the current Edge

      auto findNextBoundaryResult =
          findBoundaryThatHasIntersectionPoint(state->pointMapPtr.get(), nextEntity, state->corners);

      if (not findNextBoundaryResult.has_value())
        throw std::runtime_error("findBoundaryThatHasIntersectionPoint not successful");

      auto [boundaryIndexThatHasIntersectionPoint, edgeIntersectResult] = findNextBoundaryResult.value();
      elementBoundaries.emplace_back(state->currentNodePoint, edgeIntersectResult.globalResult);

      /// Curve Tracing

      auto traceResult =
          traceCurve(state, {boundaryIndexThatHasIntersectionPoint, edgeIntersectResult.localResult, nextEntity});

      if (not traceResult.has_value())
        throw std::runtime_error("tracCurve not successful");

      for (const auto& result : traceResult.value()) {
        auto curveToTrace    = globalBoundaries_[result.boundaryIdx].nurbsGeometry;
        Boundary newBoundary = Boundary(curveToTrace, Utilities::Domain<double>{result.beginU, result.intersectU});
        elementBoundaries.push_back(newBoundary);
      }

      auto lastResult          = traceResult.value().back();
      state->isCurrentlyOnNode = lastResult.onNode;
      state->currentNodePoint  = lastResult.intersectionPoint;

      if (lastResult.onNode)
        state->currentNode = lastResult.edge;
      else {
        state->edgeTheLoopIsOn = lastResult.edge;
        state->moreIntersectionPointsOnEdge =
            getNextIntersectionPointOnEdge(state->pointMapPtr.get(), state->edgeTheLoopIsOn).second;
      }

      return lastResult.onNode && state->currentNode == state->nodeToBegin;
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

      return node == state->nodeToBegin;
    }
    __builtin_unreachable();
  }

  ////////////
  /// Non-static Helper functions
  ////////////

  /// Returns the element corners of the tensorGrid in the IntegerDomain
  std::vector<FieldVector<intType, dim>> getElementCorners(int index) {
    auto geo = elementGeometries_[index];

    // The corners are sorted counter-clockwise to get a closed loop, also see top of this file for numbering
    return {geo.corner(0), geo.corner(1), geo.corner(3), geo.corner(2)};
  }
  Clipper2Lib::Path<intType> getElementEdges(int index) {
    auto corners = getElementCorners(index);
    return getElementEdgesFromElementCorners(corners);
  }

  std::vector<FieldVector<double, dim>> getDuneElementCorners(int index) {
    std::vector<FieldVector<double, 2>> dune_corners;
    std::ranges::transform(getElementCorners(index), std::back_inserter(dune_corners),
                           [&](const auto& c) { return toFloatDomain(c); });
    std::swap(dune_corners[2], dune_corners[3]);
    return dune_corners;
  }

  std::optional<FindNextBoundaryResult> findBoundaryThatHasIntersectionPoint(IntersectionPointMap* pointMapPtr,
                                                                             int edge, std::vector<IntPoint>& corners) {
    auto edgeCurve               = boundaryForEdge(corners, edge);
    Point intersectionPointGuess = getNextIntersectionPointOnEdge(pointMapPtr, edge).first;

    for (int i = 0; auto& boundary : globalBoundaries_) {
      if (areLinesParallel(boundary.nurbsGeometry, edgeCurve.nurbsGeometry)) {
        ++i;
        continue;
      }

      auto edgeIntersectResult = findclosestPointOnCurve(boundary.nurbsGeometry, intersectionPointGuess);

      if (edgeIntersectResult.found) {
        markIntersectionPointAsVisited(pointMapPtr, edge, intersectionPointGuess);
        assert(Dune::FloatCmp::eq(edgeIntersectResult.globalResult, intersectionPointGuess, 1e-3));
        return std::make_optional(FindNextBoundaryResult(i, edgeIntersectResult));
      }
      ++i;
    }

    return std::nullopt;
  }
  std::optional<FindNextBoundaryResult> findNextBoundaryThatHasIntersectionPointOnNode(FindNextBoundaryLoopState* state,
                                                                                       int node) {
    auto nodeCounter = state->nodeCounter;
    if (nodeCounter[node] != 1)
      return std::nullopt;

    Point intersectionPointGuess = toFloatDomain(state->corners[node]);
    for (int i = 0; auto& boundary : globalBoundaries_) {
      auto result = findclosestPointOnCurve(boundary.nurbsGeometry, intersectionPointGuess);

      if (result.found)
        return FindNextBoundaryResult(i, result);
      ++i;
    }
    return std::nullopt;
  }
  // @todo: Get rid of initializer lists
  std::optional<std::vector<TraceCurveOutput>> traceCurve(FindNextBoundaryLoopState* state,
                                                          TraceCurveInput traceCurveInput) {
    CurveGeometry curveToTrace  = globalBoundaries_[traceCurveInput.boundaryIdx].nurbsGeometry;
    auto [startPoint, endPoint] = curveStartEndPoint(curveToTrace);

    // There is an edge case where we have a closed curve (e.g. a circle) that has endPoint and startPoint at the same
    // point, and that point is in the element
    if (Dune::FloatCmp::eq(startPoint, endPoint) && pointInElement(startPoint, state->corners)) {
      auto traceCurveResult = traceCurveImpl(state, curveToTrace, traceCurveInput.startEdge);
      if (traceCurveResult.has_value())
        return std::make_optional<std::vector<TraceCurveOutput>>({
            {traceCurveInput.boundaryIdx,          traceCurveInput.startU, curveToTrace.domain()[0].right(),
             traceCurveResult->intersectionPoint,                     -1, false},
            {traceCurveInput.boundaryIdx, curveToTrace.domain()[0].left(),     traceCurveResult->intersectU,
             traceCurveResult->intersectionPoint, traceCurveResult->edge, false}
        });
      else
        return std::nullopt;
    }
    // ... or the point is on an edge, even nastier
    if (Dune::FloatCmp::eq(startPoint, endPoint) && pointInElementOrOnEdge(startPoint, state->corners)) {
      if (Dune::FloatCmp::ne(curveToTrace(traceCurveInput.startU), startPoint)) {
        auto traceCurveResult = traceCurveImpl(state, curveToTrace, traceCurveInput.startEdge);
        if (traceCurveResult.has_value())
          return std::make_optional<std::vector<TraceCurveOutput>>({
              {traceCurveInput.boundaryIdx, traceCurveInput.startU, curveToTrace.domain()[0].right(),
               traceCurveResult->intersectionPoint, traceCurveResult->edge, false}
          });
      }
    }

    // Check the boundary that was given in traceInput
    auto traceCurveResult = traceCurveImpl(state, curveToTrace, traceCurveInput.startEdge);

    if (traceCurveResult.has_value())
      return std::make_optional<std::vector<TraceCurveOutput>>({
          {traceCurveInput.boundaryIdx, traceCurveInput.startU, traceCurveResult->intersectU,
           traceCurveResult->intersectionPoint, traceCurveResult->edge, traceCurveResult.value().onNode}
      });

    // If this is unsuccessful the curve is split into two parametrisation

    int status = -1;
    // Now check if the startPoint or endPoint of the curve is in the element
    if (pointInElement(endPoint, state->corners))
      status = 1;
    else if (pointInElement(startPoint, state->corners))
      status = 0;

    // If neither start nor endpoint is in the element then traceCurveImpl most likely made a mistake
    if (status == -1)
      return std::nullopt;

    // Now we check if there is another boundary that begins at this point and trace this potential curve
    Point checkPoint      = (status == 1) ? endPoint : startPoint;
    double endUOfCurve1   = (status == 1) ? curveToTrace.domain()[0].right() : curveToTrace.domain()[0].left();
    double startUOfCurve2 = std::numeric_limits<double>::quiet_NaN();

    // Find the next BoundarySegment
    auto it = std::ranges::find_if(globalBoundaries_, [&](auto& _boundary) {
      auto [startPointOfBoundaryToCheck, endPointOfBoundaryToCheck] = curveStartEndPoint(_boundary.nurbsGeometry);
      // Check if it's the same boundary Geometry
      if (startPointOfBoundaryToCheck == startPoint && endPointOfBoundaryToCheck == endPoint)
        return false;

      if (Dune::FloatCmp::eq(checkPoint, startPointOfBoundaryToCheck, 1e-8)) {
        startUOfCurve2 = _boundary.domain.left();
        return true;
      } else if (Dune::FloatCmp::eq(checkPoint, endPointOfBoundaryToCheck, 1e-8)) {
        startUOfCurve2 = _boundary.domain.right();
        return true;
      }
      return false;
    });

    if (it != globalBoundaries_.end()) {
      int newBoundaryIdx     = static_cast<int>(std::ranges::distance(globalBoundaries_.begin(), it));
      auto curveToTrace2     = globalBoundaries_[newBoundaryIdx].nurbsGeometry;
      auto traceCurveResult2 = traceCurveImpl(state, curveToTrace2, traceCurveInput.startEdge);

      assert(Dune::FloatCmp::eq(curveToTrace(endUOfCurve1), curveToTrace2(startUOfCurve2), 1e-8));

      if (traceCurveResult2.has_value())
        return std::make_optional<std::vector<TraceCurveOutput>>({
            {traceCurveInput.boundaryIdx, traceCurveInput.startU,                  endUOfCurve1,                           checkPoint,-1, false                                 },
            {             newBoundaryIdx,         startUOfCurve2, traceCurveResult2->intersectU, traceCurveResult2->intersectionPoint,
             traceCurveResult2->edge, false}
        });
    }

    return std::nullopt;
  }

  //////////////////
  /// Static Helpers
  /////////////////

  static Clipper2Lib::Paths<intType> getClip(TrimData* trimData) {
    Clipper2Lib::Paths<intType> clipPaths;
    Clipper2Lib::Path<intType> tempPath;
    for (auto& loop : trimData->boundaryLoops) {
      tempPath.clear();
      for (auto& boundary : loop)
        for (auto& point : boundary.path<intType>(pathSamples, scale))
          tempPath.push_back(point);

      // Sanitize Path, so we don't get duplicate intersection Points on corners
      tempPath = Clipper2Lib::TrimCollinear(tempPath);
      clipPaths.push_back(tempPath);
    }

    return clipPaths;
  }

  static intType toIntDomain(double x) {
    return x * scaleFactor();
  }
  static FieldVector<intType, 2> toIntDomain(const FieldVector<double, 2>& x) {
    return x * scaleFactor();
  }

  static double toFloatDomain(intType x) {
    return x / scaleFactor();
  }
  static FieldVector<double, 2> toFloatDomain(const FieldVector<intType, 2>& x) {
    return x / scaleFactor();
  }

  static std::vector<Boundary> extractBoundaries(TrimData* trimData) {
    std::vector<Boundary> boundaries;
    boundaries.reserve(trimData->numBoundaries());

    for (auto& loop : trimData->boundaryLoops) {
      std::ranges::copy(loop, std::back_inserter(boundaries));
    }
    return boundaries;
  }

  static inline double distance(const ClipperPoint p1, const ClipperPoint p2) {
    return std::hypot(p1.x - p2.x, p1.y - p2.y);
  }

  static inline double distance(const Point p1, const Point p2) {
    return std::hypot(p1[0] - p2[0], p1[1] - p2[1]);
  }

  static inline bool samePoint(const IntPoint& a, const ClipperPoint& b) {
    return a == IntPoint{b.x, b.y};
  }

  static inline int getEdgeOrientation(const int edge) {
    return (edge == 0 || edge == 2) ? 0 : 1;
  }

  static inline int nextEntityIdx(const int i, const int x) {
    return (i + x) % 4;
  }

  template <std::integral T>
  static inline bool intCmpEq(T i1, T i2, T tol = 0) {
    return std::abs(i1 - i2) <= tol;
  }

  static std::array<Point, 2> curveStartEndPoint(const CurveGeometry& curve) {
    std::array<int, 2> indices{0, static_cast<int>(curve.patchData_.controlPoints.strideSizes()[0]) - 1};
    return {curve.patchData_.controlPoints.get({indices[0]}).p, curve.patchData_.controlPoints.get({indices[1]}).p};
  }

  /// This scales back to the floatDomain
  static Boundary boundaryForEdge(const std::vector<IntPoint>& corners, int edge) {
    return Boundary(toFloatDomain(corners[nextEntityIdx(edge, 0)]), toFloatDomain(corners[nextEntityIdx(edge, 1)]));
  }

  static std::pair<Point, bool> getNextIntersectionPointOnEdge(IntersectionPointMap* pointMapPtr, int edge) {
    if ((*pointMapPtr)[edge].empty())
      return {Point(), false};

    for (const IntersectionPoint& ip : (*pointMapPtr)[edge]) {
      if (!ip.alreadyVisited)
        return {ip.point, true};
    }
    return {Point(), false};
  }

  static void markIntersectionPointAsVisited(IntersectionPointMap* pointMapPtr, int edge, const Point& point) {
    auto it = std::find_if((*pointMapPtr)[edge].begin(), (*pointMapPtr)[edge].end(),
                           [&point](IntersectionPoint x) { return (Dune::FloatCmp::eq(x.point, point)); });
    assert(it != (*pointMapPtr)[edge].end());
    (*it).alreadyVisited = true;
  }

  static double findGoodStartingPoint(CurveGeometry& curve, Point& intersectionPoint, int N) {
    auto linSpace = Utilities::linspace(curve.domain()[0], N);
    std::vector<double> distances;
    std::ranges::transform(linSpace, std::back_inserter(distances),
                           [&](const auto& u) { return distance(curve(u), intersectionPoint); });
    auto min_it = std::ranges::min_element(distances);

    auto min_idx = std::ranges::distance(distances.begin(), min_it);
    return linSpace[min_idx];
  }

  static IntersectionResult findclosestPointOnCurve(CurveGeometry& curve, Point& intersectionPoint) {
    // There are some easy cases where the IP is at the back or at the front of the domain, check these cases first
    auto [start, end] = curveStartEndPoint(curve);

    const auto domain = curve.domain()[0];

    if (Dune::FloatCmp::eq(start, intersectionPoint, tolerance))
      return {domain.left(), start, true};
    if (Dune::FloatCmp::eq(end, intersectionPoint, tolerance))
      return {domain.right(), end, true};
    Dune::FieldVector<double, 1> uGuess{findGoodStartingPoint(curve, intersectionPoint, 10)};
    auto [u, Ru, fu, gap] = Dune::IGA::closestPointProjectionByTrustRegion(curve, intersectionPoint, uGuess);

    if (gap > gapTolerance)
      return {u[0], intersectionPoint, false};
    return {u[0], curve(u[0]), true};
  }

  static std::optional<TraceCurveResult> traceCurveImpl(FindNextBoundaryLoopState* state, CurveGeometry& curveToTrace,
                                                        const int startEdge) {
    // Look through the edges and determine the next Intersection Point (the IPs are sorted counter-clockwise)
    auto pointMapPtr = state->pointMapPtr.get();
    for (int i = 1; i < 5; ++i) {
      auto entityIdx                                 = nextEntityIdx(startEdge, i);
      auto [intersectionPoint, hasIntersectionPoint] = getNextIntersectionPointOnEdge(pointMapPtr, entityIdx);
      if (hasIntersectionPoint) {
        // Check if the curveToTrace goes through the intersection point
        auto [localResult, globalResult, found] = findclosestPointOnCurve(curveToTrace, intersectionPoint);
        if (found) {
          markIntersectionPointAsVisited(pointMapPtr, entityIdx, intersectionPoint);
          return std::make_optional<TraceCurveResult>({localResult, globalResult, entityIdx, false});
        }
      }
      // It may also be on the corresponding Node
      if (state->nodeCounter[entityIdx] > 0 and entityIdx != state->currentNode) {
        auto nodePoint                          = toFloatDomain(state->corners[entityIdx]);
        auto [localResult, globalResult, found] = findclosestPointOnCurve(curveToTrace, nodePoint);
        if (found)
          return std::make_optional<TraceCurveResult>({localResult, globalResult, entityIdx, true});
      }
    }

    return std::nullopt;
  }

  static bool areLinesParallel(const CurveGeometry& line1, const CurveGeometry& line2) {
    if (line1.degree()[0] > 1 || line2.degree()[0] > 1)
      return false;

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

  static bool pointInElement(const Point& point, const std::vector<IntPoint>& corners) {
    auto lowerLeft  = toFloatDomain(corners[0]);
    auto upperRight = toFloatDomain(corners[2]);

    return FloatCmp::gt(point[0], lowerLeft[0]) and FloatCmp::lt(point[0], upperRight[0]) and
           FloatCmp::gt(point[1], lowerLeft[1]) and FloatCmp::lt(point[1], upperRight[1]);
  }
  static bool pointInElementOrOnEdge(const Point& point, const std::vector<IntPoint>& corners) {
    auto lowerLeft  = toFloatDomain(corners[0]);
    auto upperRight = toFloatDomain(corners[2]);

    return FloatCmp::ge(point[0], lowerLeft[0], tolerance) and FloatCmp::le(point[0], upperRight[0], tolerance) and
           FloatCmp::ge(point[1], lowerLeft[1], tolerance) and FloatCmp::le(point[1], upperRight[1], tolerance);
  }

  static bool isFullElement(ClipperPath& clippedEdges, std::vector<IntPoint>& corners) {
    if (clippedEdges.size() != 4)
      return false;

    return std::ranges::all_of(clippedEdges, [&corners](const ClipperPoint& point) {
      auto it = std::ranges::find_if(corners, [&point](const IntPoint& corner) { return samePoint(corner, point); });
      return it != corners.end();
    });
  }

  /// @brief Checks if a point p is on the line spanned by the points a and b within a given tolerance
  static bool pointOnLine(const ClipperPoint& p, const ClipperPoint& a, const ClipperPoint& b,
                          const double tol = tolerance) {
    return distance(a, p) + distance(b, p) == distance(a, b);
  }

  static std::pair<bool, int> pointOnAnyEdge(const ClipperPoint& point, const auto& elementEdges) {
    for (int i = 0; i < 4; ++i)
      if (pointOnLine(point, elementEdges[i], elementEdges[nextEntityIdx(i, 1)]))
        return std::make_pair(true, i);

    return std::make_pair(false, -1);
  }

  static std::pair<bool, int> pointOnAnyVertex(const ClipperPoint& point, const std::vector<IntPoint>& corners) {
    for (size_t i = 0; i < corners.size(); ++i)
      if (samePoint(corners[i], point))
        return std::make_pair(true, static_cast<int>(i));

    return std::make_pair(false, -1);
  }

  static Clipper2Lib::Path<intType> getElementEdgesFromElementCorners(const std::vector<IntPoint>& corners) {
    Clipper2Lib::Path<intType> edges;
    for (int i = 0; i < 4; ++i)
      edges.emplace_back(corners[i][0], corners[i][1]);
    return edges;
  }

  /////////////////
  /// Definition of Result Types
  /////////////////

  struct TraceCurveInput
  {
    int boundaryIdx;
    double startU;
    int startEdge;
  };
  struct TraceCurveOutput
  {
    int boundaryIdx{};
    double beginU{};
    double intersectU{};
    Point intersectionPoint;
    int edge{};
    bool onNode{};
  };
  struct TraceCurveResult
  {
    double intersectU{};
    Point intersectionPoint;
    int edge{};
    bool onNode = false;
  };
  struct IntersectionPoint
  {
    Point point;
    bool alreadyVisited = false;

    explicit IntersectionPoint(const ClipperPoint& _clipperPoint)
        : point({toFloatDomain(_clipperPoint.x), toFloatDomain(_clipperPoint.y)}){};
  };

  struct IntersectionResult
  {
    double localResult{};
    Point globalResult;
    bool found = false;
  };

  struct FindNextBoundaryResult
  {
    int edge{};
    IntersectionResult result;
  };

  struct InnerLoops
  {
    std::vector<int> loopIndices;
  };

  struct ClippingResult
  {
    std::shared_ptr<IntersectionPointMap> pointMapPtr;
    std::array<int, 4> edgeCounter{};
    std::array<int, 4> nodeCounter{};
    std::optional<InnerLoops> innerLoops{};
  };

  struct FindNextBoundaryLoopState
  {
    // Variables that need to be set
    const int nodeToBegin;
    std::vector<IntPoint> corners;
    std::shared_ptr<IntersectionPointMap> pointMapPtr;
    const std::array<int, 4> edgeCounter;
    const std::array<int, 4> nodeCounter;

    // Loop States
    bool isCurrentlyOnNode            = true;
    bool moreIntersectionPointsOnEdge = false;
    int edgeTheLoopIsOn               = -1; // This variable is only filed when !loopIsCurrentlyOnNode

    // Current positions (are initially set through the constructor)
    int currentNode;
    Point currentNodePoint;

    // Helper Functions
    bool hasIntersectionPointOnEdgeNr(int nr) {
      return edgeCounter[nr] > 0;
    };
    bool hasIntersectionPointOnNodeNr(int nr) {
      return nodeCounter[nr] > 0;
    };

    FindNextBoundaryLoopState(int _nodeToBegin, std::vector<IntPoint>& _corners, ClippingResult& _result)
        : nodeToBegin(_nodeToBegin),
          corners(_corners),
          pointMapPtr(std::move(_result.pointMapPtr)),
          edgeCounter(_result.edgeCounter),
          nodeCounter(_result.nodeCounter) {
      currentNode      = nodeToBegin;
      currentNodePoint = corners[currentNode];
    };

    FindNextBoundaryLoopState()                                 = delete;
    FindNextBoundaryLoopState(const FindNextBoundaryLoopState&) = delete;
  };
};
} // namespace Dune::IGA::Trim

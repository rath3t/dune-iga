// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace Dune::IGA::DefaultParameterSpace {

template <typename ParameterSpace>
void LoopHandler<ParameterSpace>::processHostToHost(const Vertex& vV1, const Vertex& vV2) {
  const auto pt1 = vV1.pt;
  const auto pt2 = vV2.pt;

  if (!Util::isConsecutive(vV1.hostId(), vV2.hostId()))
    processNonConsecutiveHosts(vV1, vV2, pt1, pt2);
  else
    processConsecutiveHosts(vV1, vV2, pt2);
}

template <typename ParameterSpace>
auto LoopHandler<ParameterSpace>::findIntersection(const Vertex& vertex, const Clipper2Lib::PointD& pt) {
  for (auto&& val : getTrimmingCurveIdxFromHost(vertex)) {
    currentCurveIdx_ = val;
    auto result = Util::callFindIntersection(patchTrimData_.getCurve(currentCurveIdx_), vertex.hostId(), pt, corners_);
    if (result.has_value() && Util::approxSamePoint(pt, result.value().curvePoint, 0.01))
      return result.value();
  }
  throwGridError();
  __builtin_unreachable();
}

template <typename ParameterSpace>
void LoopHandler<ParameterSpace>::processNonConsecutiveHosts(const Vertex& vV1, const Vertex& vV2, const Point& pt1,
                                                             const Point& pt2) {
  // Process first vertex if needed
  if (foundVertices_.empty()) {
    auto [tParam, curvePoint] = findIntersection(vV1, pt1);
    foundVertices_.push_back(curvePoint);
  }

  // Process second vertex
  auto [tParam2, curvePoint2] = findIntersection(vV2, pt2);
  auto elementTrimmingCurve =
      Util::createTrimmingCurveSlice(patchTrimData_.getCurve(currentCurveIdx_), currentT_, tParam2);
  elementTrimData_.addEdgeNewNew(elementTrimmingCurve, curvePoint2, -1);
  addNewBoundarySegementIdx(currentCurveIdx_, currentT_, tParam2);
}

template <typename ParameterSpace>
void LoopHandler<ParameterSpace>::processConsecutiveHosts(const Vertex& vV1, const Vertex& vV2, const Point& pt2) {
  elementTrimData_.addEdge(Util::giveEdgeIdx(vV1.hostId(), vV2.hostId()));
  checkAndAddBoundarySegmentIdx(Util::giveEdgeIdx(vV1.hostId(), vV2.hostId()));
  foundVertices_.push_back({pt2.x, pt2.y});
}

template <typename ParameterSpace>
void LoopHandler<ParameterSpace>::processHostToNew(const Vertex& vV1, const Vertex& vV2) {
  auto pt1 = vV1.pt;
  auto pt2 = vV2.pt;

  if (vV1.hostId() != vV2.edgeId()) {
    if (foundVertices_.empty()) {
      currentCurveIdx_ = getTrimmingCurveIdxFromHost(vV1)[0];
      auto [currentT, curvePoint] =
          Util::callFindIntersection(patchTrimData_.getCurve(currentCurveIdx_), vV1.hostId(), pt1, corners_).value();
      assertPoint(pt1, curvePoint);

      foundVertices_.push_back(curvePoint);
    }
    auto [tParam, curvePoint] =
        Util::callFindIntersection(patchTrimData_.getCurve(currentCurveIdx_), vV2.edgeId(), pt2, corners_).value();
    assertPoint(pt2, curvePoint);

    auto elementTrimmingCurve =
        Util::createTrimmingCurveSlice(patchTrimData_.getCurve(currentCurveIdx_), currentT_, tParam);
    elementTrimData_.addEdgeNewNew(elementTrimmingCurve, curvePoint, vV2.edgeId());
    addNewBoundarySegementIdx(currentCurveIdx_, currentT_, tParam);

    foundVertices_.push_back(curvePoint);

  } else {
    currentCurveIdx_ = getTrimmingCurveIdx(vV2);
    auto [tParam, curvePoint] =
        Util::callFindIntersection(patchTrimData_.getCurve(currentCurveIdx_), vV2.edgeId(), pt2, corners_).value();
    assertPoint(pt2, curvePoint);

    FieldVector<ScalarType, dim> p =
        foundVertices_.empty() ? FieldVector<ScalarType, dim>{pt1.x, pt1.y} : foundVertices_.back();
    auto trimmedEdge = Util::createHostGeometry<TrimmingCurve>(p, curvePoint);
    elementTrimData_.addEdgeHostNew(vV2.edgeId(), trimmedEdge, curvePoint);
    checkAndAddBoundarySegmentIdx(vV2.edgeId());
    foundVertices_.push_back(curvePoint);

    currentT_ = tParam;
  }
}

template <typename ParameterSpace>
void LoopHandler<ParameterSpace>::processNewToNew(const Vertex& vV1, const Vertex& vV2) {
  auto pt1 = vV1.pt;
  auto pt2 = vV2.pt;

  // If there is no Vertex yet found add search for vV1 and add it to the foundVertices, also set currentCurveIdx
  if (foundVertices_.empty()) {
    currentCurveIdx_ = getTrimmingCurveIdx(vV1);
    auto [currentT, curvePoint] =
        Util::callFindIntersection(patchTrimData_.getCurve(currentCurveIdx_), vV1.edgeId(), pt1, corners_).value();
    assertPoint(pt1, curvePoint);

    foundVertices_.push_back(curvePoint);
  }
  if (getTrimmingCurveIdx(vV2).curve != currentCurveIdx_.curve and
      getTrimmingCurveIdx(vV2).loop != currentCurveIdx_.loop)
    throwGridError();

  auto [tParam, curvePoint] =
      Util::callFindIntersection(patchTrimData_.getCurve(currentCurveIdx_), vV2.edgeId(), pt2, corners_).value();
  assertPoint(pt2, curvePoint);

  // If currentT > tParam it can have 2 reasons. 1) the trim id only on one side of the element, and it has to
  // connect back or 2) if the trimming curve consist only of one curve where it has the same front and back
  // ControlPoint, so sometimes the wrong tParam (front or back) gets determined
  if (currentT_ > tParam) {
    bool success = false;
    if (currentCurveIdx_.loop > 0 && patchTrimData_.loops()[currentCurveIdx_.loop].size() == 1) {
      if (auto curve = patchTrimData_.getCurve(currentCurveIdx_); curve.isConnectedAtBoundary(0)) {
        auto elementTrimmingCurve = Util::createTrimmingCurveSlice(curve, currentT_, curve.domain()[0].back());
        elementTrimData_.addEdgeNewNew(elementTrimmingCurve, curvePoint, vV2.edgeId());
        addNewBoundarySegementIdx(currentCurveIdx_, currentT_, curve.domain()[0].back());

        success = true;
      }
    } else if (vV1.edgeId() == vV2.edgeId()) {
      auto trimmedEdge = Util::createHostGeometry<TrimmingCurve>(foundVertices_.back(), curvePoint);
      elementTrimData_.addEdgeNewNewOnHost(vV2.edgeId(), trimmedEdge, curvePoint);
      checkAndAddBoundarySegmentIdx(vV2.edgeId());
      checkAndAddBoundarySegmentIdx(vV2.edgeId());

      currentT_ = tParam;
      success   = true;
    }
    if (not success)
      throwGridError();
  } else {
    auto elementTrimmingCurve =
        Util::createTrimmingCurveSlice(patchTrimData_.getCurve(currentCurveIdx_), currentT_, tParam);
    elementTrimData_.addEdgeNewNew(elementTrimmingCurve, curvePoint, vV2.edgeId());
    addNewBoundarySegementIdx(currentCurveIdx_, currentT_, tParam);
  }

  foundVertices_.push_back(curvePoint);
}

template <typename ParameterSpace>
void LoopHandler<ParameterSpace>::processNewToHost(const Vertex& vV1, const Vertex& vV2) {
  auto pt1 = vV1.pt;
  auto pt2 = vV2.pt;

  if (not Util::isConsecutive(vV1.edgeId(), vV2.hostId())) {
    auto [tParam, curvePoint] =
        Util::callFindIntersection(patchTrimData_.getCurve(currentCurveIdx_), vV2.hostId(), pt2, corners_).value();
    assertPoint(pt2, curvePoint);

    auto elementTrimmingCurve =
        Util::createTrimmingCurveSlice(patchTrimData_.getCurve(currentCurveIdx_), currentT_, tParam);
    elementTrimData_.addEdgeNewNew(elementTrimmingCurve, curvePoint, -1);
    addNewBoundarySegementIdx(currentCurveIdx_, currentT_, tParam);

  } else {
    auto v2 = Dune::FieldVector<ScalarType, dim>{pt2.x, pt2.y};

    FieldVector<ScalarType, dim> p;
    if (foundVertices_.empty()) {
      auto findIntersectionResult =
          Util::callFindIntersection(patchTrimData_.getCurve(getTrimmingCurveIdx(vV1)), vV1.edgeId(), pt1, corners_)
              .value();
      p = findIntersectionResult.curvePoint;
    } else
      p = foundVertices_.back();

    auto trimmedEdge = Util::createHostGeometry<TrimmingCurve>(p, v2);
    elementTrimData_.addEdgeNewHost(vV1.edgeId(), trimmedEdge, vV2.hostId());
    checkAndAddBoundarySegmentIdx(vV1.edgeId());

    foundVertices_.push_back(v2);
  }
}

template <typename ParameterSpace>
bool LoopHandler<ParameterSpace>::processAdditionalCases(const Vertex& vV1, const Vertex& vV2) {
  auto pt1 = vV1.pt;
  auto pt2 = vV2.pt;

  if (vV1.isNew() and vV2.isInside()) {
    if (foundVertices_.empty()) {
      currentCurveIdx_ = getTrimmingCurveIdx(vV1);
      FieldVector<ScalarType, dim> curvePoint;
      auto findIntersectionResult =
          Util::callFindIntersection(patchTrimData_.getCurve(currentCurveIdx_), vV1.edgeId(), pt1, corners_).value();
      currentT_  = findIntersectionResult.tParam;
      curvePoint = findIntersectionResult.curvePoint;
      assertPoint(pt1, curvePoint);

      foundVertices_.push_back(curvePoint);
    }

    const auto& curve         = patchTrimData_.loops()[vV2.loop()].curves()[vV2.formerCurve()];
    double tParam             = curve.domain().front().back();
    auto elementTrimmingCurve = Util::createTrimmingCurveSlice(curve, currentT_, tParam);
    auto curvePoint           = curve.corner(1);
    elementTrimData_.addEdgeNewNew(elementTrimmingCurve, curvePoint, vV1.edgeId());
    addNewBoundarySegementIdx(
        CurveIndex{.loop = static_cast<int>(vV2.loop()), .curve = static_cast<int>(vV2.formerCurve())}, currentT_,
        tParam);
    foundVertices_.push_back(curvePoint);
  } else if (vV1.isInside() and vV2.isNew()) {
    const auto& curve         = patchTrimData_.loops()[vV1.loop()].curves()[vV1.subsequentCurve()];
    currentT_                 = curve.domain().front().front();
    auto [tParam, curvePoint] = Util::callFindIntersection(curve, vV2.edgeId(), pt2, corners_).value();
    assertPoint(pt2, curvePoint);

    auto elementTrimmingCurve = Util::createTrimmingCurveSlice(curve, currentT_, tParam);
    elementTrimData_.addEdgeNewNew(elementTrimmingCurve, curvePoint, vV2.edgeId());
    addNewBoundarySegementIdx(
        CurveIndex{.loop = static_cast<int>(vV1.loop()), .curve = static_cast<int>(vV1.subsequentCurve())}, currentT_,
        tParam);

    foundVertices_.push_back(curvePoint);
  } else {
    return false;
  }
  return true;
}

} // namespace Dune::IGA::DefaultParameterSpace
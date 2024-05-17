// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <clipper2/clipper.h>

#include <dune/iga/geometrykernel/findintersection.hh>
#include <dune/iga/geometrykernel/slicecurve.hh>
#include <dune/iga/trimmer/defaulttrimmer/elementtrimdata.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/clipelementrectangle.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/trimutils.hh>

namespace Dune::IGA::DefaultTrim {

template <int dim, int dimworld, typename ScalarType>
auto TrimmerImpl<dim, dimworld, ScalarType>::trimElement(const YASPEntity<0>& element, const auto& gv,
                                                         const PatchTrimData& patchTrimData,
                                                         bool initFlag) -> ElementTrimData {
  auto geo = element.geometry();

  using CurveIndex = Impl::CurveLoopIndexEncoder::IndexResult;

  static constexpr int numberOfCorners = 4;
  std::array<FieldVector<double, 2>, numberOfCorners> corners;
  for (const auto i : Dune::range(numberOfCorners))
    corners[i] = geo.corner(Util::vertexIndexMapping[i]);

  std::array<std::pair<bool, size_t>, numberOfCorners> boundarySegementIndices;
  for (const auto& intersection : Dune::intersections(gv, element)) {
    auto idx = Transformations::mapToTrimmer(1, intersection.indexInInside());
    if (not intersection.boundary())
      boundarySegementIndices[idx] = std::make_pair(false, std::numeric_limits<size_t>::max());
    else
      boundarySegementIndices[idx] = std::make_pair(true, intersection.boundarySegmentIndex());
  }

  auto [flag, result] = Util::clipElementRectangle(element, patchTrimData);

  ElementTrimData elementTrimData(flag, element);
  if (flag != ElementTrimFlag::trimmed)
    return elementTrimData;

  auto nextEntity                  = [&](const int i) { return (i + 1) % result.vertices_.size(); };
  auto getTrimmingCurveIdx         = [&](auto& vV) -> auto { return patchTrimData.getIndices(vV.zValue()); };
  auto getTrimmingCurveIdxFromHost = [&](auto& vV) -> auto {
    return patchTrimData.getIndices(vV.additionalZValue().value());
  };

  auto throwGridError = [] {
    DUNE_THROW(Dune::GridError, "TrimElement wasn't successfull. Could not connect element trim curves");
  };

  auto assertPoint = [&](const Clipper2Lib::PointD& ip, const FieldVector<ScalarType, dim>& curvePt,
                         double prec = 0.01) {
    if (not(FloatCmp::eq(ip.x, curvePt[0], prec) and FloatCmp::eq(ip.y, curvePt[1], prec))) {
      std::cout << "Found ClipperPoint " << ip << " does not coincide with determined point on Curve " << curvePt
                << std::endl;
      throwGridError();
    }
  };

  auto checkAndAddBoundarySegmentIdx = [&](int edgeIdx) {
    if (boundarySegementIndices[edgeIdx].first)
      elementTrimData.addBoundarySegmentIdxToLastEdge(boundarySegementIndices[edgeIdx].second);
  };
  const auto& idSet = gv.grid().globalIdSet();
  auto yaspID       = idSet.id(element);

  auto addNewBoundarySegementIdx = [&](const CurveIndex& index, double t1, double t2) {
    if (initFlag) {
      elementTrimData.addBoundarySegmentIdxToLastEdge(numBoundarySegments_);
      boundarySegmentsArchive_[yaspID].push_back(std::make_tuple(index, t1, t2, numBoundarySegments_));
      ++numBoundarySegments_;
    } else {
      auto fatherId    = idSet.id(Util::coarsestFather(element));
      auto archiveData = boundarySegmentsArchive_.at(fatherId);
      auto it          = std::ranges::find_if(archiveData, [&](const auto& item) {
        const auto& [i_, t1_, t2_, b_] = item;
        // TODO generalize for more than one segment with the same curve Idx
        return i_.loop == index.loop and i_.curve == index.curve;
      });
      assert(it != archiveData.end());
      elementTrimData.addBoundarySegmentIdxToLastEdge(std::get<3>(*it));
    }
  };

  using namespace Util;

  // State
  std::vector<FieldVector<ScalarType, dim>> foundVertices;
  CurveIndex currentCurveIdx = {std::numeric_limits<size_t>::infinity(), std::numeric_limits<size_t>::infinity(),
                                std::numeric_limits<size_t>::infinity()};
  double currentT            = std::numeric_limits<double>::infinity();

  for (const auto i : std::views::iota(0u, result.vertices_.size())) {
    auto vV1 = result.vertices_[i];
    auto vV2 = result.vertices_[nextEntity(i)];

    auto pt1 = vV1.pt;
    auto pt2 = vV2.pt;

    // First case edge is completly untrimmed
    if (vV1.isHost() and vV2.isHost()) {
      if (not isConsecutive(vV1.hostId(), vV2.hostId())) {
        // if indices are not consecutive, we add a New New Segment
        if (foundVertices.empty()) {
          currentCurveIdx = getTrimmingCurveIdxFromHost(vV1);
          auto [currentT, curvePoint] =
              Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), vV1.hostId(), pt1, corners);
          assertPoint(pt1, curvePoint);

          foundVertices.push_back(curvePoint);
        }
        auto [tParam, curvePoint] =
            Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), vV2.hostId(), pt2, corners);
        assertPoint(pt2, curvePoint);

        auto elementTrimmingCurve =
            Util::createTrimmingCurveSlice(patchTrimData.getCurve(currentCurveIdx), currentT, tParam);
        elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
        addNewBoundarySegementIdx(currentCurveIdx, currentT, tParam);

        currentT = std::numeric_limits<double>::infinity();

      } else {
        elementTrimData.addEdge(Util::giveEdgeIdx(vV1.hostId(), vV2.hostId()));
        checkAndAddBoundarySegmentIdx(Util::giveEdgeIdx(vV1.hostId(), vV2.hostId()));
        foundVertices.push_back({pt2.x, pt2.y});
      }
    }

    // Second case edge begins on a hostVertex and ends on a newVertex
    else if (vV1.isHost() and vV2.isNew()) {
      if (vV1.hostId() != vV2.edgeId()) {
        if (foundVertices.empty()) {
          currentCurveIdx = getTrimmingCurveIdxFromHost(vV1);
          auto [currentT, curvePoint] =
              Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), vV1.hostId(), pt1, corners);
          assertPoint(pt1, curvePoint);

          foundVertices.push_back(curvePoint);
        }
        auto [tParam, curvePoint] =
            Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), vV2.edgeId(), pt2, corners);
        assertPoint(pt2, curvePoint);

        auto elementTrimmingCurve =
            Util::createTrimmingCurveSlice(patchTrimData.getCurve(currentCurveIdx), currentT, tParam);
        elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
        addNewBoundarySegementIdx(currentCurveIdx, currentT, tParam);

        foundVertices.push_back(curvePoint);

        currentT = std::numeric_limits<double>::infinity();
      } else {
        currentCurveIdx = getTrimmingCurveIdx(vV2);
        auto [tParam, curvePoint] =
            Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), vV2.edgeId(), pt2, corners);
        assertPoint(pt2, curvePoint);

        FieldVector<ScalarType, dim> p =
            foundVertices.empty() ? FieldVector<ScalarType, dim>{pt1.x, pt1.y} : foundVertices.back();
        auto trimmedEdge = Util::createHostGeometry<TrimmingCurve>(p, curvePoint);
        elementTrimData.addEdgeHostNew(vV2.edgeId(), trimmedEdge, curvePoint);
        checkAndAddBoundarySegmentIdx(vV2.edgeId());
        foundVertices.push_back(curvePoint);

        currentT = tParam;
      }
    }
    // Third case newVertex - newVertex
    else if (vV1.isNew() && vV2.isNew()) {
      // If there is no Vertex yet found add search for vV1 and add it to the foundVertices, also set currentCurveIdx
      if (foundVertices.empty()) {
        currentCurveIdx = getTrimmingCurveIdx(vV1);
        auto [currentT, curvePoint] =
            Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), vV1.edgeId(), pt1, corners);
        assertPoint(pt1, curvePoint);

        foundVertices.push_back(curvePoint);
      }
      if (getTrimmingCurveIdx(vV2).curve != currentCurveIdx.curve and
          getTrimmingCurveIdx(vV2).loop != currentCurveIdx.loop)
        throwGridError();

      auto [tParam, curvePoint] =
          Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), vV2.edgeId(), pt2, corners);
      assertPoint(pt2, curvePoint);

      // If currentT > tParam it can have 2 reasons. 1) the trim id only on one side of the element, and it has to
      // connect back or 2) if the trimming curve consist only of one curve where it has the same front and back
      // Controlpoint, so sometimes the wrong tParam (front or back) gets determined
      if (currentT > tParam) {
        bool success = false;
        if (currentCurveIdx.loop > 0 && patchTrimData.loops()[currentCurveIdx.loop].size() == 1) {
          auto curve = patchTrimData.getCurve(currentCurveIdx);
          if (curve.isConnectedAtBoundary(0)) {
            auto elementTrimmingCurve = Util::createTrimmingCurveSlice(curve, currentT, curve.domain()[0].back());
            elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
            addNewBoundarySegementIdx(currentCurveIdx, currentT, curve.domain()[0].back());

            currentT = std::numeric_limits<double>::infinity();
            success  = true;
          }
        } else if (vV1.edgeId() == vV2.edgeId()) {
          auto trimmedEdge = Util::createHostGeometry<TrimmingCurve>(foundVertices.back(), curvePoint);
          elementTrimData.addEdgeNewNewOnHost(vV2.edgeId(), trimmedEdge, curvePoint);
          checkAndAddBoundarySegmentIdx(vV2.edgeId());

          currentT = tParam;
          success  = true;
        }
        if (not success)
          throwGridError();
      } else {
        auto elementTrimmingCurve =
            Util::createTrimmingCurveSlice(patchTrimData.getCurve(currentCurveIdx), currentT, tParam);
        elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
        addNewBoundarySegementIdx(currentCurveIdx, currentT, tParam);

        currentT = std::numeric_limits<double>::infinity();
      }

      foundVertices.push_back(curvePoint);
    }

    // Fourth case edge begins on a newVertex and ends in a HostVertex
    else if (vV1.isNew() and vV2.isHost()) {
      if (not isConsecutive(vV1.edgeId(), vV2.hostId())) {
        auto [tParam, curvePoint] =
            Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), vV2.hostId(), pt2, corners);
        assertPoint(pt2, curvePoint);

        auto elementTrimmingCurve =
            Util::createTrimmingCurveSlice(patchTrimData.getCurve(currentCurveIdx), currentT, tParam);
        elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
        addNewBoundarySegementIdx(currentCurveIdx, currentT, tParam);

      } else {
        auto v2 = Dune::FieldVector<ScalarType, dim>{pt2.x, pt2.y};

        FieldVector<ScalarType, dim> p;
        if (foundVertices.empty())
          std::tie(std::ignore, p) =
              Util::callFindIntersection(patchTrimData.getCurve(getTrimmingCurveIdx(vV1)), vV1.edgeId(), pt1, corners);
        else
          p = foundVertices.back();

        auto trimmedEdge = Util::createHostGeometry<TrimmingCurve>(p, v2);
        elementTrimData.addEdgeNewHost(vV1.edgeId(), trimmedEdge, vV2.hostId());
        checkAndAddBoundarySegmentIdx(vV1.edgeId());

        foundVertices.push_back(v2);
      }
    }
    // Additional cases to cover inside vertices
    else if (vV1.isNew() and vV2.isInside()) {
      if (foundVertices.empty()) {
        currentCurveIdx = getTrimmingCurveIdx(vV1);
        FieldVector<ScalarType, dim> curvePoint;
        std::tie(currentT, curvePoint) =
            Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), vV1.edgeId(), pt1, corners);
        assertPoint(pt1, curvePoint);

        foundVertices.push_back(curvePoint);
      }

      const auto& curve         = patchTrimData.loops()[vV2.loop()].curves()[vV2.formerCurve()];
      double tParam             = curve.domain().front().back();
      auto elementTrimmingCurve = Util::createTrimmingCurveSlice(curve, currentT, tParam);
      auto curvePoint           = curve.corner(1);
      elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
      addNewBoundarySegementIdx(CurveIndex{.loop = (int)vV2.loop(), .curve = (int)vV2.formerCurve()}, currentT, tParam);
      foundVertices.push_back(curvePoint);
    } else if (vV1.isInside() and vV2.isNew()) {
      const auto& curve         = patchTrimData.loops()[vV1.loop()].curves()[vV1.subsequentCurve()];
      currentT                  = curve.domain().front().front();
      auto [tParam, curvePoint] = Util::callFindIntersection(curve, vV2.edgeId(), pt2, corners);
      assertPoint(pt2, curvePoint);

      auto elementTrimmingCurve = Util::createTrimmingCurveSlice(curve, currentT, tParam);
      elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
      addNewBoundarySegementIdx(CurveIndex{.loop = (int)vV1.loop(), .curve = (int)vV1.subsequentCurve()}, currentT,
                                tParam);

      foundVertices.push_back(curvePoint);
    } else {
      throwGridError();
    }
  }
  elementTrimData.finalize();
  return std::move(elementTrimData);
}
} // namespace Dune::IGA::DefaultTrim

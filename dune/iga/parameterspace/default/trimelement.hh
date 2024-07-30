// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <clipper2/clipper.h>

#include <dune/iga/geometrykernel/findintersection.hh>
#include <dune/iga/geometrykernel/slicecurve.hh>
#include <dune/iga/parameterspace/default/elementtrimdata.hh>
#include <dune/iga/parameterspace/default/parameterspace.hh>
#include <dune/iga/parameterspace/default/trimmingutils/clipelementrectangle.hh>
#include <dune/iga/parameterspace/default/trimmingutils/trimutils.hh>

namespace Dune::IGA::DefaultParameterSpace {

// TODO move to Util namespace and rename
template <typename ParameterSpace>
class LoopHandler
{
  using CurveIndex      = Impl::CurveLoopIndexEncoder::IndexResult;
  using PatchTrimData   = typename ParameterSpace::PatchTrimData;
  using ElementTrimData = typename ParameterSpace::ElementTrimData;
  using TrimmingCurve   = typename ParameterSpace::TrimmingCurve;

  static constexpr int dim             = ParameterSpace::mydimension;
  static constexpr int numberOfCorners = 4;

  using ScalarType                    = typename ParameterSpace::ctype;
  using Vertex                        = Util::ClippingResult::Vertex;
  using BoundarySegmentIndexContainer = std::array<std::pair<bool, size_t>, numberOfCorners>;
  using Point                         = Clipper2Lib::PointD;

public:
  template <typename GV>
  LoopHandler(const PatchTrimData& data, ElementTrimData& trimData, bool init, ParameterSpace* paraSpace, GV&& gridView,
              const auto& ele)
      : patchTrimData_(data),
        elementTrimData_(trimData),
        initFlag_(init),
        parameterSpace_(paraSpace),
        idSet(gridView.grid().globalIdSet()),
        element(ele) {
    auto geo = element.geometry();

    for (const auto i : Dune::range(numberOfCorners))
      corners_[i] = geo.corner(Util::vertexIndexMapping[i]);

    setUpBoundarySegmentIndices(std::forward<GV>(gridView));
  }

  void processVertexPair(const Vertex& vV1, const Vertex& vV2) {
    if (vV1.isHost() && vV2.isHost()) {
      processHostToHost(vV1, vV2);
    } else if (vV1.isHost() && vV2.isNew()) {
      processHostToNew(vV1, vV2);
    } else if (vV1.isNew() && vV2.isNew()) {
      processNewToNew(vV1, vV2);
    } else if (vV1.isNew() && vV2.isHost()) {
      processNewToHost(vV1, vV2);
    } else {
      if (not processAdditionalCases(vV1, vV2))
        throwGridError();
    }
  }

private:
  // see file trimminutils/trimelement.inl
  void processHostToHost(const Vertex& vV1, const Vertex& vV2);
  void processHostToNew(const Vertex& vV1, const Vertex& vV2);
  void processNewToNew(const Vertex& vV1, const Vertex& vV2);
  void processNewToHost(const Vertex& vV1, const Vertex& vV2);
  bool processAdditionalCases(const Vertex& vV1, const Vertex& vV2);
  auto findIntersection(const Vertex& vertex, const Point& pt);
  void processNonConsecutiveHosts(const Vertex& vV1, const Vertex& vV2, const Point& pt1, const Point& pt2);
  void processConsecutiveHosts(const Vertex& vV1, const Vertex& vV2, const Point& pt2);

  void throwGridError() const {
    DUNE_THROW(Dune::GridError, "TrimElement wasn't successful. Could not connect element trim curves");
  }
  void checkAndAddBoundarySegmentIdx(int edgeIdx) {
    if (boundarySegementIndices_[edgeIdx].first)
      elementTrimData_.addBoundarySegmentIdxToLastEdge(boundarySegementIndices_[edgeIdx].second);
  }
  auto findSegmentWithSameCurveIdx(const auto& archiveData, const CurveIndex& index) {
    return std::ranges::find_if(archiveData, [&](const auto& item) {
      const auto& [i_, t1_, t2_, b_] = item;
      return i_.loop == index.loop and i_.curve == index.curve;
    });
  }

  void addNewBoundarySegementIdx(const CurveIndex& index, double t1, double t2) {
    if (initFlag_) {
      elementTrimData_.addBoundarySegmentIdxToLastEdge(parameterSpace_->numBoundarySegments_);
      parameterSpace_->boundarySegmentsArchive_[yaspID].push_back(
          std::make_tuple(index, t1, t2, parameterSpace_->numBoundarySegments_));
      ++parameterSpace_->numBoundarySegments_;
    } else {
      auto fatherId    = idSet.id(Util::coarsestFather(element));
      auto archiveData = parameterSpace_->boundarySegmentsArchive_.at(fatherId);
      auto it          = findSegmentWithSameCurveIdx(archiveData, index);
      assert(it != archiveData.end());
      elementTrimData_.addBoundarySegmentIdxToLastEdge(std::get<3>(*it));
    }
  };

  template <typename GV>
  void setUpBoundarySegmentIndices(GV&& gridView) {
    const auto& idSet = gridView.grid().globalIdSet();
    yaspID            = idSet.id(element);
    for (const auto& intersection : Dune::intersections(gridView, element)) {
      auto idx = Transformations::mapToParameterSpace(1, intersection.indexInInside());
      if (not intersection.boundary())
        boundarySegementIndices_[idx] = std::make_pair(false, std::numeric_limits<size_t>::max());
      else
        boundarySegementIndices_[idx] = std::make_pair(true, intersection.boundarySegmentIndex());
    }
  }

  CurveIndex getTrimmingCurveIdx(const Vertex& vV) {
    return patchTrimData_.getIndices(vV.zValue());
  }
  std::vector<CurveIndex> getTrimmingCurveIdxFromHost(const Vertex& vV) {
    std::vector<CurveIndex> indices{};
    std::ranges::transform(vV.additionalZValues(), std::back_inserter(indices),
                           [&](auto&& zVal) { return patchTrimData_.getIndices(zVal); });
    return indices;
  }

  void assertPoint(const Clipper2Lib::PointD& ip, const FieldVector<ScalarType, dim>& curvePt, double prec = 0.01) {
    if (!Util::approxSamePoint(ip, curvePt, prec)) {
      std::cout << "Found ClipperPoint " << ip << " does not coincide with determined point on Curve " << curvePt
                << std::endl;
      throwGridError();
    }
  }

  // Member variables

  // Input Variables
  const PatchTrimData& patchTrimData_;
  ElementTrimData& elementTrimData_;
  bool initFlag_;
  std::array<FieldVector<double, dim>, numberOfCorners> corners_;

  // State Variables
  std::vector<FieldVector<ScalarType, dim>> foundVertices_;
  CurveIndex currentCurveIdx_ = {std::numeric_limits<size_t>::infinity(), std::numeric_limits<size_t>::infinity(),
                                 std::numeric_limits<size_t>::infinity()};
  double currentT_            = std::numeric_limits<double>::infinity();

  // Helpers
  BoundarySegmentIndexContainer boundarySegementIndices_;

  // Pointer to parameterSpace Grid (non-owning)
  ParameterSpace* parameterSpace_;

  // Refs to identities
  const typename ParameterSpace::UntrimmedParameterSpaceGrid::GlobalIdSet& idSet;
  typename ParameterSpace::UntrimmedParameterSpaceGrid::GlobalIdSet::IdType yaspID;
  const typename ParameterSpace::template YASPEntity<0>& element;
};

template <int dim, int dimworld, typename ScalarType>
auto ParameterSpaceImpl<dim, dimworld, ScalarType>::trimElement(const YASPEntity<0>& element, const auto& gv,
                                                                const PatchTrimData& patchTrimData,
                                                                bool initFlag) -> ElementTrimData {
  auto [flag, result] = Util::clipElementRectangle(element, patchTrimData);

  ElementTrimData elementTrimData(flag, element);
  if (flag != ElementTrimFlag::trimmed)
    return elementTrimData;

  LoopHandler<ParameterSpaceImpl> handler(patchTrimData, elementTrimData, initFlag, this, gv, element);

  for (const auto i : std::views::iota(0u, result.vertices_.size())) {
    auto vV1 = result.vertices_[i];
    auto vV2 = result.vertices_[(i + 1) % result.vertices_.size()];
    handler.processVertexPair(vV1, vV2);
  }

  elementTrimData.finalize();
  return elementTrimData;
}
} // namespace Dune::IGA::DefaultParameterSpace
#include "trimmingutils/trimelement.inl"

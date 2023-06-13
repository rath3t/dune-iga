// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/vtk/datacollectors/unstructureddatacollector.hh>
#include <dune/vtk/function.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/utility/lagrangepoints.hh>

namespace Dune::Vtk {
  /// Implementation of \ref Discontinuous DataCollector for Iga cells with trimming information
  template <class GridView>
  requires(GridView::dimension == 2) class DiscontinuousIgaDataCollector
      : public UnstructuredDataCollectorInterface<GridView, DiscontinuousIgaDataCollector<GridView>, Partitions::All> {
    using Self  = DiscontinuousIgaDataCollector;
    using Super = UnstructuredDataCollectorInterface<GridView, Self, Partitions::All>;

   public:
    using Super::dim;
    using Super::partition;

   public:
    explicit DiscontinuousIgaDataCollector(GridView const& gridView, int subSample = 0) : Super(gridView) {
      // FIXME this destroys the original trimmed element representation, if someone's wants to continue calculating
      // with this grid?
      for (auto element : elements(gridView_)) {
        element.impl().trimmedElementRepresentation()->refineAndConstructGrid(subSample);
      }
    }

    /// Construct the point sets
    void updateImpl() {
      pointSets_.clear();

      auto const& indexSet        = gridView_.indexSet();
      std::uint64_t vertexCounter = 0;
      numCells_                   = 0;
      numPoints_                  = 0;
      for (auto element : elements(gridView_)) {
        const size_t elementId = indexSet.index(element);

        auto elementRepr       = element.impl().trimmedElementRepresentation();
        auto subGridView       = elementRepr->gridView();
        auto verticesInSubGrid = subGridView.size(dim);

        numCells_ += subGridView.size(0);
        const auto& subGridIndexSet = subGridView.indexSet();

        for (auto gt : subGridIndexSet.types(0))
          pointSets_.try_emplace(gt, 1);

        for (auto& [type, pointSet] : pointSets_)
          if (pointSet.size() == 0) pointSet.build(type);

        for (auto subGridVertex : vertices(subGridView)) {
          std::size_t idx = subGridIndexSet.index(subGridVertex);
          vertexIndex_.emplace(std::array<std::size_t, 2>({elementId, idx}), vertexCounter++);
        }
        numPoints_ += verticesInSubGrid;
      }
    }

    /// Return number of Lagrange nodes
    [[nodiscard]] std::uint64_t numPointsImpl() const { return numPoints_; }

    /// Return a vector of point coordinates.
    /**
     * The vector of point coordinates is composed of vertex coordinates of the untrimmed elements and
     * the vertices of the triangulated trimmed elements
     **/
    template <class T>
    [[nodiscard]] std::vector<T> pointsImpl() const {
      std::vector<T> data(this->numPoints() * 3);
      auto const& indexSet = gridView_.indexSet();
      for (auto element : elements(gridView_, partition)) {
        auto geometry          = element.geometry();
        const size_t elementId = indexSet.index(element);

        auto elementRepr = element.impl().trimmedElementRepresentation();
        auto subGridView = elementRepr->gridView();

        const auto& trimmedIndexSet = subGridView.indexSet();
        for (auto subGridVertex : vertices(subGridView)) {
          auto vecInLocal = subGridVertex.geometry().center();
          auto v          = geometry.global(vecInLocal);

          std::int64_t idx = 3 * vertexIndex_.at({elementId, trimmedIndexSet.index(subGridVertex)});

          for (std::size_t j = 0; j < v.size(); ++j)
            data[idx + j] = T(v[j]);
          for (std::size_t j = v.size(); j < 3u; ++j)
            data[idx + j] = T(0);
        }
      }

      return data;
    }

    /// Return number of grid cells
    [[nodiscard]] std::uint64_t numCellsImpl() const { return numCells_; }

    /// \brief Return cell types, offsets, and connectivity. \see Cells
    /**
     * The cell connectivity is composed of cell vertices
     **/
    [[nodiscard]] Cells cellsImpl() const {
      Cells cells;
      cells.connectivity.reserve(this->numPoints());
      cells.offsets.reserve(this->numCells());
      cells.types.reserve(this->numCells());

      auto const& indexSet = gridView_.indexSet();

      std::int64_t old_o = 0;
      for (auto const& ele : elements(gridView_, partition)) {
        const std::size_t elementId = indexSet.index(ele);

        auto elementRepr = ele.impl().trimmedElementRepresentation();
        auto subGridView = elementRepr->gridView();

        const auto& subGridIndexSet = subGridView.indexSet();
        for (auto subGridElement : elements(subGridView)) {
          Vtk::CellType cellType(subGridElement.type(), Vtk::CellType::LAGRANGE);

          auto const& pointSet = pointSets_.at(subGridElement.type());

          for (std::size_t i = 0; i < pointSet.size(); ++i) {
            auto const& p        = pointSet[i];
            auto const& localKey = p.localKey();
            std::int64_t idx
                = vertexIndex_.at({elementId, subGridIndexSet.subIndex(subGridElement, localKey.subEntity(), dim)});

            cells.connectivity.push_back(idx);
          }
          cells.types.push_back(cellType.type());
          cells.offsets.push_back(old_o += pointSet.size());
        }
      }
      return cells;
    }

    /// Evaluate the `fct` at element vertices and edge centers in the same order as the point coords.
    template <class T, class GlobalFunction>
    [[nodiscard]] std::vector<T> pointDataImpl(GlobalFunction const& fct) const {
      int nComps = fct.numComponents();
      std::vector<T> data(this->numPoints() * nComps);

      auto localFct        = localFunction(fct);
      auto const& indexSet = gridView_.indexSet();
      for (auto element : elements(gridView_, partition)) {
        localFct.bind(element);
        auto geometry          = element.geometry();
        const size_t elementId = indexSet.index(element);

        auto elementRepr = element.impl().trimmedElementRepresentation();
        auto subGridView = elementRepr->gridView();

        const auto& subGridIndexSet = subGridView.indexSet();
        for (auto subGridVertex : vertices(subGridView)) {
          auto vecInLocal  = subGridVertex.geometry().center();
          std::int64_t idx = nComps * vertexIndex_.at({elementId, subGridIndexSet.index(subGridVertex)});

          for (std::size_t comp = 0; comp < nComps; ++comp)
            data[idx + comp] = T(localFct.evaluate(comp, vecInLocal));
        }
      }

      return data;
    }
    // Evaluate `fct` in center of cell.
    template <class T, class VtkFunction>
    [[nodiscard]] std::vector<T> cellDataImpl(VtkFunction const& fct) const {
      int nComps = fct.numComponents();
      std::vector<T> data;
      data.reserve(this->numCells_ * nComps);

      auto localFct        = localFunction(fct);
      auto const& indexSet = gridView_.indexSet();

      for (auto element : elements(gridView_, partition)) {
        localFct.bind(element);
        auto geometry          = element.geometry();
        const size_t elementId = indexSet.index(element);

        auto elementRepr            = element.impl().trimmedElementRepresentation();
        auto subGridView            = elementRepr->gridView();
        const auto& trimmedIndexSet = subGridView.indexSet();

        for (auto subGridElement : elements(subGridView)) {
          auto vecInLocal = subGridElement.geometry().center();

          for (std::size_t comp = 0; comp < nComps; ++comp)
            data.push_back(localFct.evaluate(comp, vecInLocal));
        }
      }
      return data;
    }

   protected:
    using Super::gridView_;

    std::uint64_t numPoints_ = 0;
    std::uint64_t numCells_  = 0;

    using PointSet = LagrangePointSet<typename GridView::ctype, GridView::dimension>;
    std::map<GeometryType, PointSet> pointSets_;
    std::vector<std::int64_t> indexMap_;
    std::map<std::array<std::size_t, 2>, std::int64_t> vertexIndex_;
  };

}  // namespace Dune::Vtk

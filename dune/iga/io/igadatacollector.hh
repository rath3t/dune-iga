// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "igarefinedgeometries.hh"

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/vtk/datacollectors/unstructureddatacollector.hh>
#include <dune/vtk/function.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/utility/lagrangepoints.hh>

namespace Dune::Vtk {
  /// Implementation of \ref Discontinuous DataCollector for Iga cells with or without trimming information
  template <class GridView>
  requires(GridView::dimension == 2) class DiscontinuousIgaDataCollector
      : public UnstructuredDataCollectorInterface<GridView, DiscontinuousIgaDataCollector<GridView>, Partitions::All> {
    using Self  = DiscontinuousIgaDataCollector;
    using Super = UnstructuredDataCollectorInterface<GridView, Self, Partitions::All>;

   public:
    using Super::dim;
    using Super::partition;

   public:
    DiscontinuousIgaDataCollector(GridView const& gridView, int subSampleFull, int subSampleTrimmed)
        : Super(gridView), geometries_(gridView, subSampleFull, subSampleTrimmed) {}
    // Does not subsample
    explicit DiscontinuousIgaDataCollector(GridView const& gridView) : DiscontinuousIgaDataCollector(gridView, 0, 0){};

    // Sub-samples trimmed elements by creating a new grid
    DiscontinuousIgaDataCollector(GridView const& gridView, int subsample)
        : DiscontinuousIgaDataCollector(gridView, subsample, subsample){};

    /// Construct the point sets
    void updateImpl() {
      pointSets_.clear();

      auto const& indexSet        = this->gridView().indexSet();
      std::uint64_t vertexCounter = 0;
      numCells_                   = 0;
      numPoints_                  = 0;
      for (auto element : elements(this->gridView())) {
        const size_t elementId = indexSet.index(element);

        auto verticesInSubGrid = geometries_.nVertices(elementId);

        numCells_ += geometries_.nElements(elementId);

        pointSets_.try_emplace(geometries_.geometryType(elementId), 1);

        for (auto& [type, pointSet] : pointSets_)
          if (pointSet.size() == 0) pointSet.build(type);

        for (std::size_t vIdx = 0; auto& subGridVertex : geometries_.getVertices(elementId))
          vertexIndex_.emplace(std::array<std::size_t, 2>({elementId, vIdx++}), vertexCounter++);

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
      auto const& indexSet = this->gridView().indexSet();
      for (auto element : elements(this->gridView(), partition)) {
        auto geometry          = element.geometry();
        const size_t elementId = indexSet.index(element);

        for (std::uint64_t vIdx = 0; auto& subGridVertex : geometries_.getVertices(elementId)) {
          auto v = geometry.global(subGridVertex);

          std::int64_t idx = 3 * vertexIndex_.at({elementId, vIdx++});

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

      auto const& indexSet = this->gridView().indexSet();

      for (std::int64_t old_o = 0; auto const& ele : elements(this->gridView(), partition)) {
        const std::size_t elementId = indexSet.index(ele);

        for (std::size_t eIdx = 0; auto& subGridElement : geometries_.getElements(elementId)) {
          Vtk::CellType cellType(subGridElement.type(), Vtk::CellType::LAGRANGE);

          auto const& pointSet = pointSets_.at(subGridElement.type());

          for (std::size_t i = 0; i < pointSet.size(); ++i) {
            auto const& p        = pointSet[i];
            auto const& localKey = p.localKey();
            std::int64_t idx
                = vertexIndex_.at({elementId, geometries_.vertexSubIndex(elementId, eIdx, localKey.subEntity())});

            cells.connectivity.push_back(idx);
          }
          cells.types.push_back(cellType.type());
          cells.offsets.push_back(old_o += pointSet.size());
          ++eIdx;
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
      auto const& indexSet = this->gridView().indexSet();
      for (auto element : elements(this->gridView(), partition)) {
        localFct.bind(element);
        auto geometry          = element.geometry();
        const size_t elementId = indexSet.index(element);

        for (std::uint64_t vIdx = 0; auto& subGridVertex : geometries_.getVertices(elementId)) {
          std::int64_t idx = nComps * vertexIndex_.at({elementId, vIdx++});

          for (std::size_t comp = 0; comp < nComps; ++comp)
            data[idx + comp] = T(localFct.evaluate(comp, subGridVertex));
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
      auto const& indexSet = this->gridView().indexSet();

      for (auto element : elements(this->gridView(), partition)) {
        localFct.bind(element);
        auto geometry          = element.geometry();
        const size_t elementId = indexSet.index(element);

        for (auto& subGridElement : geometries_.getElements(elementId)) {
          auto vecInLocal = subGridElement.center();

          for (std::size_t comp = 0; comp < nComps; ++comp)
            data.push_back(localFct.evaluate(comp, vecInLocal));
        }
      }
      return data;
    }

   protected:
#if DUNE_VERSION_LTE(DUNE_VTK, 2, 10)
    using Super::gridView_;
    const auto& gridView() const { return gridView_; }
#endif

    std::uint64_t numPoints_ = 0;
    std::uint64_t numCells_  = 0;

    using PointSet = LagrangePointSet<typename GridView::ctype, GridView::dimension>;
    std::map<GeometryType, PointSet> pointSets_;
    std::vector<std::int64_t> indexMap_;
    std::map<std::array<std::size_t, 2>, std::int64_t> vertexIndex_;
    IGA::IGARefinedGeometries<GridView> geometries_;
  };

}  // namespace Dune::Vtk

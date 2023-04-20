#pragma once

#include <cassert>
#include <map>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/vtk/datacollectors/unstructureddatacollector.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/utility/lagrangepoints.hh>

namespace Dune {
  namespace Vtk {
    /// Implementation of \ref Discontinuous DataCollector for Iga cells with trimming information
    template <class GridView>
    class DiscontinuousIgaDataCollector
        : public UnstructuredDataCollectorInterface<GridView, DiscontinuousIgaDataCollector<GridView>,
                                                    Partitions::All> {
      using Self  = DiscontinuousIgaDataCollector;
      using Super = UnstructuredDataCollectorInterface<GridView, Self, Partitions::All>;

     public:
      using Super::dim;
      using Super::partition;

     public:
      DiscontinuousIgaDataCollector(GridView const& gridView, int subSample = 1)
          : Super(gridView), subSample_(subSample) {}

      /// Construct the point sets
      void updateImpl() {
        pointSets_.clear();

        auto const& indexSet        = gridView_.indexSet();
        std::uint64_t vertexCounter = 0;
        numCells_                   = 0;
        numPoints_                  = 0;
        for (auto element : elements(gridView_)) {
          const size_t elementId = indexSet.index(element);

          if (element.impl().isTrimmed()) {
            assert(element.impl()
                       .trimmedElementRepresentation()
                       .has_value());  // at this point there should be a trimming information
            auto elementRepr               = element.impl().trimmedElementRepresentation().value();
            auto trimmedGridView           = elementRepr->gridView();
            auto verticesInTrimmedGridView = trimmedGridView.size(dim);
            numCells_ += trimmedGridView.size(0);
            const auto& trimmedIndexSet = trimmedGridView.indexSet();
            for (auto gt : trimmedIndexSet.types(0))
              pointSets_.try_emplace(gt, 1);

            for (auto& [type, pointSet] : pointSets_)
              if (pointSet.size() == 0) pointSet.build(type);

            for (auto triangulationVertices : vertices(trimmedGridView)) {
              std::size_t idx = trimmedIndexSet.index(triangulationVertices);
              vertexIndex_.emplace(std::array<std::size_t, 2>({elementId, idx}), vertexCounter++);
            }
            numPoints_ += verticesInTrimmedGridView;
          } else {
            assert(!element.impl()
                        .trimmedElementRepresentation()
                        .has_value());  // at this point there should be a trimming information
            auto verticesInTrimmedGridView = 4;
            ++numCells_;
            pointSets_.try_emplace(Dune::GeometryTypes::quadrilateral, 1);

            for (auto& [type, pointSet] : pointSets_)
              if (pointSet.size() == 0) pointSet.build(type);

            auto const& pointSet = pointSets_.at(Dune::GeometryTypes::quadrilateral);

            for (std::size_t i = 0; i < pointSet.size(); ++i) {
              auto const& p        = pointSet[i];
              auto const& localKey = p.localKey();
              vertexIndex_.emplace(std::array<std::size_t, 2>({elementId, localKey.subEntity()}), vertexCounter++);
            }

            numPoints_ += verticesInTrimmedGridView;
          }
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
          if (element.impl().isTrimmed()) {
            assert(element.impl()
                       .trimmedElementRepresentation()
                       .has_value());  // at this point there should be a trimming information
            auto elementRepr               = element.impl().trimmedElementRepresentation().value();
            auto trimmedGridView           = elementRepr->gridView();
            auto verticesInTrimmedGridView = trimmedGridView.size(dim);
            const auto& trimmedIndexSet    = trimmedGridView.indexSet();
            for (auto triangulationVertex : vertices(trimmedGridView)) {
              auto trimmedGeometry = triangulationVertex.geometry().center();

              std::int64_t idx = 3 * vertexIndex_.at({elementId, trimmedIndexSet.index(triangulationVertex)});
              auto vecInLocal  = geometry.impl().spanToLocal(trimmedGeometry);
              auto v           = geometry.global(vecInLocal);
              for (std::size_t j = 0; j < v.size(); ++j)
                data[idx + j] = T(v[j]);
              for (std::size_t j = v.size(); j < 3u; ++j)
                data[idx + j] = T(0);
            }
          } else {
            assert(!element.impl()
                        .trimmedElementRepresentation()
                        .has_value());  // at this point there should be a trimming information

            auto const& pointSet = pointSets_.at(Dune::GeometryTypes::quadrilateral);

            for (std::size_t i = 0; i < pointSet.size(); ++i) {
              auto const& p        = pointSet[i];
              auto const& localKey = p.localKey();
              std::int64_t idx     = 3 * vertexIndex_.at({elementId, localKey.subEntity()});
              auto v               = geometry.corner(localKey.subEntity());
              for (std::size_t j = 0; j < v.size(); ++j)
                data[idx + j] = T(v[j]);
              for (std::size_t j = v.size(); j < 3u; ++j)
                data[idx + j] = T(0);
            }
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

          if (ele.impl().isTrimmed()) {
            assert(ele.impl()
                       .trimmedElementRepresentation()
                       .has_value());  // at this point there should be a trimming information
            auto elementRepr               = ele.impl().trimmedElementRepresentation().value();
            auto trimmedGridView           = elementRepr->gridView();
            auto verticesInTrimmedGridView = trimmedGridView.size(dim);
            const auto& trimmedIndexSet    = trimmedGridView.indexSet();
            for (auto triangulationElement : elements(trimmedGridView)) {
              Vtk::CellType cellType(triangulationElement.type(), Vtk::CellType::LAGRANGE);

              auto const& pointSet = pointSets_.at(triangulationElement.type());

              for (std::size_t i = 0; i < pointSet.size(); ++i) {
                auto const& p        = pointSet[i];
                auto const& localKey = p.localKey();
                std::int64_t idx     = vertexIndex_.at(
                    {elementId, trimmedIndexSet.subIndex(triangulationElement, localKey.subEntity(), dim)});

                cells.connectivity.push_back(idx);
              }
              cells.types.push_back(cellType.type());
              cells.offsets.push_back(old_o += pointSet.size());
            }
          } else {
            Vtk::CellType cellType(Dune::GeometryTypes::quadrilateral, Vtk::CellType::LAGRANGE);

            auto const& pointSet = pointSets_.at(ele.type());

            for (std::size_t i = 0; i < pointSet.size(); ++i) {
              auto const& p        = pointSet[i];
              auto const& localKey = p.localKey();
              std::int64_t idx     = vertexIndex_.at({elementId, localKey.subEntity()});

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
          if (element.impl().isTrimmed()) {
            assert(element.impl()
                       .trimmedElementRepresentation()
                       .has_value());  // at this point there should be a trimming information
            auto elementRepr               = element.impl().trimmedElementRepresentation().value();
            auto trimmedGridView           = elementRepr->gridView();
            auto verticesInTrimmedGridView = trimmedGridView.size(dim);
            const auto& trimmedIndexSet    = trimmedGridView.indexSet();
            for (auto triangulationVertex : vertices(trimmedGridView)) {
              auto trimmedGeometry = triangulationVertex.geometry().center();

              std::int64_t idx = nComps * vertexIndex_.at({elementId, trimmedIndexSet.index(triangulationVertex)});
              auto vecInLocal  = geometry.impl().spanToLocal(trimmedGeometry); //transform vertex position from iga span to [0,1]^d
              for (std::size_t comp = 0; comp < nComps; ++comp)
                data[idx + comp] = T(localFct.evaluate(comp, vecInLocal));
            }
          } else {
            assert(!element.impl()
                        .trimmedElementRepresentation()
                        .has_value());  // at this point there should be no trimming information
            auto const& pointSet = pointSets_.at(Dune::GeometryTypes::quadrilateral);

            for (std::size_t i = 0; i < pointSet.size(); ++i) {
              auto const& p        = pointSet[i];
              auto const& localKey = p.localKey();
              std::int64_t idx     = nComps * vertexIndex_.at({elementId, localKey.subEntity()});
              for (std::size_t comp = 0; comp < nComps; ++comp)
                data[idx + comp] = T(localFct.evaluate(comp, p.point()));
            }
          }
        }

        return data;
      }

     protected:
      using Super::gridView_;

      unsigned int subSample_;
      std::uint64_t numPoints_ = 0;
      std::uint64_t numCells_  = 0;

      using PointSet = LagrangePointSet<typename GridView::ctype, GridView::dimension>;
      std::map<GeometryType, PointSet> pointSets_;
      std::vector<std::int64_t> indexMap_;
      std::map<std::array<std::size_t, 2>, std::int64_t> vertexIndex_;
    };

  }  // end namespace Vtk
}  // end namespace Dune

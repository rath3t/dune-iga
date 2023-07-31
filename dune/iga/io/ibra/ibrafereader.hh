//
// Created by Henri on 28.07.2023.
//

#pragma once
#include "ibrageometry.hh"
#include "ibrareader.hh"

#include <unordered_set>

#include "dune/iga/io/igarefinedgeometries.hh"
#include "dune/iga/trim/nurbstrimboundary.hh"
#include <dune/curvedgrid/curvedgrid.hh>
#include <dune/curvedgrid/gridfunctions/analyticgridfunction.hh>
#include <dune/vtk/datacollectors/lagrangedatacollector.hh>

namespace Dune::IGA {

  // At the moment only gridDim und worldDim == 2 or 3 supported
  template <typename Grid>
  class IbraFEReader {
   public:
    static const int gridDim  = Grid::dimension;
    static const int worldDim = Grid::dimensionworld;
    using ctype               = Grid::ctype;
    using ControlPoint        = Dune::IGA::NURBSPatchData<gridDim, worldDim>::ControlPointType;
    using PatchData           = Dune::IGA::NURBSPatchData<gridDim, worldDim>;
    using PatchGeometry       = Dune::IGA::NURBSPatchGeometry<gridDim, worldDim>;
    using IGAGrid             = Dune::IGA::NURBSGrid<gridDim, worldDim, ctype>;

    using ControlPointNetType = Dune::IGA::MultiDimensionNet<gridDim, ControlPoint>;

    struct Transformer {
      PatchGeometry* patchGeometry;

      template <typename VecType>
      auto operator()(const VecType& cp) const -> VecType {
        return patchGeometry->global(cp);
      }
    };

    struct Comparator {
      bool operator()(const Dune::FieldVector<double, 2>& p1, const Dune::FieldVector<double, 2>& p2) const {
        if (p1[0] < p2[0]) return true;
        if (p1[0] > p2[0]) return false;
        return p1[1] < p2[1];
      };
    };

    class IGAProjection {
      using Domain   = FieldVector<ctype, worldDim>;
      using Jacobian = FieldMatrix<ctype, worldDim, worldDim>;

      PatchGeometry patchGeometry_;

     public:
      IGAProjection(PatchGeometry& patchGeometry) : patchGeometry_(patchGeometry) {}

      /// \brief project the coordinate from parameter space to
      Domain operator()(const Domain& x) const { return patchGeometry_.global(x); }

      /// \brief derivative of the projection
      friend auto derivative(const IGAProjection& igaProjection) {
        return [patchGeometry = igaProjection.patchGeometry_](const Domain& x) {
          return patchGeometry.jacobianTransposed(x);
        };
      }
    };

    static auto read(const std::string& fileName, int refine, bool trim = true) {
      std::ifstream ibraInputFile;
      ibraInputFile.open(fileName);

      auto brep = getBrep<decltype(ibraInputFile), 2>(ibraInputFile);

      Ibra::Surface<worldDim> _surface                           = brep.surfaces[0];
      const std::array<std::vector<double>, gridDim> knotSpans   = _surface.compileKnotVectors();
      const std::vector<std::vector<ControlPoint>> controlPoints = _surface.transformControlPoints();
      std::array<int, gridDim> dimsize                           = _surface.n_controlPoints;

      auto controlNet = ControlPointNetType(dimsize, controlPoints);
      PatchData patchData{knotSpans, controlNet, _surface.degree};
      PatchGeometry patchGeometry{std::make_shared<PatchData>(patchData)};
      Transformer transformer{&patchGeometry};

      // Create IGA Grid
      auto trimData = constructGlobalBoundaries(brep);
      auto igaGrid  = trim ? std::make_shared<IGAGrid>(patchData, trimData) : std::make_shared<IGAGrid>(patchData);
      igaGrid->globalRefine(refine);
      auto igaGridView = igaGrid->leafGridView();

      // Gather vertices & elements
      std::set<Dune::FieldVector<ctype, gridDim>, Comparator> vertices;

      const auto& indexSet = igaGridView.indexSet();
      auto geometries      = IGARefinedGeometries(igaGridView, 0, 0);

      for (auto& element : elements(igaGridView)) {
        const size_t elementId = indexSet.index(element);
        auto geometry          = element.geometry();

        for (auto& subGridVertex : geometries.getVertices(elementId))
          vertices.insert(geometry.impl().localToSpan(subGridVertex));
      }

      auto gridFactory = Dune::GridFactory<Grid>();
      // Add vertices
      for (auto& vertex : vertices)
        gridFactory.insertVertex(vertex);

      // Reconstruct elements with global Vertex Indices

      for (auto& element : elements(igaGridView)) {
        auto geometry          = element.geometry();
        const size_t elementId = indexSet.index(element);
        auto gt                = geometries.geometryType(elementId);
        unsigned int nSubI     = gt == GeometryTypes::simplex(2) ? 3 : 4;

        auto& eleVertices = geometries.getVertices(elementId);

        for (auto subEleIdx : std::views::iota(0u, geometries.nElements(elementId))) {
          std::vector<unsigned int> elementVertices;

          for (auto subEntityIndex : std::views::iota(0u, nSubI)) {
            auto localVertexIdx                       = geometries.vertexSubIndex(elementId, subEleIdx, subEntityIndex);
            Dune::FieldVector<double, gridDim> vertex = geometry.impl().localToSpan(eleVertices[localVertexIdx]);

            // Find Idx
            auto it = vertices.find(vertex);
            assert(it != vertices.end());

            elementVertices.push_back(std::distance(vertices.begin(), it));
          }
          gridFactory.insertElement(gt, elementVertices);
        }
      }

      auto gridFunction = analyticGridFunction<Grid>(IGAProjection(patchGeometry));

      // Wrap the reference curvedGrid to build a curved curvedGrid
      auto refGrid = gridFactory.createGrid();

      // pass in refGrid by reference, curvedGrid is not owner of refGrid
      CurvedGrid curvedGrid{*refGrid, gridFunction};

      // The moment the refGrid goes out of scope, the curvedGrid is invalidated
      return std::make_pair<decltype(refGrid), decltype(curvedGrid)>(std::move(refGrid), std::move(curvedGrid));
    }
  };
}  // namespace Dune::IGA

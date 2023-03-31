//
// Created by Henri on 08.03.2023.
//

#ifndef IKARUS_RECONSTRUCTEDGRIDHANDLER_H
#define IKARUS_RECONSTRUCTEDGRIDHANDLER_H

#include <mapbox/earcut.hpp>

#include <dune/grid/uggrid.hh>
#include <dune/iga/nurbspatchgeometry.h>

// Add support for Dune::Fieldvector in Earcut
namespace mapbox::util {

  template <typename T>
  struct nth<0, Dune::FieldVector<T, 2>> {
    inline static auto get(const Dune::FieldVector<T, 2>& t) { return t[0]; };
  };

  template <typename T>
  struct nth<1, Dune::FieldVector<T, 2>> {
    inline static auto get(const Dune::FieldVector<T, 2>& t) { return t[1]; };
  };
}  // namespace mapbox::util

namespace Dune::IGA {

  template <int worldDim>
  struct GridBoundarySegment : Dune::BoundarySegment<2, worldDim, double> {
    // Types

    // Constructor
    explicit GridBoundarySegment(Boundary& _boundary, auto _transferToGlobal)
        : boundary(_boundary), transferToGlobal{_transferToGlobal} {}

    Dune::FieldVector<double, worldDim> operator()(const Dune::FieldVector<double, 1>& local) const override {
      // u muss auf die Domäne gemappt werden, local[0] ist zwischen 0 und 1
      double u_local = local[0];
      double u       = Utilities::map<double>(u_local, 0.0, 1.0, boundary.domain[0], boundary.domain[1]);
      auto res       = boundary.nurbsGeometry(u);

      // Transfer to global Surface Coordinates
      return transferToGlobal(res);
    }

    Boundary boundary;
    NURBSPatchGeometry<2, worldDim> transferToGlobal;
  };

  template <int worldDim>
  class TrimmedElementRepresentation {
   public:
    using Grid     = Dune::UGGrid<worldDim>;
    using GridView = Grid::LeafGridView;

    // Abkürzungen
    using BoundaryVector = std::vector<Boundary>;
    using LocalPoint     = Dune::FieldVector<double, 2>;
    using GlobalPoint    = Dune::FieldVector<double, worldDim>;

    using LocalPointVector  = std::vector<LocalPoint>;
    using GlobalPointVector = std::vector<GlobalPoint>;

    // Grid
    std::unique_ptr<Grid> grid;

   private:
    BoundaryVector boundaries;
    NURBSPatchGeometry<2, worldDim> patchGeometry;
    std::array<int, 2> gridDegree;

   public:
    // Construct with igaGrid and boundaries
    TrimmedElementRepresentation(BoundaryVector& _boundaries, auto _patchGeometry, std::array<int, 2> _gridDegree)
        : boundaries(_boundaries), patchGeometry(_patchGeometry), gridDegree(_gridDegree) {
      reconstructTrimmedElement();
    }

   private:
    void reconstructTrimmedElement() {
      // Getting global Parameters
      auto parameters = Utilities::getParameters();

      if (parameters.preSample > 0) splitBoundaries();

      Dune::GridFactory<Grid> gridFactory;
      auto [indices, vertices] = createMesh();

      for (auto& vertex : vertices)
        gridFactory.insertVertex(vertex);

      int n_ele = (int)indices.size() / 3;
      for (int i = 0; i < n_ele; ++i) {
        gridFactory.insertElement(Dune::GeometryTypes::triangle, {indices[3 * i], indices[3 * i + 1], indices[3 * i + 2]});
      }

      for (auto& boundary : boundaries) {
        auto idx = getControlPointIndices(vertices, boundary);
        gridFactory.insertBoundarySegment({idx[0], idx[1]},
                                          std::make_shared<GridBoundarySegment<worldDim>>(boundary, patchGeometry));
      }

      auto boundaryToRefineMap = determineBoundariesToRefine();

      // Create Grid
      grid = gridFactory.createGrid();

      GridView gridView = grid->leafGridView();

      if (parameters.preGlobalRefine > 0) grid->globalRefine(parameters.preGlobalRefine);

      for (int i = 0; i < parameters.edgeRefinements; ++i) {
        for (const auto& ele : elements(grid->leafGridView())) {
          bool mark = false;
          if (ele.hasBoundaryIntersections())
            for (auto& intersection : intersections(gridView, ele))
              if (intersection.boundary())
                if (boundaryToRefineMap[intersection.boundarySegmentIndex()]) mark = true;
          if (mark) grid->mark(1, ele);
        }

        grid->preAdapt();
        grid->adapt();
        grid->postAdapt();
      }

      GridView gridViewRefined = grid->leafGridView();

      double areaRefinedGrid
          = std::accumulate(elements(gridViewRefined).begin(), elements(gridViewRefined).end(), 0.0,
                            [](double rhs, const auto& element) { return rhs + element.geometry().volume(); });

      std::cout << "Reconstructed Grid with " << gridViewRefined.size(0)
                << " elements. Area of elements: " << areaRefinedGrid
                << ". Approx area of trimmed element: " << calculateArea() << std::endl;
    }


    std::pair<std::vector<unsigned int>, std::vector<GlobalPoint>> createMesh() {
      // Construct mesh with Earcut
      // C.f. https://github.com/mapbox/earcut.hpp

      // Create array of points
      GlobalPointVector vertices;
      for (auto& boundary : boundaries)
        vertices.push_back(patchGeometry({boundary.endPoints.front()}));

      std::vector<std::vector<LocalPoint>> polygonInput;
      polygonInput.push_back(vertices);

      std::vector<unsigned int> indices = mapbox::earcut<unsigned int>(polygonInput);

      return {indices, vertices};
    }

    double calculateArea(unsigned int div = 200) {
      Clipper2Lib::PathD polygon;
      for (auto& boundary : boundaries) {
        auto path = boundary.path(div, false);
        for (Clipper2Lib::PointD point : path) {
          auto globalPoint = patchGeometry({point.x, point.y});
          polygon.emplace_back(globalPoint[0], globalPoint[1]);
        }
      }
      return std::fabs(Clipper2Lib::Area(polygon));
    }

    bool approxSamePoint(const LocalPoint& a, const LocalPoint& b) { return Dune::FloatCmp::eq(a, b, 1e-8); };

    void splitBoundaries() {
      auto parameters = Utilities::getParameters();

      std::map<int, bool> refineMap = determineBoundariesToRefine();

      using Domain            = std::array<double, 2>;
      using DomainInformation = std::pair<Domain, int>;

      // Get all domains
      std::vector<DomainInformation> domains;
      for (int i = 0; i < boundaries.size(); ++i) {
        domains.emplace_back(boundaries[i].domain, i);
      }

      std::vector<DomainInformation> tempDomains;

      for (int i = 0; i < parameters.preSample; ++i) {
        tempDomains.clear();
        tempDomains.insert(tempDomains.end(), domains.begin(), domains.end());
        domains.clear();
        for (auto& domain : tempDomains) {
          bool refine = !parameters.preSampleOnlyCurvedEdges;
          if (parameters.preSampleOnlyCurvedEdges && refineMap[domain.second]) refine = true;

          if (refine) {
            std::array<Domain, 2> newDomains = Utilities::splitDomains(domain.first);
            domains.emplace_back(newDomains[0], domain.second);
            domains.emplace_back(newDomains[1], domain.second);
          } else
            domains.push_back(domain);
        }
      }
      tempDomains.clear();

      // Create Splittet trimData
      std::vector<Boundary> newBoundaries;
      for (const auto& domain_info : domains) {
        const int i = domain_info.second;
        newBoundaries.emplace_back(boundaries[i].nurbsGeometry, domain_info.first);
      }
      boundaries = newBoundaries;
    }

    std::array<unsigned int, 2> getControlPointIndices(GlobalPointVector& vertices, Boundary& boundary) {
      auto it1 = std::find_if(vertices.begin(), vertices.end(),
                              [&](const GlobalPoint& v) { return approxSamePoint(v, patchGeometry(boundary.endPoints.front())); });
      if (it1 == vertices.end()) throw std::runtime_error("Reconstruction of Grid failed: Vertex 1 not found");

      auto it2 = std::find_if(vertices.begin(), vertices.end(),
                              [&](const LocalPoint& v) { return approxSamePoint(v, patchGeometry(boundary.endPoints.back())); });
      if (it2 == vertices.end()) throw std::runtime_error("Reconstruction of Grid failed: Vertex 2 not found");

      const unsigned int idx1 = std::distance(vertices.begin(), it1);
      const unsigned int idx2 = std::distance(vertices.begin(), it2);

      assert(idx1 != idx2);
      return {idx1, idx2};
    }

    std::map<int, bool> determineBoundariesToRefine() {
      std::map<int, bool> boundaryToRefineMap;

      for (int i = 0; const auto& boundary : boundaries) {
        // Determine whether to refine this boundary based on its degree and orientation
        bool shouldRefine = false;
        if (boundary.degree() > 1) {
          shouldRefine = true;
        } else {
          Boundary::EdgeOrientation orientation = boundary.getEdgeOrientation();
          int degree = (orientation == Boundary::EdgeOrientation::u) ? gridDegree[0] : gridDegree[1];
          if (degree > 1 || orientation == Boundary::EdgeOrientation::Unknown) {
            shouldRefine = true;
          }
          if (gridDegree[0] == 1 && gridDegree[1] == 1)
            shouldRefine = false;
        }
        boundaryToRefineMap[i] = shouldRefine;
        ++i;
      }
      return boundaryToRefineMap;
    }
  };
}  // namespace Dune::IGA

#endif  // IKARUS_RECONSTRUCTEDGRIDHANDLER_H

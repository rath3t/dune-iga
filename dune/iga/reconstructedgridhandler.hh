//
// Created by Henri on 08.03.2023.
//

#ifndef IKARUS_RECONSTRUCTEDGRIDHANDLER_H
#define IKARUS_RECONSTRUCTEDGRIDHANDLER_H

#include <mapbox/earcut.hpp>
#include <dune/iga/nurbspatchgeometry.h>

#include <dune/grid/uggrid.hh>

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
  class ReconstructedGridHandler {
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
    ReconstructedGridHandler(BoundaryVector& _boundaries, auto _patchGeometry, std::array<int, 2> _gridDegree)
        : boundaries(_boundaries), patchGeometry(_patchGeometry), gridDegree(_gridDegree) {
      reconstructTrimmedElement();
    }

   private:
    void reconstructTrimmedElement() {
      // Getting global Parameters
      auto parameters = Utilities::getParameters();

      if (parameters.preSample > 0) splitBoundaries();

      LocalPointVector vertices;
      for (auto& boundary : boundaries) {
        vertices.push_back(boundary.controlPoints[0]);
      }

      GlobalPointVector globalVertices;
      for (auto& point : vertices) {
        globalVertices.push_back(patchGeometry({point[0], point[1]}));
      }

      // Grid Factory
      Dune::GridFactory<Grid> gridFactory;

      for (auto& vertex : globalVertices)
        gridFactory.insertVertex(vertex);

      // Make Grid Template
      auto meshTemplate = createOriginalMeshTemplate<unsigned int>();

      int n_ele = (int)meshTemplate.size() / 3;

      for (int i = 0; i < n_ele; ++i) {
        gridFactory.insertElement(Dune::GeometryTypes::triangle,
                                  {meshTemplate[3 * i], meshTemplate[3 * i + 1], meshTemplate[3 * i + 2]});
      }

      for (auto& boundary : boundaries) {
        auto idx = getControlPointIndices(vertices, boundary);
        gridFactory.insertBoundarySegment({idx[0], idx[1]}, std::make_shared<GridBoundarySegment<worldDim>>(boundary, patchGeometry));
      }

      // Setze Idx boundarySegmentIdx welche refined werden müssen
      auto boundaryToRefineMap = determineBoundariesToRefine();

      // Calculate the area of the trimmed element
      // auto area = calculateArea<double>();

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

      std::cout << "This gridview contains: ";
      std::cout << gridViewRefined.size(0) << " elements" << std::endl;

      double areaRefinedGrid = 0.0;
      for (auto& element : elements(gridViewRefined))
        areaRefinedGrid += element.geometry().volume();
      std::cout << "The area of all elements in the refined grid is: " << areaRefinedGrid << std::endl;
    }

    template <typename N>
    std::vector<N> createOriginalMeshTemplate() {
      // Construct mesh with Earcut
      // C.f. https://github.com/mapbox/earcut.hpp

      // Create array of points
      LocalPointVector polygon;
      for (auto& boundary : boundaries)
        polygon.push_back(boundary.controlPoints[0]);

      std::vector<std::vector<LocalPoint>> polygonInput;
      polygonInput.push_back(polygon);

      std::vector<N> indices = mapbox::earcut<N>(polygonInput);

      return indices;
    }

    template <typename T>
    T calculateArea(unsigned int div = 200) {
      using namespace Clipper2Lib;

      Path<T> polygon;
      for (auto& boundary : boundaries) {
        auto path = boundary.path<T>(div);
        polygon.insert(polygon.end(), path.begin(), path.end());
      }

      return Area(polygon);
    }

    bool approximatelySamePoint(const LocalPoint& a, const LocalPoint& b) {
      const double tolerance = 16.0 * std::numeric_limits<double>::epsilon();

      return std::abs(a[0] - b[0]) < tolerance && std::abs(a[1] - b[1]) < tolerance;
    };

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

      // Create Splittet boundaries
      std::vector<Boundary> newBoundaries;
      for (const auto& domain_info : domains) {
        const int i = domain_info.second;
        newBoundaries.emplace_back(boundaries[i].nurbsGeometry, domain_info.first);
      }
      boundaries = newBoundaries;
    }

    std::array<unsigned int, 2> getControlPointIndices(LocalPointVector& vertices, Boundary& boundary) {
      auto it1 = std::find_if(vertices.begin(), vertices.end(), [&](const LocalPoint& v) {
        return approximatelySamePoint(v, boundary.controlPoints.front());
      });
      if (it1 == vertices.end()) {
        std::cerr << "Vertex 1 not found" << std::endl;
        return {0, 0};
      }
      const unsigned int idx1 = std::distance(vertices.begin(), it1);

      auto it2 = std::find_if(vertices.begin(), vertices.end(), [&](const LocalPoint& v) {
        return approximatelySamePoint(v, boundary.controlPoints.back());
      });
      if (it2 == vertices.end()) {
        std::cerr << "Vertex 2 not found" << std::endl;
        return {0, 0};
      }
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
          Boundary::EdgeOrientation orientation = boundary.getOrientation();
          int degree = (orientation == Boundary::EdgeOrientation::u) ? gridDegree[0] : gridDegree[1];
          if (degree > 1 || orientation == Boundary::EdgeOrientation::Unknown) {
            shouldRefine = true;
          }
        }
        boundaryToRefineMap[i] = shouldRefine;
        ++i;
      }
      return boundaryToRefineMap;
    }
  };
}  // namespace Dune::IGA

#endif  // IKARUS_RECONSTRUCTEDGRIDHANDLER_H

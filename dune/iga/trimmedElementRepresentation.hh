//
// Created by Henri on 08.03.2023.
//

#pragma once

#include <bits/ranges_algo.h>
#include <bits/ranges_algobase.h>
#include <clipper2/clipper.core.h>
#include <cstddef>
#include <iterator>
#include <mapbox/earcut.hpp>
#include <vector>

#include <dune/alugrid/grid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/iga/nurbstrimboundary.hh>
#include <dune/iga/nurbstrimfunctionality.hh>

// Add support for Dune::FieldVector in Earcut
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

  template <int dim>
  struct GridBoundarySegment : Dune::BoundarySegment<2, dim, double> {
    explicit GridBoundarySegment(Boundary& _boundary) : boundary(_boundary) {}

    Dune::FieldVector<double, dim> operator()(const Dune::FieldVector<double, 1>& local) const override {
      // u has to be mapped on the domain of 0 to 1
      // Interestingly enough ALUGrid asks for values less than 0, e.g. -1e-16, and greater than 1 ??
      double u_local = local[0];
      if (!Dune::FloatCmp::ge(u_local, 0.0)) u_local = 0.0;
      if (!Dune::FloatCmp::le(u_local, 1.0)) u_local = 1.0;

      double u = Utilities::mapToRange(u_local, 0.0, 1.0, boundary.domain[0], boundary.domain[1]);
      return boundary.nurbsGeometry(u);
    }

    Boundary boundary;
  };

  /** \brief representation of the trimmed element in the parameter space */
  template <int dim, typename Grid = Dune::UGGrid<dim>>
  class TrimmedElementRepresentation {
   public:
    using GridView = Grid::LeafGridView;
    using Point    = Dune::FieldVector<double, dim>;

    // Grid
    std::unique_ptr<Grid> grid{};

   private:
    std::vector<Boundary> boundaries;
    std::vector<Boundary> originalBoundaries;
    bool trimmed;
    bool verbose = false;

   public:
    // Construct with boundaries
    explicit TrimmedElementRepresentation(std::vector<Boundary>& _boundaries) : boundaries(_boundaries), trimmed(true) {
      std::ranges::copy(boundaries, std::back_inserter(originalBoundaries));
      reconstructTrimmedElement();
    }
    /// brief: Constructs an untrimmed elementRepresentation, expects indices in Dune ordering
    explicit TrimmedElementRepresentation(std::vector<FieldVector<double, 2>>& corners) : trimmed(false) {
        constructFromCorners(corners);
    }

    // Accessors
    GridView gridView() const { return grid->leafGridView(); }
    [[nodiscard]] bool isTrimmed() const {
      return trimmed;
    }
    bool alreadyRefined = false;
    GridView refinedGridView(unsigned int refinement = 0) {
      if (refinement > 0 and not alreadyRefined) {
        grid->globalRefine(refinement);
        alreadyRefined = true;
      }
      return grid->leafGridView();
    }

   private:
    void reconstructTrimmedElement() {
      // Getting global Parameters
      auto parameters = Utilities::getParameters();

      if (parameters.preSample > 0) splitBoundaries();

      Dune::GridFactory<Grid> gridFactory;
      auto [indices, vertices] = triangulate();

      for (auto& vertex : vertices)
        gridFactory.insertVertex(vertex);

      // The element indices are stored as a flat vector, 3 indices always constitute 1 triangle
      for (auto it = indices.begin(); it < indices.end(); it += 3) {
        gridFactory.insertElement(Dune::GeometryTypes::triangle, std::vector<unsigned int>(it, std::next(it, 3)));

        #if 0
        // Check for counter-clockwise orientation
        auto eleIdx = std::vector<unsigned int>(it, std::next(it, 3));
        Clipper2Lib::PathD path;
        std::ranges::transform(eleIdx, std::back_inserter(path), [&](unsigned int idx) -> Clipper2Lib::PointD{
          return {vertices[idx][0], vertices[idx][1]};
        });
        assert(Clipper2Lib::IsPositive(path));
        #endif
      }
      for (auto& boundary : boundaries)
        gridFactory.insertBoundarySegment(getControlPointIndices(vertices, boundary),
                                          std::make_shared<GridBoundarySegment<dim>>(boundary));

      // Create Grid
      grid = gridFactory.createGrid();
      if (parameters.preGlobalRefine > 0) grid->globalRefine(parameters.preGlobalRefine);

      refineGridOnEdges(parameters.edgeRefinements);

      GridView gridViewRefined = grid->leafGridView();

      if (verbose)
        std::cout << "Reconstructed Grid with " << gridViewRefined.size(0)
                << " elements. Area of elements: " << calculateArea()
                << ". Approx area of trimmed element: " << calculateTargetArea() << std::endl;
    }

    void refineGridOnEdges(int refCount) {
      GridView gridView       = grid->leafGridView();
      auto boundariesToRefine = determineBoundariesToRefine();

      for (int i = 0; i < refCount; ++i) {
        for (const auto& ele : elements(grid->leafGridView())) {
          bool mark = false;
          if (ele.hasBoundaryIntersections())
            for (auto& intersection : intersections(gridView, ele))
              if (intersection.boundary())
                if (boundariesToRefine[intersection.boundarySegmentIndex()]) mark = true;
          if (mark) grid->mark(1, ele);
        }

        grid->preAdapt();
        grid->adapt();
        grid->postAdapt();
      }
    }

    std::pair<std::vector<unsigned int>, std::vector<Point>> triangulate() {
      // Construct mesh with Earcut
      // C.f. https://github.com/mapbox/earcut.hpp

      // Create array of points
      std::vector<Point> vertices;
      for (auto& boundary : boundaries)
        vertices.push_back({boundary.endPoints.front()});

      std::vector<std::vector<Point>> polygonInput;
      polygonInput.push_back(vertices);

      std::vector<unsigned int> indices = mapbox::earcut<unsigned int>(polygonInput);

      return {indices, vertices};
    }

    void constructFromCorners(std::vector<FieldVector<double, 2>>& corners) {
      assert(corners.size() == 4);
      assert(corners[2][0] < corners[3][0]);

      Dune::GridFactory<Grid> gridFactory;

      for (auto& vertex : corners)
        gridFactory.insertVertex(vertex);

      gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 1, 2, 3});
      grid = gridFactory.createGrid();
    }

   public:
    // Obtain boundaries and vertexes
    [[nodiscard]] std::vector<Boundary> getBoundaries() const {
      return originalBoundaries;
    }
    [[nodiscard]] std::size_t numBoundaries() const {
      return originalBoundaries.size();
    }
    [[nodiscard]] Boundary getBoundaryAt(std::size_t index) const {
      return originalBoundaries.at(index);
    }
    [[nodiscard]] std::vector<Point> getVertexes() const {
      std::vector<Point> vertexes;
      std::ranges::transform(originalBoundaries, std::back_inserter(vertexes), [](const auto& boundary){
        return boundary.endPoints.front();
      });
    }
    [[nodiscard]] std::size_t numVertexes() const {
      return originalBoundaries.size();
    }
    [[nodiscard]] Point getVertexAt(std::size_t index) const {
      return originalBoundaries.at(index).endPoints.front();
    }

    /// \brief Calculates the area from the actual trim paths (might be expensive)
    double calculateTargetArea(unsigned int div = 200) {
      Clipper2Lib::PathD polygon;
      polygon.reserve(div*boundaries.size());


      for (auto& boundary : boundaries) {



        auto path = boundary.path(div);
        for (Clipper2Lib::PointD point : path) {
          polygon.emplace_back(point);
        }
      }
      return std::fabs(Clipper2Lib::Area(polygon));
    }

    /// \brief Calculates the area from the simplex elements in the current grid
    double calculateArea() {
      GridView gv = grid->leafGridView();
      return std::accumulate(elements(gv).begin(), elements(gv).end(), 0.0,
                             [](double rhs, const auto& element) { return rhs + element.geometry().volume(); });
    }

    /// \brief Checks the grid for overlaps, returns true if any element overlaps another
    bool checkGridForOverlappingElements() {
      GridView gV = grid->leafGridView();

      std::vector<Clipper2Lib::PointD> centers;
      std::ranges::transform(elements(gV), std::back_inserter(centers), [](const auto& ele) -> Clipper2Lib::PointD {
        auto c = ele.geometry().center();
        return {c[0], c[1]};
      });

      auto getElementEdges = [](const auto& geo) -> Clipper2Lib::PathD {
        std::vector<Point> corners{geo.corner(0), geo.corner(1), geo.corner(2)};
        Clipper2Lib::PathD edges;
        for (auto& c : corners)
          edges.emplace_back(c[0], c[1]);
        return edges;
      };
      std::vector<Clipper2Lib::PathD> elementEdges;
      std::ranges::transform(elements(gV), std::back_inserter(elementEdges),
                             [&getElementEdges](const auto& ele) { return getElementEdges(ele.geometry()); });

      int centerCounter = -1;
      int edgeCounter   = -1;
      auto hasOverlap   = std::ranges::any_of(centers, [&](const Clipper2Lib::PointD& center) {
        ++centerCounter;
        edgeCounter = -1;
        return std::ranges::any_of(elementEdges, [&](const Clipper2Lib::PathD& edges) {
          ++edgeCounter;
          if (edgeCounter == centerCounter) return false;
          auto res = Clipper2Lib::PointInPolygon(center, edges);
          return (res == Clipper2Lib::PointInPolygonResult::IsInside || res == Clipper2Lib::PointInPolygonResult::IsOn);
        });
      });

      return hasOverlap;
    }

   private:
    bool approxSamePoint(const Point& a, const Point& b) { return Dune::FloatCmp::eq(a, b, 1e-8); };

    void splitBoundaries() {
      auto parameters = Utilities::getParameters();

      auto refineMap = determineBoundariesToRefine();

      using Domain            = std::array<double, 2>;
      using DomainInformation = std::pair<Domain, int>;

      // Get all domains
      std::vector<DomainInformation> domains;
      domains.reserve(boundaries.size());
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

      // Create split boundaries from DomainInformation
      std::vector<Boundary> newBoundaries;
      std::ranges::transform(domains, std::back_inserter(newBoundaries), [&](const auto& domain_info) -> Boundary{
        return {boundaries[domain_info.second].nurbsGeometry, domain_info.first};
      });
      boundaries = std::move(newBoundaries);
    }

    std::vector<unsigned int> getControlPointIndices(std::vector<Point>& vertices, Boundary& boundary) {
      auto it1 = std::ranges::find_if(vertices,
                              [&](const Point& v) { return approxSamePoint(v, boundary.endPoints.front()); });

      auto it2 = std::ranges::find_if(vertices,
                              [&](const Point& v) { return approxSamePoint(v, boundary.endPoints.back()); });

      assert(it1 != vertices.end());
      assert(it2 != vertices.end());

      const unsigned int idx1 = std::distance(vertices.begin(), it1);
      const unsigned int idx2 = std::distance(vertices.begin(), it2);

      assert(idx1 != idx2);
      return {idx1, idx2};
    }

    std::vector<bool> determineBoundariesToRefine() {
      std::vector<bool> result;

      for (const auto& boundary : boundaries) {
        if (boundary.degree() > 1)
          result.push_back(true);
        else
          result.push_back(false);
      }
      return result;
    }
  };
}  // namespace Dune::IGA

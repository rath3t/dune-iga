//
// Created by Henri on 14.03.2023.
//
#pragma once

#include "ibraGeometry.hh"
#include "nurbstrimutils.hh"

#include <clipper2/clipper.core.h>

#include <dune/iga/nurbspatchgeometry.h>

namespace Dune::IGA {

  class Boundary {
   public:
    static constexpr int worldDim = 2;
    static constexpr int dim      = 1;

    using Geometry            = NURBSPatchGeometry<dim, worldDim>;
    using ControlPoint        = Geometry::ControlPointType;
    using ControlPointNetType = Dune::IGA::MultiDimensionNet<dim, ControlPoint>;
    using PatchData           = Dune::IGA::NURBSPatchData<dim, worldDim>;

    using Point  = Dune::FieldVector<double, worldDim>;

    // Member variables
    Geometry nurbsGeometry;
    std::array<double, 2> domain{};
    std::array<Point, 2> endPoints;

    Boundary() = default;

    Boundary(Geometry& _nurbsGeometry, std::array<double, 2> _domain)
        : nurbsGeometry(_nurbsGeometry),
          domain(_domain),
          endPoints({nurbsGeometry(domain[0]), nurbsGeometry(domain[1])}) {
      assert(domain[1] > domain[0]);
    }

    explicit Boundary(Ibra::BrepTrim& _trim)
        : nurbsGeometry(geometryFromTrim(_trim.geometry)),
          domain(_trim.domain),
          endPoints({nurbsGeometry(domain[0]), nurbsGeometry(domain[1])}) {}

    explicit Boundary(const Point& a, const Point& b)
        : nurbsGeometry(lineGeometryFromPoints(a, b)), domain(nurbsGeometry.domain()[0]), endPoints({a, b}) {}

    // Helper classes for construction of nurbsGeometry
   private:
    static Geometry geometryFromTrim(Ibra::Curve2D& _curve) {
      const auto _cp = _curve.transformControlPoints()[0];
      std::array<int, dim> dimSize{static_cast<int>(_cp.size())};

      // Construct patch Data
      std::array<std::vector<double>, dim> knotSpans {_curve.compileKnotVectors()};
      auto controlNet {ControlPointNetType(dimSize, _cp)};
      std::array<int, 1> degree {_curve.degree};

      return Geometry(std::make_shared<PatchData>(knotSpans, controlNet, degree));
    }

    /** \brief creates a line geometry from a to b */
    static Geometry lineGeometryFromPoints(const Point& a, const Point b) {
      std::vector<Point> _controlPoints{a, b};
      std::array<ControlPoint, worldDim> _cp{ControlPoint{.p{_controlPoints.front()}, .w = 1},
                                             {.p{_controlPoints.back()}, .w = 1}};
      std::array<int, 1> dimSize{2};

      // Construct patch Data
      std::array<std::vector<double>, dim> knotSpans {std::vector<double>{0.0, 0.0, 1.0, 1.0}};
      auto controlNet {ControlPointNetType(dimSize, _cp)};
      std::array<int, 1> degree {1};

      return Geometry(std::make_shared<PatchData>(knotSpans, controlNet, degree));
    }

   public:
    [[nodiscard]] int degree() const { return nurbsGeometry.degree()[0]; };

    enum class EdgeOrientation { u, v, Unknown };
    [[nodiscard]] EdgeOrientation getEdgeOrientation(const double tolerance = 1e-6) const {
      if (Dune::FloatCmp::eq(endPoints[0][0], endPoints[1][0], tolerance))
        return EdgeOrientation::v;
      else if (Dune::FloatCmp::eq(endPoints[0][1], endPoints[1][1], tolerance))
        return EdgeOrientation::u;
      else
        return EdgeOrientation::Unknown;
    }

    [[nodiscard]] Clipper2Lib::Path<double> path(unsigned int samples = 200, bool getOnlyTwoPointsIfStraight = true) const {
      Clipper2Lib::Path<double> path;

      if (getOnlyTwoPointsIfStraight && degree() == 1 && endPoints.size() == 2) {
        path.emplace_back(endPoints.front()[0], endPoints.front()[1]);
        path.emplace_back(endPoints.back()[0], endPoints.back()[1]);
        return path;
      }

      auto linS = Utilities::linspace<double>(domain[0], domain[1], samples);
      for (auto u : linS) {
        auto vertex = nurbsGeometry(u);
        path.emplace_back(vertex[0], vertex[1]);
      }
      return path;
    }

  };
  struct BoundaryLoop {
    enum class Orientation {ClockWise, CounterClockWise};
    std::vector<Boundary> boundaries;
    Orientation orientation;
  };

  class TrimData {
   public:
    std::vector<BoundaryLoop> boundaryLoops;

    TrimData() = default;

    void addLoop(std::vector<Boundary> _boundaries) {
      auto orientation = determineOrientation(_boundaries);
      boundaryLoops.push_back({_boundaries, orientation});
    }

    [[nodiscard]] auto getOrientationForBoundary(const Boundary& _boundary) const ->BoundaryLoop::Orientation {
      for (const auto& loop : boundaryLoops) {
        // Find boundary in loop
        auto sameBoundary = [&_boundary](auto checkBoundary){
          return _boundary.endPoints == checkBoundary.endPoints;
        };
        if (std::ranges::find_if(loop.boundaries, sameBoundary) != loop.boundaries.end())
          return loop.orientation;
      }
    }

   private:
    // TODO Maybe use Clipper::isPositive
    static auto determineOrientation(std::vector<Boundary>& _boundaries) -> BoundaryLoop::Orientation {
      // extract some vertices to test
      std::vector<Dune::FieldVector<double, 2>> vertices;
      if (_boundaries.size() == 1) {
        auto boundary = _boundaries[0];
        auto linSpace = Utilities::linspace(boundary.domain, 4);
        linSpace.pop_back();
        for (auto u : linSpace)
          vertices.push_back(boundary.nurbsGeometry(u));
      } else {
        for (const auto& boundary : _boundaries)
          vertices.push_back(boundary.endPoints[0]);
      }

      // C.f. https://stackoverflow.com/a/1180256 and http://www.faqs.org/faqs/graphics/algorithms-faq/ Subject 2.07
      // Find lowest and rightmost vertex
      auto it = std::ranges::min_element(vertices, [](auto v1, auto v2) {
        if (Dune::FloatCmp::eq(v1[1], v2[1]))
          return v1[0] > v2[0];
        return v1[1] < v2[1];
      });

      auto idxA = std::ranges::distance(vertices.begin(), it);
      auto idxB = (idxA + 1) % vertices.size();
      auto idxC = (idxA - 1) % vertices.size();

      auto a = vertices[idxA];
      auto b = vertices[idxB];
      auto c = vertices[idxC];

      auto ab = b-a;
      auto ac = c-a;

      auto isCrossPositive = [](Dune::FieldVector<double, 2>& a, Dune::FieldVector<double, 2>& b) {
        double cross =  a[0] * b[1] - b[0] * a[1];
        if (cross > 0)
          return true;
        else return false;
      };
      if (isCrossPositive(ab, ac))
        return BoundaryLoop::Orientation::CounterClockWise;
      else
        return BoundaryLoop::Orientation::ClockWise;
    }
  };

}  // namespace Dune::IGA

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
          endPoints({nurbsGeometry(domain[0]), nurbsGeometry(domain[1])}) {}

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

    /** \brief creates a line IbraBase from a to b */
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
    [[nodiscard]] EdgeOrientation getOrientation(const double tolerance = 1e-6) const {
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
}  // namespace Dune::IGA

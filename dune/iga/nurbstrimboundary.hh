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
    static constexpr int gridDim = 1;

    using Geometry = NURBSPatchGeometry<gridDim, worldDim>;
    using ControlPoint = Geometry::ControlPointType;
    using ControlPointNetType = Dune::IGA::MultiDimensionNet<gridDim, ControlPoint>;
    using PatchData = Dune::IGA::NURBSPatchData<gridDim, worldDim>;

    using Point       = Dune::FieldVector<double, worldDim>;
    using PointVector = std::vector<Point>;
    using Domain      = std::array<double, 2>;


    // Member variables
    Geometry nurbsGeometry;

    Domain domain{};

    // Control Points are the start and endpoint (Rename???)
    PointVector controlPoints;

    Boundary() = default;
    Boundary(Geometry& _nurbsGeometry, Domain _domain, PointVector& _controlPoints)
        : nurbsGeometry(_nurbsGeometry), domain(_domain), controlPoints(_controlPoints) {}

    Boundary(Geometry& _nurbsGeometry, Domain _domain) : nurbsGeometry(_nurbsGeometry), domain(_domain) {
      controlPoints = {nurbsGeometry(domain[0]), nurbsGeometry(domain[1])};
    }

    explicit Boundary(Ibra::BrepTrim& _trim) {
      domain        = _trim.domain;
      nurbsGeometry = geometryFromTrim(_trim.geometry);
      controlPoints = {nurbsGeometry(domain[0]), nurbsGeometry(domain[1])};
    }

    Boundary(Ibra::BrepTrim& _trim, Domain _domain) : domain(_domain) {
      nurbsGeometry = geometryFromTrim(_trim.geometry);
      controlPoints = {nurbsGeometry(domain[0]), nurbsGeometry(domain[1])};
    }

    explicit Boundary(PointVector& _controlPoints) : controlPoints(_controlPoints) {
      assert(controlPoints.size() == 2);

      nurbsGeometry = geometryFromPoints(_controlPoints);
      domain        = nurbsGeometry.domain()[0];
    }
    // Helper classes for construction of nurbsGeometry
   private:
    static Geometry geometryFromTrim(Ibra::Curve& _curve) {

      const auto _cp = _curve.compileControlPoints();

      std::array<std::vector<double>, gridDim> knotSpans;
      knotSpans[0] = _curve.compileKnotVectors();

      // Dim Size is the amount of controlPoints
      std::array<int, gridDim> dimsize  {(int)_cp.size()};

      // Construct patch Data
      auto controlNet = ControlPointNetType(dimsize, _cp);
      PatchData _patchData;

      _patchData.knotSpans = knotSpans;
      _patchData.degree = {_curve.degree};
      _patchData.controlPoints = controlNet;

      return Geometry(std::make_shared<PatchData>(_patchData));
    }

    static Geometry geometryFromPoints(PointVector& _controlPoints) {
      std::array<ControlPoint, worldDim> _cp {ControlPoint{.p{_controlPoints.front()}, .w = 1},
                                         {.p{_controlPoints.back()}, .w = 1}};

      std::array<std::vector<double>, gridDim> knotSpans;
      knotSpans[0] = {0, 0, 1, 1};

      std::array<int, 1> dimsize = {2};

      // Construct patch Data
      auto controlNet = ControlPointNetType(dimsize, _cp);
      PatchData _patchData;

      _patchData.knotSpans = knotSpans;
      _patchData.degree = {1};
      _patchData.controlPoints = controlNet;

      return Geometry(std::make_shared<PatchData>(_patchData));
    }

   public:
    [[nodiscard]] int degree() const { return nurbsGeometry.degree()[0]; };

    enum EdgeOrientation { u, v, Unknown };
    [[nodiscard]] EdgeOrientation getOrientation() const {
      // const double tolerance = double(16) * std::numeric_limits<double>::epsilon();
      const double tolerance = 1e-6;
      if (std::fabs(controlPoints[0][0] - controlPoints[1][0]) < tolerance)
        return EdgeOrientation::v;
      else if (std::fabs(controlPoints[0][1] - controlPoints[1][1]) < tolerance)
        return EdgeOrientation::u;
      else
        return EdgeOrientation::Unknown;
    }

    template <typename T = double>
    Clipper2Lib::Path<T> path(unsigned int samples = 200, bool getOnlyTwoPointsIfStraight = true) {
      using namespace Clipper2Lib;

      Path<T> path;
      if (getOnlyTwoPointsIfStraight && degree() == 1 && controlPoints.size() == 2) {
        path.emplace_back(controlPoints.front()[0], controlPoints.front()[1]);
        path.emplace_back(controlPoints.back()[0], controlPoints.back()[1]);
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

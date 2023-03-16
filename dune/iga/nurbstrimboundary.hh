//
// Created by Henri on 14.03.2023.
//
#pragma once

#include "ibraGeometry.hh"
#include "nurbstrimutils.hh"

#include <clipper2/clipper.core.h>

namespace Dune::IGA {

  /** \brief this is a dummy of the geometry implementation for NURBS. Doenst depend on the NURBSGrid as template
   * parameter*/
  template <std::integral auto mydim, std::integral auto dimworld>
  class NURBSCurveGeometry {
   public:
    static constexpr std::integral auto mydimension = mydim;

    static constexpr std::integral auto coorddimension = dimworld;
    static constexpr std::integral auto griddim        = 1;

    using ctype               = typename DuneLinearAlgebraTraits<double>::value_type;
    using LinearAlgebraTraits = DuneLinearAlgebraTraits<double>;
    using LocalCoordinate     = typename LinearAlgebraTraits::template FixedVectorType<mydimension>;
    using GlobalCoordinate    = typename LinearAlgebraTraits::template FixedVectorType<coorddimension>;
    using JacobianTransposed  = typename LinearAlgebraTraits::template FixedMatrixType<mydimension, coorddimension>;
    using JacobianInverseTransposed =
        typename DuneLinearAlgebraTraits<double>::template FixedMatrixType<coorddimension, mydimension>;

    using ControlPointType = typename NURBSPatchData<griddim, dimworld, LinearAlgebraTraits>::ControlPointType;

   private:
    /* Helper class to compute a matrix pseudo inverse */
    using MatrixHelper = typename MultiLinearGeometryTraits<ctype>::MatrixHelper;

   public:
    /** \brief Constructor from NURBSPatchData and an iterator to a specific knot
     *
     *  \param Patchdata shared pointer to an object where the all the data of the NURBSPatch is stored
     *  \param fixedOrVaryingDirections indicates if the direction free or fixed. This means that the geometry does not
     * "run" in the fixed direction, e.g. an edge in the first direction is fixed in the second direction and a vertex
     * is fixed in all directions
     */
    NURBSCurveGeometry(std::shared_ptr<NURBSPatchData<griddim, dimworld, LinearAlgebraTraits>> patchData,
                       const std::array<Impl::FixedOrFree, griddim>& fixedOrVaryingDirections,
                       const std::array<int, griddim>& thisSpanIndices)
        : patchData_(patchData), fixedOrVaryingDirections_{fixedOrVaryingDirections} {
      for (int i = 0; i < griddim; ++i) {
        if (thisSpanIndices[i] + 1 < patchData_->knotSpans[i].size())
          scaling_[i] = patchData_->knotSpans[i][thisSpanIndices[i] + 1] - patchData_->knotSpans[i][thisSpanIndices[i]];
        offset_[i] = patchData_->knotSpans[i][thisSpanIndices[i]];
      }
      for (int i = 0; i < griddim; ++i)
        thisSpanIndices_[i] = (thisSpanIndices[i] == patchData->knotSpans[i].size() - 1)
                                  ? thisSpanIndices[i] - patchData->degree[i] - 1
                                  : thisSpanIndices[i];

      // If we are a vertex and on the rightmost end of the knotspan, we receive here the last index,
      //  For the proper construction of the nurbs and controlpoint net we need the indices end -degree -2
      //  To properly extract for the last span and not for the span after the last one
      if constexpr (mydim == 0)
        for (int i = 0; i < griddim; ++i)
          if (thisSpanIndices[i] == patchData->knotSpans[i].size() - 1)
            thisSpanIndices_[i] = patchData->knotSpans[i].size() - patchData->degree[i] - 2;

      nurbs_ = Dune::IGA::Nurbs<griddim, LinearAlgebraTraits>(*patchData, thisSpanIndices_);
      cpCoordinateNet_
          = netOfSpan(thisSpanIndices_, patchData_->degree, extractControlCoordinates(patchData_->controlPoints));
    }

    NURBSCurveGeometry() = default;

    [[nodiscard]] const auto& nurbs() const { return nurbs_; }

    [[nodiscard]] const auto& controlPoints() const { return cpCoordinateNet_; }

    std::shared_ptr<NURBSPatchData<griddim, dimworld, LinearAlgebraTraits>> patchData_;
    std::array<int, griddim> thisSpanIndices_;
    std::array<Impl::FixedOrFree, griddim> fixedOrVaryingDirections_{Impl::FixedOrFree::free};
    Dune::IGA::Nurbs<griddim, LinearAlgebraTraits> nurbs_;
    std::array<ctype, griddim> offset_;
    std::array<ctype, griddim> scaling_;
    MultiDimensionNet<griddim, typename ControlPointType::VectorType> cpCoordinateNet_;
  };

  class NurbsCurveHandler {
   public:
    static constexpr int worldDim = 2;
    static constexpr int gridDim  = 1;

    using PatchData    = NURBSPatchData<gridDim, worldDim>;
    using ControlPoint = PatchData::ControlPointType;

    using KnotSpanGeometry    = NURBSCurveGeometry<gridDim, worldDim>;
    using ControlPointNetType = Dune::IGA::MultiDimensionNet<gridDim, ControlPoint>;

    using JacobianTransposed = typename KnotSpanGeometry::JacobianTransposed;

    Dune::IGA::Nurbs<gridDim> nurbsCurve;
    std::shared_ptr<PatchData> patchData;

    NurbsCurveHandler() = default;

    explicit NurbsCurveHandler(Ibra::Curve2D& _curve) {
      const std::vector<ControlPoint> controlPoints = _curve.template compileControlPoints<ControlPoint>();

      std::array<std::vector<double>, gridDim> knotSpans;
      knotSpans[0] = _curve.compileKnotVectors();

      // Dim Size is the amount of controlPoints
      std::array<int, gridDim> dimsize = {(int)controlPoints.size()};

      // Construct patch Data
      auto controlNet = ControlPointNetType(dimsize, controlPoints);
      PatchData _patchData;

      _patchData.knotSpans     = knotSpans;
      _patchData.degree        = {_curve.degree};
      _patchData.controlPoints = controlNet;

      patchData  = std::make_shared<PatchData>(_patchData);
      nurbsCurve = Dune::IGA::Nurbs(_patchData);
    }

    /// Creates a linear curve with 2 ControlPoints
    explicit NurbsCurveHandler(std::array<ControlPoint, 2>& cp) {
      std::array<std::vector<double>, 1> knotSpans;
      knotSpans[0] = {0, 0, 1, 1};

      std::array<int, 1> dimsize = {2};

      // Construct patch Data
      auto controlNet = ControlPointNetType(dimsize, cp);
      PatchData _patchData;

      _patchData.knotSpans     = knotSpans;
      _patchData.degree        = {1};
      _patchData.controlPoints = controlNet;

      patchData  = std::make_shared<PatchData>(_patchData);
      nurbsCurve = Dune::IGA::Nurbs(_patchData);
    }

    [[nodiscard]] std::array<double, 2> domain() const {
      double min = std::min_element(patchData->knotSpans[0].begin(), patchData->knotSpans[0].end())[0];
      double max = std::max_element(patchData->knotSpans[0].begin(), patchData->knotSpans[0].end())[0];

      return {min, max};
    }

    [[nodiscard]] int degree() const { return patchData->degree[0]; }

    [[nodiscard]] double domainMidPoint() const {
      auto dom = domain();
      return (dom[0] + dom[1]) / 2;
    }

    [[nodiscard]] int n_controlPoints() const { return patchData->controlPoints.size()[0]; }

    Dune::FieldVector<double, 2> controlPointAt(int index) {
      ControlPoint cp = patchData->controlPoints.get({index});
      return cp.p;
    }

    [[nodiscard]] int findSpanIndex(double u) const {
      return static_cast<int>(Dune::IGA::findSpanUncorrected(patchData->degree[0], u, patchData->knotSpans[0]));
    }

    [[nodiscard]] auto geometry(int spanIndex) const {
      return KnotSpanGeometry(patchData, {Dune::IGA::Impl::FixedOrFree::free}, {spanIndex});
    }

    Dune::FieldVector<double, 2> operator()(double u) const {
      auto spanIndex = findSpanIndex(u);
      auto geo       = geometry(spanIndex);
      auto basis     = nurbsCurve.basisFunctionNet({u});

      return Dune::IGA::dot(basis, geo.cpCoordinateNet_);
    }

    [[nodiscard]] JacobianTransposed jacobianTransposed(double _local) const {
      std::array<double, 1> u = {_local};
      JacobianTransposed result;
      std::array<unsigned int, gridDim> subDirs{};

      for (int subI = 0, i = 0; i < gridDim; ++i) {
        // if (fixedOrVaryingDirections_[i] == Grid::FixedOrFree::fixed) continue;
        subDirs[subI++] = i;
      }
      auto spanIndex = findSpanIndex(_local);
      auto geo       = geometry(spanIndex);

      const auto basisFunctionDerivatives = nurbsCurve.basisFunctionDerivatives(u, 1);

      for (int dir = 0; dir < gridDim; ++dir) {
        std::array<unsigned int, gridDim> ithVec{};
        ithVec[subDirs[dir]] = 1;
        result[dir]          = dot(basisFunctionDerivatives.get(ithVec), geo.cpCoordinateNet_);
      }
      return result;
    }

    Dune::FieldVector<double, 2> operator()(double u, int deriv) const {
      assert(deriv > 0);
      assert(deriv < 2);

      auto jacobian = jacobianTransposed(u);

      Dune::FieldVector<double, 2> v = {jacobian[0][0], jacobian[0][1]};
      return v;
    }

    double local(Dune::FieldVector<double, 2> _global) {
      const double tolerance = double(16) * std::numeric_limits<double>::epsilon();

      Dune::FieldVector<double, 1> u{domainMidPoint()};
      Dune::FieldVector<double, 1> du;
      do {
        const Dune::FieldVector<double, 2> dGlobal = (*this)(u)-_global;

        using MatrixHelper = typename Dune::MultiLinearGeometryTraits<double>::MatrixHelper;

        MatrixHelper::template xTRightInvA<gridDim, worldDim>(jacobianTransposed(u), dGlobal, du);

        // Update
        u -= du;

        // Clamp result, this might be already an indicator that there is no solution
        if (Dune::FloatCmp::gt(u[0], patchData->knotSpans[0].back())) u[0] = patchData->knotSpans[0].back();
        if (Dune::FloatCmp::lt(u[0], patchData->knotSpans[0].front())) u[0] = patchData->knotSpans[0].front();

      } while (du.two_norm() > tolerance);

      return u[0];
    }
  };

  class Boundary {
    // Abkürzungen
    using Point       = Dune::FieldVector<double, 2>;
    using PointVector = std::vector<Point>;
    using Domain      = std::array<double, 2>;

   public:
    NurbsCurveHandler nurbsGeometry;

    Domain domain{};

    // Control Points are the start and endpoint (Rename???)
    PointVector controlPoints;

    Boundary() = default;
    Boundary(NurbsCurveHandler& _nurbsGeometry, Domain _domain, PointVector& _controlPoints)
        : nurbsGeometry(_nurbsGeometry), domain(_domain), controlPoints(_controlPoints) {}

    Boundary(NurbsCurveHandler& _nurbsGeometry, Domain _domain) : nurbsGeometry(_nurbsGeometry), domain(_domain) {
      controlPoints = {nurbsGeometry(domain[0]), nurbsGeometry(domain[1])};
    }

    explicit Boundary(Ibra::BrepTrim& _trim) {
      domain        = _trim.domain;
      nurbsGeometry = NurbsCurveHandler(_trim.geometry);
      controlPoints = {nurbsGeometry(domain[0]), nurbsGeometry(domain[1])};
    }

    Boundary(Ibra::BrepTrim& _trim, Domain _domain) : domain(_domain) {
      nurbsGeometry = NurbsCurveHandler(_trim.geometry);
      controlPoints = {nurbsGeometry(domain[0]), nurbsGeometry(domain[1])};
    }

    explicit Boundary(PointVector& _controlPoints) : controlPoints(_controlPoints) {
      assert(controlPoints.size() == 2);
      using ControlPoint = NurbsCurveHandler::ControlPoint;

      std::array<ControlPoint, 2> points{ControlPoint{.p{controlPoints.front()}, .w = 1},
                                         {.p{controlPoints.back()}, .w = 1}};
      nurbsGeometry = NurbsCurveHandler(points);
      domain        = nurbsGeometry.domain();
    }

    [[nodiscard]] int degree() const { return nurbsGeometry.degree(); };

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

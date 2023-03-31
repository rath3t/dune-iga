//
// Created by Henri on 17.03.2023.
//

#pragma once

#include <dune/iga/igaalgorithms.hh>

namespace Dune::IGA {

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl =  DuneLinearAlgebraTraits<double>>
  class NURBSPatchGeometry {
   public:
    static constexpr std::integral auto patchDim = dim;
    static constexpr std::integral auto coorddimension = dimworld;

    using LinearAlgebraTraits = NurbsGridLinearAlgebraTraitsImpl;
    using ctype               = typename LinearAlgebraTraits::value_type;
    using LocalCoordinate     = typename LinearAlgebraTraits::template FixedVectorType<patchDim>;
    using GlobalCoordinate    = typename LinearAlgebraTraits::template FixedVectorType<coorddimension>;
    using JacobianTransposed  = typename LinearAlgebraTraits::template FixedMatrixType<patchDim, coorddimension>;
    using JacobianInverseTransposed =
        typename LinearAlgebraTraits::template FixedMatrixType<coorddimension, patchDim>;

    using ControlPointType = typename NURBSPatchData<patchDim, dimworld, LinearAlgebraTraits>::ControlPointType;

   private:
    /* Helper class to compute a matrix pseudo inverse */
    using MatrixHelper = typename MultiLinearGeometryTraits<ctype>::MatrixHelper;

   public:

    NURBSPatchGeometry() = default;

    explicit NURBSPatchGeometry(std::shared_ptr<NURBSPatchData<patchDim, dimworld, LinearAlgebraTraits>> patchData) : patchData_(patchData) {
      nurbs_ = Nurbs<patchDim, LinearAlgebraTraits>(*patchData);
    }

    /** \brief Map the center of the element to the geometry */
    [[nodiscard]] GlobalCoordinate center() const {
      LocalCoordinate localcenter(0.5);
      return global(localcenter);
    }

    /** \brief Computes the volume of the element with an integration rule for order max(order)*elementdim */
    [[nodiscard]] double volume() const {
      const auto rule = Dune::QuadratureRules<ctype, patchDim>::rule(
          this->type(), patchDim * (*std::ranges::max_element(patchData_->degree)));
      ctype vol = 0.0;
      for (auto& gp : rule)
        vol += integrationElement(gp.position()) * gp.weight();
      return vol;
    }

    /** \brief Return the number of corners of the element */
    [[nodiscard]] int corners() const { return 1 << patchDim; }

    /** \brief Return world coordinates of the k-th corner of the element */
    [[nodiscard]] GlobalCoordinate corner(int k) const {
      LocalCoordinate localcorner;
      for (size_t i = 0; i < patchDim; i++) {
        localcorner[i] = (k & (1 << i)) ? 1 : 0;
      }
      return global(localcorner);
    }

    /** \brief evaluates the geometric position
     *
     *  \param[in] local local coordinates for each dimension in [0,1] domain
     */
    [[nodiscard]] FieldVector<ctype, dimworld> global(const LocalCoordinate& local) const {
      std::array<double, dim> u{};
      std::ranges::copy(local, std::begin(u));

      auto spanIndex = findSpanUncorrected(patchData_->degree, u, patchData_->knotSpans);
      auto cpCoordinateNet = netOfSpan(u, patchData_->knotSpans, patchData_->degree, extractControlCoordinates(patchData_->controlPoints));
      auto basis     = nurbs_.basisFunctionNet(u);

      return Dune::IGA::dot(basis, cpCoordinateNet);
    }

    FieldVector<ctype, dimworld> operator()(const LocalCoordinate& local) const {
      return global(local);
    }

    /** \brief Inverse of global this function gets a point defined in the world space and return
     * the closest point in local coordinates, i.e. in [0,1] domain for each grid dimension
     *
     *  \param global global coordinates for the point where the local coordinates are searched for
     */
    [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
      const ctype tolerance = ctype(16) * std::numeric_limits<ctype>::epsilon();
      LocalCoordinate x{};
      std::ranges::copy(domainMidPoint(), std::begin(x));

      LocalCoordinate dx{};
      do {  // from multilinearGeometry
        const GlobalCoordinate dglobal = (*this).global(x) - global;
        MatrixHelper::template xTRightInvA<patchDim, coorddimension>(jacobianTransposed(x), dglobal, dx);
        const bool invertible
            = MatrixHelper::template xTRightInvA<patchDim, coorddimension>(jacobianTransposed(x), dglobal, dx);

        if (!invertible) return LocalCoordinate(std::numeric_limits<ctype>::max());
        x -= dx;
        // if local is outside the maximum knot vector span bound, thus we clamp it to it and hope for convergence
        for (int i = 0; i < patchDim; ++i) {
          if (Dune::FloatCmp::gt(x[i], patchData_->knotSpans[i].back())) x[i] = patchData_->knotSpans[i].back();
          if (Dune::FloatCmp::lt(x[i], patchData_->knotSpans[i].front())) x[i] = patchData_->knotSpans[i].front();
        }

      } while (dx.two_norm2() > tolerance);
      return x;
    }

    /** \brief compute the Jacobian transposed matrix
     *
     *  \param local local coordinates for each dimension
     */
    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
      std::array<double, dim> u{};
      std::ranges::copy(local, std::begin(u));

      JacobianTransposed result;
      std::array<unsigned int, patchDim> subDirs;
      for (int subI = 0, i = 0; i < patchDim; ++i) {
        // if (fixedOrVaryingDirections_[i] == Impl::FixedOrFree::fixed) continue;
        subDirs[subI++] = i;
      }

      const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(u, 1);
      auto cpCoordinateNet = netOfSpan(u, patchData_->knotSpans, patchData_->degree, extractControlCoordinates(patchData_->controlPoints));

      for (int dir = 0; dir < patchDim; ++dir) {
        std::array<unsigned int, patchDim> ithVec{};
        ithVec[subDirs[dir]] = 1;
        result[dir]          = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet);

      }
      return result;
    }

    [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
      auto j = jacobianTransposed(local);
      return MatrixHelper::template sqrtDetAAT<patchDim, coorddimension>(j);
    }

    /** \brief Type of the element: a hypercube of the correct dimension */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(patchDim); }

    // The following are functions are not part of the IbraBase Interface
    [[nodiscard]] std::array<std::array<double, 2>, patchDim> domain() const {
      std::array<std::array<double, 2>, patchDim> result{};
      for (int i = 0; i < patchDim; ++i) {
        auto [min, max] = std::ranges::minmax_element(patchData_->knotSpans[i]);
        result[i] = {*min, *max};
      }

      return result;
    }

    [[nodiscard]] std::array<int, patchDim> degree() const {
      return patchData_->degree;
    }

    [[nodiscard]] std::array<ctype, patchDim> domainMidPoint() const {
      auto dom = domain();
      std::array<ctype, patchDim> result{};
      for (int i = 0; i < patchDim; ++i)
        result[i] = (dom[i][0] + dom[i][1]) / 2;

      return result;
    }
   public:
    std::shared_ptr<NURBSPatchData<patchDim, dimworld, LinearAlgebraTraits>> patchData_;

   private:
    Dune::IGA::Nurbs<patchDim, LinearAlgebraTraits> nurbs_;

  };

}



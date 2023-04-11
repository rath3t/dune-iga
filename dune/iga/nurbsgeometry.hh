// SPDX-FileCopyrightText: 2022 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

//
// Created by lex on 16.11.21.
//

#pragma once
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/trimmedElementRepresentation.hh>

namespace Dune::IGA {
  namespace Impl {
    enum class FixedOrFree { fixed, free };
  }

  /** \brief a geometry implementation for NURBS*/
  template <std::integral auto mydim, std::integral auto dimworld, class GridImpl>
  class NURBSGeometry {
   public:
    static constexpr std::integral auto mydimension = mydim;

    static constexpr std::integral auto coorddimension = dimworld;
    static constexpr std::integral auto griddim        = GridImpl::dimension;

    using ctype               = typename GridImpl::LinearAlgebraTraits::value_type;
    using LinearAlgebraTraits = typename GridImpl::LinearAlgebraTraits;
    using LocalCoordinate     = typename LinearAlgebraTraits::template FixedVectorType<mydimension>;
    using GlobalCoordinate    = typename LinearAlgebraTraits::template FixedVectorType<coorddimension>;
    using JacobianTransposed  = typename LinearAlgebraTraits::template FixedMatrixType<mydimension, coorddimension>;
    using JacobianInverseTransposed =
        typename LinearAlgebraTraits::template FixedMatrixType<coorddimension, mydimension>;

    using ControlPointType = typename NURBSPatchData<griddim, dimworld, LinearAlgebraTraits>::ControlPointType;

    using TrimmedElementRepresentationType = GridImpl::Traits::TrimmedElementRepresentationType;

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
    NURBSGeometry(std::shared_ptr<NURBSPatchData<griddim, dimworld, LinearAlgebraTraits>> patchData,
                  const std::array<Impl::FixedOrFree, griddim>& fixedOrVaryingDirections,
                  const std::array<int, griddim>& thisSpanIndices,
                  const std::optional<std::shared_ptr<TrimmedElementRepresentationType>>& trimRepr = std::nullopt)
        : patchData_(patchData), fixedOrVaryingDirections_{fixedOrVaryingDirections}, trimRepr_(trimRepr) {
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

    NURBSGeometry() = default;

    /** \brief Map the center of the element to the geometry */
    [[nodiscard]] GlobalCoordinate center() const {
      LocalCoordinate localcenter(0.5);
      return global(localcenter);
    }

    // TODO Is this the volume in the parameter Space or in the physical space

    /** \brief Computes the volume of the element with an integration rule for order max(order)*elementdim */
    [[nodiscard]] double volume() const {
      if (trimRepr_.has_value())
        return trimRepr_.value()->calculateArea();


      const auto rule = Dune::QuadratureRules<ctype, mydimension>::rule(
          this->type(), mydimension * (*std::ranges::max_element(patchData_->degree)));
      // TODO Maybe use std::accumulate
      ctype vol = 0.0;
      for (auto& gp : rule)
        vol += integrationElement(gp.position()) * gp.weight();
      return vol;
    }

    /** \brief Return the number of corners of the element */
    [[nodiscard]] int corners() const { return 1 << mydimension; }

    /** \brief Return world coordinates of the k-th corner of the element */
    [[nodiscard]] GlobalCoordinate corner(int k) const {
      LocalCoordinate localcorner;
      for (size_t i = 0; i < mydimension; i++) {
        localcorner[i] = (k & (1 << i)) ? 1 : 0;
      }
      return global(localcorner);
    }

    [[nodiscard]] bool affine() const { return false; }

    /** \brief evaluates the geometric position
     *
     *  \param[in] local local coordinates for each dimension in [0,1] domain
     */
    [[nodiscard]] GlobalCoordinate global(const LocalCoordinate& local) const {
      const auto localInSpan = localToSpan(local);

      auto basis = nurbs_.basisFunctionNet(localInSpan);
      return dot(basis, cpCoordinateNet_);
    }

    /** \brief Inverse of global this function gets a point defined in the world space and return
     * the closest point in local coordinates, i.e. in [0,1] domain for each grid dimension
     *
     *  \param global global coordinates for the point where the local coordinates are searched for
     */
    [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
      const ctype tolerance = ctype(16) * std::numeric_limits<ctype>::epsilon();
      LocalCoordinate x     = LocalCoordinate(0.5);
      LocalCoordinate dx{};
      do {  // from multilinearGeometry
        const GlobalCoordinate dglobal = (*this).global(x) - global;
        MatrixHelper::template xTRightInvA<mydimension, coorddimension>(jacobianTransposed(x), dglobal, dx);
        const bool invertible
            = MatrixHelper::template xTRightInvA<mydimension, coorddimension>(jacobianTransposed(x), dglobal, dx);

        if (!invertible) return LocalCoordinate(std::numeric_limits<ctype>::max());
        x -= dx;
        // if local is outside the maximum knot vector span bound, thus we clamp it to it and hope for convergence
        for (int i = 0; i < mydim; ++i)
          if (Dune::FloatCmp::gt(x[i], patchData_->knotSpans[i].back())) x[i] = patchData_->knotSpans[i].back();

      } while (dx.two_norm2() > tolerance);
      return x;
    }

    /** \brief compute the Jacobian transposed matrix
     *
     *  \param local local coordinates for each dimension
     */
    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
      JacobianTransposed result;
      std::array<unsigned int, mydimension> subDirs;
      for (int subI = 0, i = 0; i < griddim; ++i) {
        if (fixedOrVaryingDirections_[i] == Impl::FixedOrFree::fixed) continue;
        subDirs[subI++] = i;
      }

      const auto localInSpan              = localToSpan(local);
      const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(localInSpan, 1);

      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 1;
        result[dir]          = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet_);
        result[dir] *= scaling_[subDirs[dir]];  // transform back to 0..1 domain
      }
      return result;
    }

    [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
      auto j = jacobianTransposed(local);
      return MatrixHelper::template sqrtDetAAT<mydimension, coorddimension>(j);
    }

    [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
      JacobianInverseTransposed jacobianInverseTransposed1;
      MatrixHelper::template rightInvA<mydimension, coorddimension>(jacobianTransposed(local),
                                                                    jacobianInverseTransposed1);
      return jacobianInverseTransposed1;
    }

    [[nodiscard]] GlobalCoordinate unitNormal(const LocalCoordinate& local) const requires(mydimension == 2)
        && (coorddimension == 3) {
      auto N = normal(local);
      return N / N.two_norm();
    }

    [[nodiscard]] GlobalCoordinate normal(const LocalCoordinate& local) const requires(mydimension == 2)
        && (coorddimension == 3) {
      auto J = jacobianTransposed(local);
      return cross(J[0], J[1]);
    }

    [[nodiscard]] ctype gaussianCurvature(const LocalCoordinate& local) const requires(mydimension == 2) {
      auto metricDet = metric(local).determinant();
      auto secondF   = secondFundamentalForm(local).determinant();
      return secondF / metricDet;
    }

    auto metric(const LocalCoordinate& local) const {
      const auto J = jacobianTransposed(local);
      FieldMatrix<ctype, mydimension, mydimension> metric;
      MatrixHelper::AAT(J, metric);
      return metric;
    }

    auto secondFundamentalForm(const LocalCoordinate& local) const requires(mydimension == 2) && (coorddimension == 3) {
      const auto secDerivatives = secondDerivativeOfPosition(local);
      const auto unitnormal     = unitNormal(local);
      FieldMatrix<ctype, mydimension, mydimension> b;
      b[0][0] = secDerivatives[0] * unitnormal;
      b[1][1] = secDerivatives[1] * unitnormal;
      b[0][1] = b[1][0] = secDerivatives[2] * unitnormal;
      return b;
    }

    auto secondDerivativeOfPosition(const LocalCoordinate& local) const {
      FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, coorddimension> result;
      std::array<unsigned int, mydimension> subDirs;
      for (int subI = 0, i = 0; i < griddim; ++i) {
        if (fixedOrVaryingDirections_[i] == Dune::IGA::Impl::FixedOrFree::fixed) continue;
        subDirs[subI++] = i;
      }
      const auto localInSpan              = localToSpan(local);
      const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(localInSpan, 2);
      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 2;  // second derivative in dir direction
        result[dir]          = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet_);
        result[dir] *= Dune::power(scaling_[subDirs[dir]], 2);  // transform back to 0..1
      }
      std::array<int, mydimension> mixeDerivs;
      std::ranges::fill_n(mixeDerivs.begin(), 2, 1);  // first mixed derivatives
      int mixedDireCounter = mydimension;
      do {
        result[mixedDireCounter++] = dot(basisFunctionDerivatives.get(mixeDerivs), cpCoordinateNet_);
        for (int dir = 0; dir < mixeDerivs.size(); ++dir) {
          if (mixeDerivs[dir] == 0) continue;
          result[mixedDireCounter - 1] *= scaling_[subDirs[dir]];
        }
      } while (std::ranges::next_permutation(mixeDerivs, std::greater()).found);

      return result;
    }

    /** \brief Type of the element: a hypercube of the correct dimension */
    [[nodiscard]] GeometryType type() const {
      if (trimRepr_.has_value())
        return GeometryTypes::none(mydimension);
      else
        return GeometryTypes::cube(mydimension);
    }

    /** \brief Return from the 0 to 1 domain the position in the current knot span */
    template <typename ReturnType = std::array<typename LocalCoordinate::value_type, griddim>>
    auto localToSpan(const LocalCoordinate& local) const {
      ReturnType localInSpan;
      if constexpr (LocalCoordinate::dimension != 0) {
        for (int loci = 0, i = 0; i < griddim; ++i) {
          localInSpan[i] = (fixedOrVaryingDirections_[i] == Impl::FixedOrFree::free)
                               ? local[loci++] * scaling_[i] + offset_[i]
                               : offset_[i];
        }
      } else
        for (int i = 0; i < griddim; ++i)
          localInSpan[i] = offset_[i];
      return localInSpan;
    }

    [[nodiscard]] const auto& nurbs() const { return nurbs_; }

    [[nodiscard]] const auto& controlPoints() const { return cpCoordinateNet_; }

    std::shared_ptr<NURBSPatchData<griddim, dimworld, LinearAlgebraTraits>> patchData_;
    std::array<int, griddim> thisSpanIndices_;
    std::array<Impl::FixedOrFree, griddim> fixedOrVaryingDirections_{Impl::FixedOrFree::free};
    Dune::IGA::Nurbs<griddim, LinearAlgebraTraits> nurbs_;
    std::array<ctype, griddim> offset_;
    std::array<ctype, griddim> scaling_;
    MultiDimensionNet<griddim, typename ControlPointType::VectorType> cpCoordinateNet_;
    std::optional<std::shared_ptr<TrimmedElementRepresentationType>> trimRepr_;
  };

  template <std::integral auto mydim, std::integral auto dimworld, class GridImpl>
  auto referenceElement(const NURBSGeometry<mydim, dimworld, GridImpl>& geo) {
    return Dune::ReferenceElements<typename GridImpl::LinearAlgebraTraits::value_type, mydim>::cube();
  };
}  // namespace Dune::IGA

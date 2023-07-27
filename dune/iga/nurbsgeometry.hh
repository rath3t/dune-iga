// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "dune/iga/geometry/closestpointprojection.hh"
#include "dune/iga/nurbsalgorithms.hh"
#include "dune/iga/trim/subgrid.hh"
#include "dune/iga/utils/fillquadraturerule.hh"
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/geometry.hh>

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

    using ctype                     = typename GridImpl::ctype;
    using LocalCoordinate           = Dune::FieldVector<ctype, mydimension>;
    using GlobalCoordinate          = Dune::FieldVector<ctype, coorddimension>;
    using JacobianTransposed        = Dune::FieldMatrix<ctype, mydimension, coorddimension>;
    using Jacobian                  = Dune::FieldMatrix<ctype, coorddimension, mydimension>;
    using JacobianInverseTransposed = Dune::FieldMatrix<ctype, coorddimension, mydimension>;
    using JacobianInverse           = Dune::FieldMatrix<ctype, mydimension, coorddimension>;

    using ControlPointType = typename NURBSPatchData<griddim, dimworld, ctype>::ControlPointType;

    using SubGridType = GridImpl::Traits::SubGridType;

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
    NURBSGeometry(std::shared_ptr<NURBSPatchData<griddim, dimworld, ctype>> patchData,
                  const std::array<Impl::FixedOrFree, griddim>& fixedOrVaryingDirections,
                  const std::array<int, griddim>& thisSpanIndices, const std::shared_ptr<SubGridType> subGrid = nullptr)
        : patchData_(patchData), fixedOrVaryingDirections_{fixedOrVaryingDirections}, subgrid_(subGrid) {
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

      nurbs_ = Dune::IGA::Nurbs<griddim, ctype>(*patchData, thisSpanIndices_);
      cpCoordinateNet_
          = netOfSpan(thisSpanIndices_, patchData_->degree, extractControlCoordinates(patchData_->controlPoints));
    }

    NURBSGeometry() = default;

    /** \brief Map the center of the element to the geometry */
    [[nodiscard]] GlobalCoordinate center() const {
      LocalCoordinate localcenter(0.5);
      return global(localcenter);
    }

    /** \brief Computes the volume of the element with an integration rule for order max(order)*elementdim */
    [[nodiscard]] double volume() const {
      if constexpr (mydimension == 2)
        if (subgrid_) {
          Dune::QuadratureRule<double, mydimension> rule;
          fillQuadratureRuleImpl(rule, *subgrid_.get(), (*std::ranges::max_element(patchData_->degree)));
          ctype vol = 0.0;
          for (auto& gp : rule)
            vol += integrationElement(gp.position()) * gp.weight();
          return vol;
        }

      const auto& rule = Dune::QuadratureRules<ctype, mydimension>::rule(
          GeometryTypes::cube(mydimension), (*std::ranges::max_element(patchData_->degree)));
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
      if constexpr (mydimension == 0)
        return {};
      else if constexpr (mydimension != coorddimension) {
        auto [u, Ru, fu, gap] = Dune::IGA::closestPointProjectionByTrustRegion(*this, global);
        return u;
      } else {
        LocalCoordinate x = LocalCoordinate(0.5);
        LocalCoordinate dx{};

        GlobalCoordinate dglobal;
        do {  // from multilinearGeometry
          dglobal = (*this).global(x) - global;
          MatrixHelper::template xTRightInvA<mydimension, coorddimension>(jacobianTransposed(x), dglobal, dx);
          const bool invertible
              = MatrixHelper::template xTRightInvA<mydimension, coorddimension>(jacobianTransposed(x), dglobal, dx);

          if (!invertible) return LocalCoordinate(std::numeric_limits<ctype>::max());
          x -= dx;
          // if local is outside of maximum knot vector span bound, thus we clamp it to it and return
          // clamp result into boundaries
          if (Utilities::clampToBoundaryAndCheckIfIsAtAllBoundaries(x, domain())) {
            break;
          }

        } while (dx.two_norm2() > tolerance);

        return x;
      }
    }

    /** \brief compute the Jacobian transposed matrix
     *
     *  \param local local coordinates for each dimension
     */
    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
      JacobianTransposed result;
      std::array<unsigned int, mydimension> subDirs = getDirectionsOfSubEntityInParameterSpace();

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

    [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const { return transpose(jacobianTransposed(local)); }

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

    [[nodiscard]] JacobianInverse jacobianInverse(const LocalCoordinate& local) const {
      return transpose(jacobianInverseTransposed(local));
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
      std::array<unsigned int, mydimension> subDirs = getDirectionsOfSubEntityInParameterSpace();
      const auto localInSpan                        = localToSpan(local);
      const auto basisFunctionDerivatives           = nurbs_.basisFunctionDerivatives(localInSpan, 2);
      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 2;  // second derivative in dir direction
        result[dir]          = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet_);
        result[dir] *= Dune::power(scaling_[subDirs[dir]], 2);  // transform back to 0..1
      }
      if constexpr (mydimension > 1 and griddim > 1) {
        std::array<int, griddim> mixeDerivs;
        std::ranges::fill(mixeDerivs, 0);  // first mixed derivatives
        if constexpr (mydimension == 2)
          for (int dir = 0; dir < mydimension; ++dir) {
            mixeDerivs[subDirs[dir]] = 1;
          }
        else
          std::ranges::fill_n(mixeDerivs.begin(), 2, 1);  // first mixed derivatives
        int mixedDireCounter = mydimension;
        do {
          result[mixedDireCounter++] = dot(basisFunctionDerivatives.get(mixeDerivs), cpCoordinateNet_);
          for (int dir = 0; dir < mydimension; ++dir) {
            if (mixeDerivs[dir] == 0) continue;
            result[mixedDireCounter - 1] *= scaling_.at(subDirs[dir]);
          }
          if constexpr (mydimension == 2) break;

        } while (std::ranges::next_permutation(mixeDerivs, std::greater()).found);
      }

      return result;
    }

    auto getDirectionsOfSubEntityInParameterSpace() const {
      std::array<unsigned int, mydimension> subDirs;
      for (int subI = 0, i = 0; i < griddim; ++i) {
        if constexpr (mydimension != griddim)
          if (fixedOrVaryingDirections_[i] == Dune::IGA::Impl::FixedOrFree::fixed) continue;
        subDirs[subI++] = i;
      }
      return subDirs;
    }
    auto zeroFirstAndSecondDerivativeOfPosition(const LocalCoordinate& uL) const {
      const auto u = localToSpan(uL);
      FieldVector<ctype, dimworld> pos;

      FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, coorddimension> H;
      std::array<unsigned int, mydimension> subDirs = getDirectionsOfSubEntityInParameterSpace();

      const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(u, 2);

      std::array<unsigned int, griddim> ithVecZero{};
      pos = Dune::IGA::dot(basisFunctionDerivatives.get(ithVecZero), cpCoordinateNet_);
      JacobianTransposed J;
      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 1;
        J[dir]               = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet_);
        J[dir] *= scaling_.at(subDirs[dir]);  // transform back to 0..1 domain
      }

      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 2;  // second derivative in subDirs[dir] direction
        H[dir]               = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet_);

        H[dir] *= Dune::power(scaling_.at(subDirs[dir]), 2);  // transform back to 0..1
      }
      if constexpr (mydimension > 1 and griddim > 1) {
        std::array<int, griddim> mixeDerivs;
        std::ranges::fill(mixeDerivs, 0);  // first mixed derivatives
        if constexpr (mydimension == 2)
          for (int dir = 0; dir < mydimension; ++dir) {
            mixeDerivs[subDirs[dir]] = 1;
          }
        else
          std::ranges::fill_n(mixeDerivs.begin(), 2, 1);  // first mixed derivatives
        int mixedDireCounter = mydimension;
        do {
          H[mixedDireCounter++] = dot(basisFunctionDerivatives.get(mixeDerivs), cpCoordinateNet_);
          for (int dir = 0; dir < mydimension; ++dir) {
            if (mixeDerivs[dir] == 0) continue;
            H[mixedDireCounter - 1] *= scaling_.at(subDirs[dir]);
          }
          if constexpr (mydimension == 2) break;

        } while (std::ranges::next_permutation(mixeDerivs, std::greater()).found);
      }

      return std::make_tuple(pos, J, H);
    }

    /** \brief Type of the element: a hypercube of the correct dimension */
    [[nodiscard]] GeometryType type() const {
      if (subgrid_)
        return GeometryTypes::none(mydimension);
      else
        return GeometryTypes::cube(mydimension);
    }

    /** \brief Return from the 0 to 1 domain the position in the current knot span */
    template <typename ReturnType = Dune::FieldVector<typename LocalCoordinate::value_type, griddim>>
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

    /** \brief Return from the 0 to 1 domain the position in the current knot span */
    template <typename ReturnType = LocalCoordinate>
    auto spanToLocal(const LocalCoordinate& inSpan) const {
      ReturnType localL;
      for (int i = 0; i < griddim; ++i) {
        localL[i] = (inSpan[i] - offset_[i]) / scaling_[i];
        localL[i] = clampToDomain(localL[i], domain());
      }

      return localL;
    }

    // The following are functions are not part of the Geometry Interface
    [[nodiscard]] std::array<Utilities::Domain<double>, mydimension> domain() const { return {}; }

    [[nodiscard]] std::array<int, griddim> degree() const { return patchData_->degree; }

    [[nodiscard]] Dune::FieldVector<ctype, mydimension> domainMidPoint() const {
      auto dom = domain();
      Dune::FieldVector<ctype, mydimension> result{};
      for (int i = 0; i < mydimension; ++i)
        result[i] = dom[i].center();

      return result;
    }

    [[nodiscard]] const auto& nurbs() const { return nurbs_; }

    [[nodiscard]] const auto& controlPoints() const { return cpCoordinateNet_; }

    std::shared_ptr<NURBSPatchData<griddim, dimworld, ctype>> patchData_;
    std::array<int, griddim> thisSpanIndices_;
    std::array<Impl::FixedOrFree, griddim> fixedOrVaryingDirections_{Impl::FixedOrFree::free};
    Dune::IGA::Nurbs<griddim, ctype> nurbs_;
    std::array<ctype, griddim> offset_;
    std::array<ctype, griddim> scaling_;
    MultiDimensionNet<griddim, typename ControlPointType::VectorType> cpCoordinateNet_;
    std::shared_ptr<SubGridType> subgrid_;
  };

  template <std::integral auto mydim, std::integral auto dimworld, class GridImpl>
  auto referenceElement(const NURBSGeometry<mydim, dimworld, GridImpl>& geo) {
    return Dune::ReferenceElements<typename GridImpl::ctype, mydim>::cube();
  };
}  // namespace Dune::IGA

// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/geometry/closestpointprojection.hh"
#include "dune/iga/geometry/geohelper.hh"
#include "dune/iga/nurbsalgorithms.hh"
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/hierarchicpatch/nurbspatchgeometrylocalview.hh>
#include <dune/iga/hierarchicpatch/geometryalgorithms.hh>

namespace Dune::IGANEW::GeometryKernel {

  template <int dim_, int dimworld_, typename ScalarType = double>
  class NURBSPatchGeometry {
   public:
    static constexpr std::integral auto worlddimension = dimworld_;
    static constexpr std::integral auto mydimension = dim_;

    using ctype                     = ScalarType;
    using LocalCoordinate           = FieldVector<ctype, mydimension>;
    using GlobalCoordinate          = FieldVector<ctype, worlddimension>;
    using JacobianTransposed        = FieldMatrix<ctype, mydimension, worlddimension>;
    using Hessian                   = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, worlddimension>;
    using Jacobian                  = FieldMatrix<ctype, worlddimension, mydimension>;
    using JacobianInverseTransposed = FieldMatrix<ctype, worlddimension, mydimension>;
    using JacobianInverse           = FieldMatrix<ctype, mydimension, worlddimension>;
    using Volume                    = ctype;

    template <int codim>
    using ParameterSpaceGeometry
        = YaspGeometry<mydimension - codim, mydimension,
                       const YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>>;

    using ControlPointType = typename IGA::NURBSPatchData<mydimension, worlddimension, ScalarType>::ControlPointType;
    using ControlPointNetType =
        typename IGA::NURBSPatchData<mydimension, worlddimension, ScalarType>::ControlPointNetType;
    using ControlPointCoordinateNetType = IGA::MultiDimensionNet<
        mydimension, typename IGA::NURBSPatchData<mydimension, worlddimension, ScalarType>::GlobalCoordinateType>;
    using Nurbs          = Dune::IGA::Nurbs<mydimension, ScalarType>;
    using NurbsLocalView = typename Nurbs::LocalView;
    template <int codim, Trimming trim>
    using GeometryLocalView = PatchGeometryLocalView<codim, GeometryKernel::NURBSPatch, trim>;

    template <int codim, typename NURBSPatchGeometry, Trimming trim>
    friend struct PatchGeometryLocalView;

   private:
    /* Helper class to compute a matrix pseudo inverse */
    using MatrixHelper = typename MultiLinearGeometryTraits<ctype>::MatrixHelper;

   public:
    NURBSPatchGeometry() = default;

    template <int codim, Trimming trim>
    auto localView() const {
      return GeometryLocalView<codim, trim>(*this);
    }

    explicit NURBSPatchGeometry(const IGA::NURBSPatchData<mydimension, worlddimension, ScalarType>& patchData)
        : patchData_(patchData), nurbs_{patchData}, localView_(localView<0, Trimming::Disabled>()) {}

    /** \brief Map the center of the element to the geometry */
    [[nodiscard]] GlobalCoordinate center() const {
      LocalCoordinate localcenter(0.5);
      return global(localcenter);
    }

    [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
      return computeParameterSpaceCoordinate(*this, global);
    }

    /** \brief Computes the volume of the element with an integration rule for order max(order)*elementdim */
    [[nodiscard]] double volume(int scaleOrder = 1) const {
      const auto rule = Dune::QuadratureRules<ctype, mydimension>::rule(
          this->type(), scaleOrder * mydimension * (*std::ranges::max_element(patchData_.degree)));
      ctype vol = 0.0;
      for (auto& gp : rule) {
        localView_.bind(gp.position());
        vol += localView_.integrationElement(gp.position()) * gp.weight();
      }
      return vol;
    }

    [[nodiscard]] bool affine() const { return false; }

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

    FieldVector<ctype, worlddimension> operator()(const LocalCoordinate& local) const { return global(local); }

    /** \brief evaluates the geometric position
     *
     *  \param[in] u coordinates for each dimension in the [knotSpan.front(), knotSpan.back() ] domain
     */
    [[nodiscard]] FieldVector<ctype, worlddimension> global(const LocalCoordinate& u) const {
      auto [nurbsLocalView, cpNet, subNetStart] = calculateNurbsAndControlPointNet(u);
      return this->global(u, nurbsLocalView, cpNet);
    }

    auto zeroFirstAndSecondDerivativeOfPosition(const LocalCoordinate& u) const {
      auto [nurbsLocalView, cpNet, subNetStart] = calculateNurbsAndControlPointNet(u);
      return zeroFirstAndSecondDerivativeOfPositionImpl(u, nurbsLocalView, cpNet);
    }

    /** \brief compute the Jacobian transposed matrix
     *
     *  \param local coordinates for each dimension
     */
    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& u) const {
      auto [nurbsLocalView, cpNet, subNetStart] = calculateNurbsAndControlPointNet(u);
      return jacobianTransposedImpl(u, nurbsLocalView, cpNet);
    }

    [[nodiscard]] Hessian hessian(const LocalCoordinate& u) const {
      auto [nurbsLocalView, cpNet, subNetStart] = calculateNurbsAndControlPointNet(u);
      return GeometryKernel::hessian<GeometryKernel::NURBSPatch>(u, nurbsLocalView, cpNet);
    }

    [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const { return transpose(jacobianTransposed(local)); }

    [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
      auto j = jacobianTransposed(local);
      return MatrixHelper::template sqrtDetAAT<mydimension, worlddimension>(j);
    }

    /** \brief Type of the element: a hypercube of the correct dimension */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }

    // The following are functions are not part of the Geometry Interface
    [[nodiscard]] std::array<IGA::Utilities::Domain<double>, mydimension> domain() const {
      std::array<IGA::Utilities::Domain<double>, mydimension> result{};
      for (int i = 0; i < mydimension; ++i)
        result[i] = {patchData_.knotSpans[i].front(), patchData_.knotSpans[i].back()};

      return result;
    }

    [[nodiscard]] std::array<int, mydimension> degree() const { return patchData_.degree; }

    [[nodiscard]] FieldVector<ctype, mydimension> domainMidPoint() const {
      auto dom = domain();
      Dune::FieldVector<ctype, mydimension> result{};
      for (int i = 0; i < mydimension; ++i)
        result[i] = dom[i].center();

      return result;
    }

   private:
    auto calculateNurbsAndControlPointNet(const LocalCoordinate& u) const {
      auto subNetStart = IGA::findSpanCorrected(patchData_.degree, u, patchData_.knotSpans);

      auto cpCoordinateNet
          = netOfSpan(subNetStart, patchData_.degree, IGA::extractControlCoordinates(patchData_.controlPoints));
      auto nurbsLocalView = nurbs_.localView();
      nurbsLocalView.bind(subNetStart);
      return std::make_tuple(nurbsLocalView, cpCoordinateNet, subNetStart);
    }

    [[nodiscard]] static FieldVector<ctype, worlddimension> globalImpl(
        const LocalCoordinate& u, const NurbsLocalView& nurbsLocalView,
        const ControlPointCoordinateNetType& localControlPointNet) {
      auto basis = nurbsLocalView.basisFunctions(u);

      return IGA::dot(basis, localControlPointNet);
    }

    [[nodiscard]] static JacobianTransposed jacobianTransposedImpl(
        const LocalCoordinate& u, const NurbsLocalView& nurbsLocalView,
        const ControlPointCoordinateNetType& localControlPointNet) {
      JacobianTransposed result;

      const auto basisFunctionDerivatives = nurbsLocalView.basisFunctionDerivatives(u, 1);

      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, mydimension> ithVec{};
        ithVec[dir] = 1;
        result[dir] = dot(basisFunctionDerivatives.get(ithVec), localControlPointNet);
      }
      return result;
    }


    static auto zeroFirstAndSecondDerivativeOfPositionImpl(const LocalCoordinate& u,
                                                           const NurbsLocalView& nurbsLocalView,
                                                           const ControlPointCoordinateNetType& localControlPointNet) {
      // GlobalCoordinate pos;
      //
      // Hessian H;
      // const auto basisFunctionDerivatives = nurbsLocalView.basisFunctionDerivatives(u, 2);
      //
      // std::array<unsigned int, mydimension> ithVecZero{};
      // pos = Dune::IGA::dot(basisFunctionDerivatives.get(ithVecZero), localControlPointNet);
      // JacobianTransposed J;
      // for (int dir = 0; dir < mydimension; ++dir) {
      //   std::array<unsigned int, mydimension> ithVec{};
      //   ithVec[dir] = 1;
      //   J[dir]               = dot(basisFunctionDerivatives.get(ithVec), localControlPointNet);
      // }
      //
      // for (int dir = 0; dir < mydimension; ++dir) {
      //   std::array<unsigned int, mydimension> ithVec{};
      //   ithVec[dir] = 2;  // second derivative in subDirs[dir] direction
      //   H[dir]               = dot(basisFunctionDerivatives.get(ithVec), localControlPointNet);
      // }
      // if constexpr (mydimension > 1) {
      //   std::array<int, mydimension> mixeDerivs;
      //   std::ranges::fill(mixeDerivs, 0);  // first mixed derivatives
      //   if constexpr (mydimension == 2)
      //     for (int dir = 0; dir < mydimension; ++dir) {
      //       mixeDerivs[dir] = 1;
      //     }
      //   else
      //     std::ranges::fill_n(mixeDerivs.begin()+1, 2, 1);  // first mixed derivatives
      //   int mixedDireCounter = mydimension;
      //   do {
      //     H[mixedDireCounter++] = dot(basisFunctionDerivatives.get(mixeDerivs), localControlPointNet);
      //
      //     if constexpr (mydimension == 2) break;
      //
      //   } while (std::ranges::next_permutation(mixeDerivs, std::less()).found);
      // }

      // TODO the above code is more efficient since there is only on call to the derivatives

      return std::make_tuple(globalImpl(u,nurbsLocalView,localControlPointNet), jacobianTransposedImpl(u,nurbsLocalView,localControlPointNet), GeometryKernel::hessian<GeometryKernel::NURBSPatch>(u,nurbsLocalView,localControlPointNet));
    }

   public:
    IGA::NURBSPatchData<mydimension, worlddimension, ScalarType> patchData_;

   private:
    IGA::Nurbs<mydimension, ScalarType> nurbs_;
    GeometryLocalView<0, Trimming::Disabled> localView_;
  };

}  // namespace Dune::IGANEW

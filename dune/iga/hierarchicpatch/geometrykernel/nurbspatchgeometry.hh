// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/hierarchicpatch/geometrykernel/geohelper.hh>
#include <dune/iga/hierarchicpatch/geometrykernel/geometryalgorithms.hh>
#include <dune/iga/hierarchicpatch/geometrykernel/nurbspatchgeometrylocalview.hh>

namespace Dune::IGANEW::GeometryKernel {

  template <int dim_, int dimworld_, typename ScalarType = double>
  class NURBSPatch {
   public:
    static constexpr std::integral auto worlddimension = dimworld_;
    static constexpr std::integral auto mydimension    = dim_;

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

    using ControlPointType    = typename NURBSPatchData<mydimension, worlddimension, ScalarType>::ControlPointType;
    using ControlPointNetType = typename NURBSPatchData<mydimension, worlddimension, ScalarType>::ControlPointNetType;
    using ControlPointCoordinateNetType
        = MultiDimensionNet<mydimension,
                            typename NURBSPatchData<mydimension, worlddimension, ScalarType>::GlobalCoordinateType>;
    using Nurbs          = Nurbs<mydimension, ScalarType>;
    using NurbsLocalView = typename Nurbs::LocalView;
    template <int codim, Trimming trim>
    using GeometryLocalView = PatchGeometryLocalView<codim, NURBSPatch, trim>;

    template <int codim, typename NURBSPatch, Trimming trim>
    friend struct PatchGeometryLocalView;

   private:
    /* Helper class to compute a matrix pseudo inverse */
    using MatrixHelper = typename MultiLinearGeometryTraits<ctype>::MatrixHelper;

   public:
    NURBSPatch() = default;

    template <int codim, Trimming trim>
    auto localView() const {
      return GeometryLocalView<codim, trim>(*this);
    }

    explicit NURBSPatch(const NURBSPatchData<mydimension, worlddimension, ScalarType>& patchData)
        : patchData_(patchData), nurbs_{patchData}, localView_(localView<0, Trimming::Disabled>()) {}

    /** \brief Map the center of the element to the geometry */
    [[nodiscard]] GlobalCoordinate center() const { return global(domainMidPoint()); }

    [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
      return computeParameterSpaceCoordinate(*this, global);
    }

    /** \brief Computes the volume of the element with an integration rule for order max(order)*elementdim */
    [[nodiscard]] double volume(int scaleOrder = 1) const {
      const auto rule = Dune::QuadratureRules<ctype, mydimension>::rule(
          this->type(), scaleOrder * mydimension * (*std::ranges::max_element(patchData_.degree)));
      ctype vol = 0.0;
      for (auto& gp : rule) {
        vol += integrationElement(gp.position()) * gp.weight();
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
      return GeometryKernel::position(u, nurbsLocalView, cpNet);
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
      return GeometryKernel::jacobianTransposed(u, nurbsLocalView, cpNet);
    }

    [[nodiscard]] Hessian hessian(const LocalCoordinate& u) const {
      auto [nurbsLocalView, cpNet, subNetStart] = calculateNurbsAndControlPointNet(u);
      return GeometryKernel::hessian(u, nurbsLocalView, cpNet);
    }

    [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const { return transpose(jacobianTransposed(local)); }

    [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
      auto j = jacobianTransposed(local);
      return MatrixHelper::template sqrtDetAAT<mydimension, worlddimension>(j);
    }

    /** \brief Type of the element: a hypercube of the correct dimension */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }

    // The following are functions are not part of the Geometry Interface
    [[nodiscard]] std::array<Utilities::Domain<double>, mydimension> domain() const {
      std::array<Utilities::Domain<double>, mydimension> result{};
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
      auto subNetStart = findSpan(patchData_.degree, u, patchData_.knotSpans);

      auto cpCoordinateNet
          = netOfSpan(subNetStart, patchData_.degree, extractControlCoordinates(patchData_.controlPoints));
      auto nurbsLocalView = nurbs_.localView();
      nurbsLocalView.bind(subNetStart);
      return std::make_tuple(nurbsLocalView, cpCoordinateNet, subNetStart);
    }

    static auto zeroFirstAndSecondDerivativeOfPositionImpl(const LocalCoordinate& u,
                                                           const NurbsLocalView& nurbsLocalView,
                                                           const ControlPointCoordinateNetType& localControlPointNet) {
      // TODO the above code is more efficient since there is only on call to the derivatives

      return std::make_tuple(GeometryKernel::position(u, nurbsLocalView, localControlPointNet),
                             GeometryKernel::jacobianTransposed(u, nurbsLocalView, localControlPointNet),
                             GeometryKernel::hessian(u, nurbsLocalView, localControlPointNet));
    }

   public:
    NURBSPatchData<mydimension, worlddimension, ScalarType> patchData_;

   private:
    Nurbs nurbs_;
    GeometryLocalView<0, Trimming::Disabled> localView_;
  };

}  // namespace Dune::IGANEW::GeometryKernel

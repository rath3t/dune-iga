// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/geometry/closestpointprojection.hh"
#include "dune/iga/geometry/geohelper.hh"
#include "dune/iga/nurbsalgorithms.hh"
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/iga/hierarchicpatch/nurbspatchgeometrylocalview.hh>

namespace Dune::IGANEW {

  template <int dim, int dimworld,bool trim=false, typename ScalarType = double>
  class NURBSPatchGeometry {
   public:

    static constexpr std::integral auto coorddimension = dimworld;
    static constexpr std::integral auto patchDimension = dim;

    using ctype                     = ScalarType;
    using LocalCoordinate           = FieldVector<ctype, patchDimension>;
    using GlobalCoordinate          = FieldVector<ctype, coorddimension>;
    using JacobianTransposed        = FieldMatrix<ctype, patchDimension, coorddimension>;
    using Hessian                   = FieldMatrix<ctype, patchDimension*(patchDimension + 1) / 2, coorddimension>;
    using Jacobian                  = FieldMatrix<ctype, coorddimension, patchDimension>;
    using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, patchDimension>;
    using JacobianInverse           = FieldMatrix<ctype, patchDimension, coorddimension>;
    using Volume                    = ctype;

    template<int codim>
    using ParameterSpaceGeometry= YaspGeometry<patchDimension-codim,patchDimension,const YaspGrid<patchDimension,TensorProductCoordinates<ctype,patchDimension>>>;

    using ControlPointType = typename IGA::NURBSPatchData<patchDimension, dimworld, ScalarType>::ControlPointType;
    using ControlPointNetType = typename IGA::NURBSPatchData<patchDimension, dimworld, ScalarType>::ControlPointNetType;
    using ControlPointCoordinateNetType = IGA::MultiDimensionNet<dim, typename IGA::NURBSPatchData<patchDimension, dimworld, ScalarType>::GlobalCoordinateType>;
    using Nurbs = Dune::IGA::Nurbs<patchDimension, ScalarType>;
    using NurbsLocalView = typename Nurbs::LocalView;
    template<int codim>
    using GeometryLocalView= PatchGeometryLocalView<codim,NURBSPatchGeometry>;

    template<int codim, typename NURBSPatchGeometry>
    friend struct PatchGeometryLocalView;


   private:
    /* Helper class to compute a matrix pseudo inverse */
    using MatrixHelper = typename MultiLinearGeometryTraits<ctype>::MatrixHelper;

   public:
    NURBSPatchGeometry() = default;

    template<int codim>
    auto localView()const {
      return GeometryLocalView<codim>(*this);
    }


    explicit NURBSPatchGeometry(const IGA::NURBSPatchData<patchDimension, dimworld, ScalarType>& patchData)
        : patchData_(patchData), nurbs_{patchData},localView_(localView<0>()) {}

    /** \brief Map the center of the element to the geometry */
    [[nodiscard]] GlobalCoordinate center() const {
      LocalCoordinate localcenter(0.5);
      return global(localcenter);
    }

    [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
      return computeParameterSpaceCoordinate(*this,global);
    }


    /** \brief Computes the volume of the element with an integration rule for order max(order)*elementdim */
    [[nodiscard]] double volume(int scaleOrder = 1) const {
      const auto rule = Dune::QuadratureRules<ctype, patchDimension>::rule(
          this->type(), scaleOrder * patchDimension * (*std::ranges::max_element(patchData_.degree)));
      ctype vol = 0.0;
      for (auto& gp : rule) {
        localView_.bind(gp.position());
        vol += localView_.integrationElement(gp.position()) * gp.weight();
      }
      return vol;
    }

    [[nodiscard]] bool affine() const { return false; }

    /** \brief Return the number of corners of the element */
    [[nodiscard]] int corners() const { return 1 << patchDimension; }

    /** \brief Return world coordinates of the k-th corner of the element */
    [[nodiscard]] GlobalCoordinate corner(int k) const {
      LocalCoordinate localcorner;
      for (size_t i = 0; i < patchDimension; i++) {
        localcorner[i] = (k & (1 << i)) ? 1 : 0;
      }
      return global(localcorner);
    }

    FieldVector<ctype, dimworld> operator()(const LocalCoordinate& local) const { return global(local); }


    /** \brief evaluates the geometric position
     *
     *  \param[in] u coordinates for each dimension in the [knotSpan.front(), knotSpan.back() ] domain
     */
    [[nodiscard]] FieldVector<ctype, dimworld> global(const LocalCoordinate& u) const {
      auto [nurbsLocalView,cpNet,subNetStart] = calculateNurbsAndControlPointNet(u);
      return this->global(u, nurbsLocalView,cpNet);
    }

    auto zeroFirstAndSecondDerivativeOfPosition(const LocalCoordinate& u) const {
      FieldVector<ctype, dimworld> pos;
      JacobianTransposed J;
      FieldMatrix<ctype, dim*(dim + 1) / 2, coorddimension> H;
      std::array<unsigned int, dim> subDirs;
      for (int subI = 0, i = 0; i < patchDimension; ++i)
        subDirs[subI++] = i;

      const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(u, 2);
      auto cpCoordinateNet                = netOfSpan(u, patchData_.knotSpans, patchData_.degree,
                                                      IGA::extractControlCoordinates(patchData_.controlPoints));
      std::array<unsigned int, patchDimension> ithVecZero{};
      pos = Dune::IGA::dot(basisFunctionDerivatives.get(ithVecZero), cpCoordinateNet);

      for (int dir = 0; dir < patchDimension; ++dir) {
        std::array<unsigned int, patchDimension> ithVec{};
        ithVec[subDirs[dir]] = 1;
        J[dir]               = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet);
      }
      for (int dir = 0; dir < dim; ++dir) {
        std::array<unsigned int, patchDimension> ithVec{};
        ithVec[subDirs[dir]] = 2;  // second derivative in dir direction
        H[dir]               = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet);
      }
      if constexpr (dim > 1) {
        std::array<int, dim> mixeDerivs;
        std::ranges::fill_n(mixeDerivs.begin(), 2, 1);  // first mixed derivatives
        int mixedDireCounter = dim;
        do {
          H[mixedDireCounter++] = dot(basisFunctionDerivatives.get(mixeDerivs), cpCoordinateNet);
        } while (std::ranges::next_permutation(mixeDerivs, std::greater()).found);
      }
      return std::make_tuple(pos, J, H);
    }


    /** \brief compute the Jacobian transposed matrix
     *
     *  \param local coordinates for each dimension
     */
    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& u) const {
      auto [nurbsLocalView,cpNet,subNetStart] = calculateNurbsAndControlPointNet(u);
      return jacobianTransposedImpl(u, nurbsLocalView,cpNet);
    }

    [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const {
      return transpose(jacobianTransposed(local));
    }


    [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
      auto j = jacobianTransposed(local);
      return MatrixHelper::template sqrtDetAAT<patchDimension, coorddimension>(j);
    }

    /** \brief Type of the element: a hypercube of the correct dimension */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(patchDimension); }

    // The following are functions are not part of the Geometry Interface
    [[nodiscard]] std::array<IGA::Utilities::Domain<double>, patchDimension> domain() const {
      std::array<IGA::Utilities::Domain<double>, patchDimension> result{};
      for (int i = 0; i < patchDimension; ++i)
        result[i] = {patchData_.knotSpans[i].front(), patchData_.knotSpans[i].back()};

      return result;
    }

    [[nodiscard]] std::array<int, patchDimension> degree() const { return patchData_.degree; }

    [[nodiscard]] FieldVector<ctype, patchDimension> domainMidPoint() const {
      auto dom = domain();
      Dune::FieldVector<ctype, patchDimension> result{};
      for (int i = 0; i < patchDimension; ++i)
        result[i] = dom[i].center();

      return result;
    }

  private:
    auto calculateNurbsAndControlPointNet(const LocalCoordinate& u) const {
      auto subNetStart = IGA::findSpanCorrected(patchData_.degree, u, patchData_.knotSpans);

      auto cpCoordinateNet = netOfSpan(subNetStart,patchData_.degree,
                                       IGA::extractControlCoordinates(patchData_.controlPoints));
      auto nurbsLocalView = nurbs_.localView();
      nurbsLocalView.bind(subNetStart);
      return std::make_tuple(nurbsLocalView,cpCoordinateNet,subNetStart);
    }

    [[nodiscard]] static FieldVector<ctype, dimworld> global(const LocalCoordinate& u,const NurbsLocalView& nurbsLocalView,
      const ControlPointCoordinateNetType& localControlPointNet)  {
      auto basis           = nurbsLocalView.basisFunctions(u);

      return IGA::dot(basis, localControlPointNet);
    }

    [[nodiscard]] static JacobianTransposed jacobianTransposedImpl(const LocalCoordinate& u,const NurbsLocalView& nurbsLocalView,
      const ControlPointCoordinateNetType& localControlPointNet)  {

      JacobianTransposed result;

      const auto basisFunctionDerivatives = nurbsLocalView.basisFunctionDerivatives(u, 1);

      for (int dir = 0; dir < patchDimension; ++dir) {
        std::array<unsigned int, patchDimension> ithVec{};
        ithVec[dir] = 1;
        result[dir] = dot(basisFunctionDerivatives.get(ithVec), localControlPointNet);
      }
      return result;
    }

   public:
    IGA::NURBSPatchData<patchDimension, dimworld, ScalarType> patchData_;

   private:
    IGA::Nurbs<patchDimension, ScalarType> nurbs_;
    GeometryLocalView<0> localView_;
  };

}  // namespace Dune::IGANEW

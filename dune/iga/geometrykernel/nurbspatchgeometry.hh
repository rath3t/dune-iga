// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file nurbspatchgeometry.hh
 * @brief Definition of the NURBSPatch class.
 * @author Alexander Müller <mueller@ibb.uni-stuttgart.de>
 * @date 2022
 */

#pragma once

#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid/yaspgridgeometry.hh>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"
#include <dune/iga/geometrykernel/geohelper.hh>
#include <dune/iga/geometrykernel/higherorderalgorithms.hh>
#include <dune/iga/geometrykernel/nurbspatchgeometrylocalview.hh>

namespace Dune {
  template <int dim, class Coordinates>
  class YaspGrid;

  template <class ct, int dim>
  class TensorProductCoordinates;
}  // namespace Dune

namespace Dune::IGANEW::GeometryKernel {

  /**
   * @brief NURBSPatch class representing a Non-Uniform Rational B-Spline patch.
   * @tparam dim_ Dimension of the patch in parameter space.
   * @tparam dimworld_ Dimension of the patch in the world space.
   * @tparam ScalarType Scalar type for the patch.
   */
  template <int dim_, int dimworld_, typename ScalarType = double>
  class NURBSPatch {
   public:
    static constexpr std::integral auto worlddimension = dimworld_;
    static constexpr std::integral auto mydimension    = dim_;

    using ctype                        = ScalarType;
    using LocalCoordinate              = FieldVector<ctype, mydimension>;
    using GlobalCoordinate             = FieldVector<ctype, worlddimension>;
    using JacobianTransposed           = FieldMatrix<ctype, mydimension, worlddimension>;
    using Hessian                      = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, worlddimension>;
    using Jacobian                     = FieldMatrix<ctype, worlddimension, mydimension>;
    using JacobianInverseTransposed    = FieldMatrix<ctype, worlddimension, mydimension>;
    using JacobianInverse              = FieldMatrix<ctype, mydimension, worlddimension>;
    using Volume                       = ctype;
    using TensorProductCoordinatesType = std::array<std::vector<ctype>, mydimension>;

    template <int codim>
    using ParameterSpaceGeometry
        = YaspGeometry<mydimension - codim, mydimension,
                       const YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>>;

    using ControlPointType    = typename NURBSPatchData<mydimension, worlddimension, ScalarType>::ControlPointType;
    using ControlPointNetType = typename NURBSPatchData<mydimension, worlddimension, ScalarType>::ControlPointNetType;
    using ControlPointCoordinateNetType
        = MultiDimensionalNet<mydimension,
                              typename NURBSPatchData<mydimension, worlddimension, ScalarType>::GlobalCoordinateType>;
    using Nurbs          = Splines::Nurbs<mydimension, ScalarType>;
    using NurbsLocalView = typename Nurbs::LocalView;
    template <int codim,  typename  TrimmerType_ >
    using GeometryLocalView = PatchGeometryLocalView<codim, NURBSPatch, TrimmerType_>;

    template <int codim, typename NURBSPatch, typename TrimmerType_>
    friend struct PatchGeometryLocalView;

   private:
    /* @brief Helper class to compute a matrix pseudo-inverse. */
    using MatrixHelper = typename MultiLinearGeometryTraits<ctype>::MatrixHelper;

   public:
    /* @brief Default constructor for NURBSPatch.*/
    NURBSPatch() = default;

    /**
     * @brief Get a local view of the NURBS patch.
     * @tparam codim Codimension of the patch.
     * @tparam TrimmerType Type of the trimmer.
     * @return Local view of the patch.
     */
    template <int codim, typename TrimmerType>
    auto localView() const {
      return GeometryLocalView<codim, TrimmerType>(*this);
    }

    /**
     * @brief Explicit constructor for NURBSPatch.
     * @param patchData Patch data for the NURBS patch.
     */
    explicit NURBSPatch(const NURBSPatchData<dim_, dimworld_, ScalarType>& patchData)
        : patchData_(patchData),
          uniqueKnotSpans_{Splines::createUniqueKnotSpans(patchData.knotSpans)},
          nurbs_{patchData_} {}

    /**
     * @brief Explicit constructor for NURBSPatch with unique knot spans.
     * @param patchData Patch data for the NURBS patch.
     * @param uniqueKnotSpans Unique knot spans for the NURBS patch.
     */
    explicit NURBSPatch(const NURBSPatchData<dim_, dimworld_, ScalarType>& patchData,
                        const std::array<std::vector<ctype>, dim_>& uniqueKnotSpans)
        : patchData_(patchData), uniqueKnotSpans_{uniqueKnotSpans}, nurbs_{patchData_} {}

    /**
     * @brief Get the center of the element mapped to the geometry.
     * @return Global coordinate of the center.
     */
    [[nodiscard]] GlobalCoordinate center() const { return global(domainMidPoint()); }

    /**
     * @brief Get the local coordinate corresponding to a given global coordinate.
     * @param global Global coordinate in world space.
     * @return Local coordinate in parameter space.
     */
    [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
      return findClosestParameterSpaceCoordinate(*this, global);
    }

    /**
     * @brief Compute the volume of the element with an integration rule.
     * @return Volume of the patch.
     */
    [[nodiscard]] double volume() const {
      const auto rule = Dune::QuadratureRules<ctype, mydimension>::rule(
          this->type(), mydimension * (*std::ranges::max_element(patchData_.degree)));
      ctype vol = 0.0;
      for (auto& gp : rule) {
        vol += integrationElement(gp.position()) * gp.weight();
      }
      return vol;
    }

    /**
     * @brief Get the number of corners of the patch.
     * @return Number of corners.
     */
    [[nodiscard]] int corners() const { return 1 << mydimension; }

    /**
     * @brief Get the world coordinates of the k-th corner of the patch.
     * @param k Index of the corner.
     * @return Global coordinate of the corner.
     */
    [[nodiscard]] GlobalCoordinate corner(int k) const {
      LocalCoordinate localcorner;
      for (size_t i = 0; i < mydimension; i++) {
        localcorner[i] = (k & (1 << i)) ? uniqueKnotVector()[i].back() : uniqueKnotVector()[i].front();
      }
      return global(localcorner);
    }

    /**
     * @brief Evaluate the geometric position for a given set of coordinates in the parameter space.
     * @param u Local coordinates for each dimension in the [knotSpan.front(), knotSpan.back()] domain.
     * @return Global coordinate of the evaluated position.
     */
    [[nodiscard]] FieldVector<ctype, worlddimension> global(const LocalCoordinate& u) const {
      auto [nurbsLocalView, cpNet, subNetStart] = calculateNurbsAndControlPointNet(u);
      return GeometryKernel::position(u, nurbsLocalView, cpNet);
    }

    /**
     * @brief Calculate zero first and second derivatives of the position for a given local coordinate.
     * @param u Local coordinates for each dimension in the [knotSpan.front(), knotSpan.back()] domain.
     * @return Tuple containing position, transposed Jacobian, and Hessian matrices.
     */
    auto zeroFirstAndSecondDerivativeOfPosition(const LocalCoordinate& u) const {
      auto [nurbsLocalView, cpNet, subNetStart] = calculateNurbsAndControlPointNet(u);
      return zeroFirstAndSecondDerivativeOfPositionImpl(u, nurbsLocalView, cpNet);
    }

    /**
     * @brief Compute the Jacobian transposed matrix for a given local coordinate.
     * @param u Local coordinates for each dimension in the [knotSpan.front(), knotSpan.back()] domain.
     * @return Jacobian transposed matrix.
     */
    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& u) const {
      auto [nurbsLocalView, cpNet, subNetStart] = calculateNurbsAndControlPointNet(u);
      return GeometryKernel::jacobianTransposed(u, nurbsLocalView, cpNet);
    }

    /**
     * @brief Compute the Hessian matrix for a given local coordinate.
     * @param u Local coordinates for each dimension in the [knotSpan.front(), knotSpan.back()] domain.
     * @return Hessian matrix.
     */
    [[nodiscard]] Hessian hessian(const LocalCoordinate& u) const {
      auto [nurbsLocalView, cpNet, subNetStart] = calculateNurbsAndControlPointNet(u);
      return GeometryKernel::hessian(u, nurbsLocalView, cpNet);
    }

    /**
     * @brief Compute the Jacobian matrix for a given local coordinate.
     * @param u Local coordinates for each dimension in the [knotSpan.front(), knotSpan.back()] domain.
     * @return Jacobian matrix.
     */
    [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const { return transpose(jacobianTransposed(local)); }

    /**
     * @brief Compute the integration element for a given local coordinate.
     * @param u Local coordinates for each dimension in the [knotSpan.front(), knotSpan.back()] domain.
     * @return Integration element.
     */
    [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
      auto j = jacobianTransposed(local);
      return MatrixHelper::template sqrtDetAAT<mydimension, worlddimension>(j);
    }

    /**
     * @brief Get the type of the element: a hypercube of the correct dimension.
     */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }

    /**
     * @brief Get the domain of the element in parameter space.
     * @return Array of Utilities::Domain for each dimension.
     */
    [[nodiscard]] std::array<Utilities::Domain<double>, mydimension> domain() const {
      std::array<Utilities::Domain<double>, mydimension> result{};
      for (int i = 0; i < mydimension; ++i)
        result[i] = {patchData_.knotSpans[i].front(), patchData_.knotSpans[i].back()};
      return result;
    }

    /**
     * @brief Get the polynomial degree of the patch in each dimension.
     * @return Array of degrees for each dimension.
     */
    [[nodiscard]] std::array<int, mydimension> degree() const { return patchData_.degree; }

    /**
     * @brief Get the domain midpoint of the patch.
     * @return FieldVector representing the midpoint of the element domain.
     */
    [[nodiscard]] FieldVector<ctype, mydimension> domainMidPoint() const {
      auto dom = domain();
      Dune::FieldVector<ctype, mydimension> result{};
      for (int i = 0; i < mydimension; ++i)
        result[i] = dom[i].center();
      return result;
    }

    /**
     * @brief Get the number of control points in each direction.
     * @return Array containing the number of control points in each direction.
     */
    auto numberOfControlPoints() const { return patchData_.controlPoints.strideSizes(); }

    /**
     * @brief Get the number of non-zero measure spans in each direction.
     * @return Array containing the number of spans in each direction.
     */
    auto numberOfSpans() const {
      std::array<int, mydimension> spansPerDirection;
      for (int i = 0; i < mydimension; ++i)
        spansPerDirection[i] = uniqueKnotSpans_[i].size() - 1;
      return spansPerDirection;
    }

    /**
     * @brief Get the patch data of the NURBS patch.
     * @return Reference to the NURBSPatchData.
     */
    const auto& patchData() const { return patchData_; }

    /* @brief Get the unique knot vector spans of the NURBS patch. */
    const TensorProductCoordinatesType& uniqueKnotVector() const { return uniqueKnotSpans_; }

    /* @brief Get the patch data of the NURBS patch. */
    auto& patchData() { return patchData_; }

   private:
    /* @brief Calculate NURBS and control point net for a given local coordinate. */
    auto calculateNurbsAndControlPointNet(const LocalCoordinate& u) const {
      auto subNetStart     = Splines::findSpan(patchData_.degree, u, patchData_.knotSpans);
      auto cpCoordinateNet = Splines::netOfSpan(subNetStart, patchData_.degree,
                                                Splines::extractControlCoordinates(patchData_.controlPoints));
      auto nurbsLocalView  = nurbs_.localView();
      nurbsLocalView.bind(subNetStart);
      return std::make_tuple(nurbsLocalView, cpCoordinateNet, subNetStart);
    }

    /* @brief Calculate zero first and second derivatives of the position for a given local coordinate. */
    static auto zeroFirstAndSecondDerivativeOfPositionImpl(const LocalCoordinate& u,
                                                           const NurbsLocalView& nurbsLocalView,
                                                           const ControlPointCoordinateNetType& localControlPointNet) {
      return std::make_tuple(GeometryKernel::position(u, nurbsLocalView, localControlPointNet),
                             GeometryKernel::jacobianTransposed(u, nurbsLocalView, localControlPointNet),
                             GeometryKernel::hessian(u, nurbsLocalView, localControlPointNet));
    }

    NURBSPatchData<mydimension, worlddimension, ScalarType> patchData_;
    TensorProductCoordinatesType uniqueKnotSpans_;
    Nurbs nurbs_;
  };

}  // namespace Dune::IGANEW::GeometryKernel

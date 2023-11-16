// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/geometry/closestpointprojection.hh"
#include "dune/iga/geometry/geohelper.hh"
#include "dune/iga/nurbsalgorithms.hh"
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/hierarchicpatch/geometryalgorithms.hh>

namespace Dune::IGANEW {

  template <int codim, typename PatchGeometry >
    struct PatchGeometryLocalView {

      using ctype                     = typename PatchGeometry::ctype;
      static constexpr int gridDimension = PatchGeometry::patchDimension;
      static constexpr int mydimension = gridDimension-codim;
      static constexpr bool trim = PatchGeometry::trim;

      static constexpr std::integral auto coorddimension = PatchGeometry::coorddimension;

      using LocalCoordinate           = FieldVector<ctype, mydimension>;
      using GlobalCoordinate          = typename PatchGeometry::GlobalCoordinate;
      using JacobianTransposed        = FieldMatrix<ctype, mydimension, coorddimension>;
      using JacobianTransposedInPatch        = FieldMatrix<ctype, gridDimension, coorddimension>;
      // TODO trim ParameterSpaceGeometry
      using ParameterSpaceGeometry        = typename PatchGeometry::template ParameterSpaceGeometry<codim>;

    //! if we have codim==0, then the Jacobian in the parameter space of the grid entity itself is a DiagonalMatrix, and
    // Coordinates in a single knot span differ from coordinates on the B-spline patch
    // by an affine transformation.  This transformation is stored in the diagonal entries.
    // If trimming is disabled the Jacobian in the parameter space of subentities (edges or surfaces) is a DiagonalMatrixBlock but is treated as a FieldMatrix, since
    // is also just cubes and an axis-aligned geometry,
    // if trimming is enabled the Jacobian of the parameterspace for subentities is a potentially fully populated FieldMatrix,
    // Think about an arbitrary curve in the axis-aligned cube parameter grid
    //                      ------->
    // `____________________C_______D_______`
    // ::````````````````::``````````'|````::
    // ::                ::          `|    ::
    // ::              ^F::          `|    ::
    // ::              | ::          `|    ::
    // ::              | ::          `|    ::
    // ::              |E::          `|    ::
    // ::                ::          .|  x ::
    // ::                ::          _:    ::
    // ::................::..........|:....::
    // ::''''''''''''''''::''''''''''/'''''::
    // ::                ::         _:B    ::
    // ::                ::        .|      ::
    // ::                ::      `::A      ::
    // ::           '::::||::::::'         ::
    // ::         ::`    ::          x     ::
    // ::        |'      ::                ::
    // ::        |   x   ::                ::
    // ::````````|```````::````````````````::
    // `____________________________________`
    // at point A the tangent of the curve (the JacobianInParameterSpace) points into the oblique direction (to B)
    // for the axis-aligned geometry these tangents are from E to F or from C to D, depending on the orientation of the subentity, which are aligned with the axis
      using JacobianTransposedInParameterSpace        = typename ParameterSpaceGeometry::JacobianTransposed;
      using GlobalInParameterSpace        = typename ParameterSpaceGeometry::GlobalCoordinate;

      using Hessian                   = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, coorddimension>;
      using Jacobian                  = FieldMatrix<ctype, coorddimension, mydimension>;
      using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;
      using JacobianInverse           = FieldMatrix<ctype, mydimension, coorddimension>;
      using MatrixHelper = typename PatchGeometry::MatrixHelper;
      using Volume                    = ctype;

    using Nurbs = IGA::Nurbs<gridDimension, ctype>;
    using NurbsLocalView = typename Nurbs::LocalView;

    using ControlPointCoordinateNetType= typename PatchGeometry::ControlPointCoordinateNetType;

    PatchGeometryLocalView() =default;
      explicit PatchGeometryLocalView(const PatchGeometry& patchGeometry):patchGeometry_{&patchGeometry}{}

      [[nodiscard]] GlobalCoordinate center() const {
        LocalCoordinate localcenter(0.5);
        return global(localcenter);
      }

      void bind(const ParameterSpaceGeometry& lGeo) {
        parameterSpaceGeometry=std::make_optional<ParameterSpaceGeometry>(lGeo);
        std::tie(nurbsLocalView_,localControlPointNet,spanIndices_)= patchGeometry_->calculateNurbsAndControlPointNet(lGeo.center());
      }

      /** \brief evaluates the geometric position
 *
 *  \param[in] u coordinates for each dimension in the [knotSpan.front(), knotSpan.back() ] domain
 */
      [[nodiscard]] GlobalCoordinate global(const LocalCoordinate& u) const {
        assert(parameterSpaceGeometry && "Bind the local view first!");
        return patchGeometry_->global(globalInParameterSpace(u), nurbsLocalView_,localControlPointNet);
      }

      [[nodiscard]] JacobianTransposedInParameterSpace jacobianTransposedInParameterSpace(const LocalCoordinate& u) const {
        assert(parameterSpaceGeometry && "Bind the local view first!");
        return parameterSpaceGeometry.value().jacobianTransposed(u);
      }

    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& u) const {
        assert(parameterSpaceGeometry && "Bind the local view first!");
        // using JacobianTransposed        = FieldMatrix<ctype, mydimension, coorddimension>;
        const auto ouInPatch = globalInParameterSpace(u);
        const JacobianTransposedInPatch JTinPatch = patchGeometry_->jacobianTransposedImpl(ouInPatch, nurbsLocalView_,localControlPointNet);
        const JacobianTransposedInParameterSpace jInparaMeterSpace = jacobianTransposedInParameterSpace(u);
        JacobianTransposed jacobianOfEntity= jInparaMeterSpace*JTinPatch;
        return jacobianOfEntity;
      }

      [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const {
        assert(parameterSpaceGeometry && "Bind the local view first!");

        return transpose(jacobianTransposed(local));
      }

      [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
        assert(parameterSpaceGeometry && "Bind the local view first!");
        auto j = jacobianTransposed(local);
        return MatrixHelper::template sqrtDetAAT<mydimension, coorddimension>(j);
      }

      [[nodiscard]] double volume(int scaleOrder = 1) const {
        assert(parameterSpaceGeometry && "Bind the local view first!");

        if constexpr(trim) {
          // if constexpr (mydimension == 2)
            // if (subgrid_) {
            //   Dune::QuadratureRule<double, mydimension> rule;
            //   fillQuadratureRuleImpl(rule, *subgrid_.get(), (*std::ranges::max_element(patchData().degree)));
            //   Volume vol = 0.0;
            //   for (auto& gp : rule)
            //     vol += integrationElement(gp.position()) * gp.weight();
            //   return vol;
            // }
        }

        const auto& rule = Dune::QuadratureRules<ctype, mydimension>::rule(
            GeometryTypes::cube(mydimension), (*std::ranges::max_element(patchData().degree)));
        ctype vol = 0.0;
        for (auto& gp : rule)
          vol += integrationElement(gp.position()) * gp.weight();
        return vol;
      }

      [[nodiscard]] GeometryType type () const {
        return GeometryTypes::cube(mydimension);
      }

      /** \brief Return the number of corners of the element */
      [[nodiscard]] int corners() const { return 1 << mydimension; }

      /** \brief Return world coordinates of the k-th corner of the element */
      [[nodiscard]] GlobalCoordinate corner(int k) const {
        assert(parameterSpaceGeometry && "Bind the local view first!");

        LocalCoordinate localcorner;
        for (size_t i = 0; i < mydimension; i++) {
          localcorner[i] = (k & (1 << i)) ? 1 : 0;
        }
        return global(localcorner);
      }

      /** \brief Inverse of global this function gets a point defined in the world space and return
 * the closest point in local coordinates, i.e. in [0,1] domain for each grid dimension
 *
 *  \param global global coordinates for the point where the local coordinates are searched for
 */
      [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
        return computeParameterSpaceCoordinate(*this,global);
      }

      [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
        JacobianInverseTransposed jacobianInverseTransposed1;
        const JacobianTransposed Jt = jacobianTransposed(local);
        MatrixHelper::template rightInvA<mydimension, coorddimension>(Jt,
                                                                      jacobianInverseTransposed1);
        return jacobianInverseTransposed1;
      }

    /**
 * @brief Compute the first and second derivatives of the composition function h(g_1(t), g_2(t)) with respect to the parameter t.
 *
 * This function calculates the Jacobian and Hessian matrices for the composition function h(g_1(t), g_2(t)) where g(t) is the original curve
 * and h(u, v) is the parametrization function. The computation involves chain and product rules between g(t) and h(u, v).
 *
 * @param[in] t The parameter value for the curve g(t).
 * @param[in] g1_t The first component of the curve g(t).
 * @param[in] g2_t The second component of the curve g(t).
 *
 * @return std::pair<Eigen::Matrix<double, 2, 1>, Eigen::Matrix<double, 2, 2>> A pair containing the Jacobian matrix and the Hessian matrix.
 *
 * @details
 * The Jacobian matrix J is given by:
 * \f[
 * J = \begin{bmatrix}
 * \frac{\partial h}{\partial g_1} \cdot \frac{\partial g_1}{\partial t} + \frac{\partial h}{\partial g_2} \cdot \frac{\partial g_2}{\partial t} \\
 * \frac{\partial h}{\partial g_1} \cdot \frac{\partial g_1}{\partial t} \\
 * \end{bmatrix}
 * \f]
 *
 * The Hessian matrix H is given by:
 * \f[
 * H = \begin{bmatrix}
 * \frac{\partial^2 h}{\partial g_1^2} \cdot \left(\frac{\partial g_1}{\partial t}\right)^2 + 2 \cdot \frac{\partial^2 h}{\partial g_1 \partial g_2} \cdot \frac{\partial g_1}{\partial t} \cdot \frac{\partial g_2}{\partial t} + \frac{\partial^2 h}{\partial g_2^2} \cdot \left(\frac{\partial g_2}{\partial t}\right)^2 & \frac{\partial^2 h}{\partial g_1^2} \cdot \frac{\partial g_1}{\partial t} + \frac{\partial^2 h}{\partial g_1 \partial g_2} \cdot \frac{\partial g_2}{\partial t} \\
 * \frac{\partial^2 h}{\partial g_1^2} & \frac{\partial^2 h}{\partial g_1 \partial g_2} \\
 * \end{bmatrix}
 * \f]
 */
      std::tuple<GlobalCoordinate,JacobianTransposed,Hessian> zeroFirstAndSecondDerivativeOfPosition(const LocalCoordinate& u) const {
        assert(parameterSpaceGeometry && "Bind the local view first!");
        const auto ouInPatch = globalInParameterSpace(u);
        parameterSpaceGeometry.value().zeroFirstAndSecondDerivativeOfPosition(ouInPatch);


      return {};
      }

      [[nodiscard]] LocalCoordinate domainMidPoint() const {
        auto dom = domain();
        LocalCoordinate result{};
        for (int i = 0; i < mydimension; ++i)
          result[i] = dom[i].center();

        return result;
      }

      [[nodiscard]] std::array<IGA::Utilities::Domain<double>, mydimension> domain() const { return {}; }

      [[nodiscard]] bool affine() const { return false; }
      [[nodiscard]] const std::array<int, gridDimension>&  spanIndices() const { return spanIndices_; }
      [[nodiscard]] const PatchGeometry&  patchGeometry() const { return *patchGeometry_; }
      [[nodiscard]] const auto&  patchData() const { return patchGeometry_->patchData_; }
      [[nodiscard]] const NurbsLocalView&  nurbs() const { return nurbsLocalView_; }
      [[nodiscard]] const ControlPointCoordinateNetType&  controlPointCoordinates() const { return localControlPointNet; }
    private:

      GlobalInParameterSpace globalInParameterSpace(const LocalCoordinate& local)const {
        assert(parameterSpaceGeometry && "Bind the local view first!");
        return  parameterSpaceGeometry.value().global(local);
      }
      ControlPointCoordinateNetType localControlPointNet;
      NurbsLocalView nurbsLocalView_;
      std::array<int, gridDimension> spanIndices_;
      const PatchGeometry* patchGeometry_;
      std::optional<ParameterSpaceGeometry> parameterSpaceGeometry;
    };


}  // namespace Dune::IGANEW

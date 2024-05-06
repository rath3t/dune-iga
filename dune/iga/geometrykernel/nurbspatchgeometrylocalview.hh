// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file nurbspatchgeometrylocalview.hh
 * @brief Definition of the NURBSPatch geometry local view class.
 * @author Alexander Müller <mueller@ibb.uni-stuttgart.de>
 * @date 2022
 */

#pragma once

#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/geometrykernel/algorithms.hh>
#include <dune/iga/hierarchicpatch/enums.hh>
// #include <dune/iga/hierarchicpatch/hierachicpatchgridlocalgeometry.hh>
// #include <dune/iga/hierarchicpatch/patchgridentity.hh>
#include <dune/iga/splines/nurbsalgorithms.hh>
// #include <dune/iga/trimmer/identitytrimmer/trimmer.hh>

namespace Dune::IGANEW {
namespace GeometryKernel {
  namespace Impl {
    struct IndexPair
    {
      int row;
      int col;
    };

    template <int dim>
    auto& voigtIndices() {
      if constexpr (dim == 1) {
        static std::array<IndexPair, dim*(dim + 1) / 2> voigt1D = {
            {0, 0}
        };
        return voigt1D;
      } else if constexpr (dim == 2) {
        static std::array<IndexPair, dim*(dim + 1) / 2> voigt2D = {
            {{0, 0}, {1, 1}, {0, 1}}
        };
        return voigt2D;
      } else if constexpr (dim == 3) {
        static std::array<IndexPair, dim*(dim + 1) / 2> voigt3D = {
            {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}}
        };
        return voigt3D;
      }
    }
  } // namespace Impl

  /**
   * @brief Represents a local view of a patch geometry.
   *
   * @tparam codim Codimension of the patch geometry.
   * @tparam PatchGeometry Type of the patch geometry.
   * @tparam TrimmerType_ Type of the trimmer.
   */
  template <int codim, typename PatchGeometry, typename TrimmerType_, typename LocalParameterSpaceGeometry = void>
  struct PatchGeometryLocalView
  {
    using ctype                        = typename PatchGeometry::ctype; ///< Scalar type for coordinates.
    static constexpr int gridDimension = PatchGeometry::mydimension;    ///< Dimension of the underlying geometry grid.
    static constexpr int mydimension   = gridDimension - codim;         ///< Effective dimension of the local geometry.
    static constexpr int numberOfSecondDerivatives =
        mydimension * (mydimension + 1) / 2; ///< Number of second derivatives of the local view.
    static constexpr int patchNumberOfSecondDerivatives =
        gridDimension * (gridDimension + 1) / 2;                     ///< Number of second derivatives for the patch.
    using Trimmer            = TrimmerType_;                         ///< Type of the associated trimmer.
    using ParameterSpaceGrid = typename Trimmer::ParameterSpaceGrid; ///< Type of the parameter space grid.

    static constexpr std::integral auto worlddimension = PatchGeometry::worlddimension; ///< Dimension of the world.
    static constexpr std::integral auto coorddimension = worlddimension;                ///< Dimension of the world.

    using LocalCoordinate  = FieldVector<ctype, mydimension>;          ///< Type for local coordinates.
    using GlobalCoordinate = typename PatchGeometry::GlobalCoordinate; ///< Type for global coordinates.
    using JacobianTransposed =
        FieldMatrix<ctype, mydimension, coorddimension>; ///< Type for the transposed Jacobian matrix.
    using PatchJacobianTransposed =
        typename PatchGeometry::JacobianTransposed;       ///< Type for the transposed Jacobian matrix of the patch.
    using PatchHessian = typename PatchGeometry::Hessian; ///< Type for the Hessian matrix of the patch.

  private:
    static constexpr bool isParameterSpaceGeometryProvided = not std::is_same_v<LocalParameterSpaceGeometry, void>;

  public:
    using ParameterSpaceGeometry =
        std::conditional_t<isParameterSpaceGeometryProvided, LocalParameterSpaceGeometry,
                           typename Trimmer::template Codim<codim>::LocalParameterSpaceGeometry>;

    // if we have codim==0, then the Jacobian in the parameter space of the grid entity itself is a DiagonalMatrix,
    // and
    // Coordinates in a single knot span differ from coordinates on the B-spline patch
    // by an affine transformation.  This transformation is stored in the diagonal entries.
    // If trimming is disabled the Jacobian in the parameter space of subentities (edges or surfaces) is a
    // DiagonalMatrixBlock but is treated as a FieldMatrix, since is also just cubes and an axis-aligned geometry, if
    // trimming is enabled the Jacobian of the parameterspace for subentities is a potentially fully populated
    // FieldMatrix, Think about an arbitrary curve in the axis-aligned cube parameter grid
    //                     ------->
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
    // for the axis-aligned geometry these tangents are from E to F or from C to D, depending on the orientation of
    // the subentity, which are aligned with the axis
    using JacobianTransposedInParameterSpace = typename ParameterSpaceGeometry::JacobianTransposed;
    using GlobalInParameterSpace             = typename ParameterSpaceGeometry::GlobalCoordinate;

    using Hessian                   = FieldMatrix<ctype, numberOfSecondDerivatives, coorddimension>;
    using Jacobian                  = FieldMatrix<ctype, coorddimension, mydimension>;
    using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;
    using JacobianInverse           = FieldMatrix<ctype, mydimension, coorddimension>;
    using MatrixHelper              = typename PatchGeometry::MatrixHelper;
    using Volume                    = ctype;

    using Nurbs          = Splines::Nurbs<gridDimension, ctype>;
    using NurbsLocalView = typename Nurbs::LocalView;

    using ControlPointCoordinateNetType = typename PatchGeometry::ControlPointCoordinateNetType;

    PatchGeometryLocalView() = default;
    explicit PatchGeometryLocalView(const PatchGeometry& patchGeometry)
        : patchGeometry_{&patchGeometry} {}

    /**
     * @brief Get the center of the patch geometry.
     * @return Global coordinates of the center.
     */
    [[nodiscard]] GlobalCoordinate center() const {
      return global(LocalCoordinate(0.5));
    }

    /**
     * @brief Bind the local view to a parameter space geometry.
     * @param lGeo Parameter space geometry.
     */
    void bind(const ParameterSpaceGeometry& lGeo) {
      // static_assert(std::is_same_v<GlobalInParameterSpace,FieldVector<ctype,gridDimension>>);
      parameterSpaceGeometry = std::make_shared<ParameterSpaceGeometry>(lGeo);
      std::tie(nurbsLocalView_, localControlPointNet, spanIndices_) =
          patchGeometry_->calculateNurbsAndControlPointNet(lGeo.center());
    }

    /**
     * @brief Get the global position at a local coordinate.
     * @param local Local coordinate, i.e. a tuple where each coordinate is in [0,1] domain for each local view
     * dimension
     * @return Global coordinates of the position.
     */
    [[nodiscard]] GlobalCoordinate global(const LocalCoordinate& local) const {
      checkState();
      return GeometryKernel::position(globalInParameterSpace(local), nurbsLocalView_, localControlPointNet);
    }

    /**
     * @brief Get the transposed Jacobian matrix at a local coordinate.
     * @param local Local coordinate, i.e. a tuple where each coordinate is in [0,1] domain for each local view
     * dimension
     * @return Transposed Jacobian matrix.
     */
    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
      checkState();
      const PatchJacobianTransposed dfdg =
          GeometryKernel::jacobianTransposed(globalInParameterSpace(local), nurbsLocalView_, localControlPointNet);
      const JacobianTransposedInParameterSpace dgdt = jacobianTransposedInParameterSpace(local);
      return dgdt * dfdg;
    }

    // todo check
    auto secondFundamentalForm(const LocalCoordinate& local) const {
      return patchGeometry_->secondFundamentalForm(local);
    }

    // // todo check
    // auto zeroFirstAndSecondDerivativeOfPosition(const LocalCoordinate& local) const {
    //   return patchGeometry_->zeroFirstAndSecondDerivativeOfPosition(local);
    // }

    /**
     * @brief Get the Jacobian matrix at a local coordinate.
     * @param local Local coordinate, i.e. a tuple where each coordinate is in [0,1] domain for each local view
     * dimension
     * @return Jacobian matrix.
     */
    [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const {
      return transpose(jacobianTransposed(local));
    }

    /**
     * @brief Compute the integration element at a local coordinate.
     * @param local Local coordinate, i.e. a tuple where each coordinate is in [0,1] domain for each local view
     * dimension
     * dimension
     * @return Integration element.
     */
    [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
      auto jT = jacobianTransposed(local);
      return MatrixHelper::template sqrtDetAAT<mydimension, coorddimension>(jT);
    }

    /**
     * @brief Compute the volume of the patch.
     * @return Volume of the patch.
     */
    [[nodiscard]] double volume() const {
      if constexpr (not Trimmer::isAlwaysTrivial) {
        // @todo: Implement integration of trimmed quantities and new edge geometries.
      } else {
        const auto& rule = QuadratureRules<ctype, mydimension>::rule(
            GeometryTypes::cube(mydimension), (*std::ranges::max_element(patchGeometry_->patchData_.degree)));
        Volume vol = 0.0;
        for (auto& gp : rule)
          vol += integrationElement(globalInParameterSpace(gp.position())) * gp.weight() *
                 parameterSpaceGeometry->volume();
        return vol;
      }
    }

    /** @brief Type of the type of the parameter space element */
    [[nodiscard]] GeometryType type() const {
      return parameterSpaceGeometry->type();
    }

    /** @brief Return the number of corners of the element */
    [[nodiscard]] int corners() const {
      return parameterSpaceGeometry->corners();
    }

    /** @brief Return world coordinates of the k-th corner of the element */
    [[nodiscard]] GlobalCoordinate corner(int k) const {
      checkState();
      return GeometryKernel::position(parameterSpaceGeometry->corner(k), nurbsLocalView_, localControlPointNet);
    }

    /**
     * @brief Compute the local coordinate corresponding to a global position.
     * @param global Global coordinates for which the local coordinates are sought.
     * @return Local coordinates in the [0,1] domain for each grid dimension.
     */
    [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
      return GeometryKernel::findClosestParameterSpaceCoordinate(*this, global);
    }

    /**
     * @brief Get the transposed inverse Jacobian matrix at a local coordinate.
     * @param local Local coordinate, i.e. a tuple where each coordinate is in [0,1] domain for each local view
     * dimension
     * @return Transposed inverse Jacobian matrix.
     */
    [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
      JacobianTransposed Jt = jacobianTransposed(local);
      JacobianInverseTransposed jacobianInverseTransposed;
      MatrixHelper::template rightInvA<mydimension, coorddimension>(Jt, jacobianInverseTransposed);
      return jacobianInverseTransposed;
    }

    /**
     * @brief Get the inverse Jacobian matrix at a local coordinate.
     * @param local Local coordinate, i.e. a tuple where each coordinate is in [0,1] domain for each local view
     * dimension
     * @return inverse Jacobian matrix.
     */
    [[nodiscard]] JacobianInverse jacobianInverse(const LocalCoordinate& local) const {
      return transpose(jacobianInverseTransposed(local));
    }

    /**
     * @brief Compute the zeroth, first, and second derivatives of the position at a local coordinate.
     * @param local Local coordinate, i.e. a tuple where each coordinate is in [0,1] domain for each local view
     * dimension
     * @return Tuple of the global position, transposed Jacobian matrix, and Hessian matrix.
     */
    std::tuple<GlobalCoordinate, JacobianTransposed, Hessian> zeroFirstAndSecondDerivativeOfPosition(
        const LocalCoordinate& local) const {
      checkState();

      return std::make_tuple(global(local), jacobianTransposed(local), hessian(local));
    }

    /**
     * @brief Computes the second derivatives of the local view.
     *
     * @param local Local coordinate, i.e. a tuple where each coordinate is in [0,1] domain for each local view
     * dimension
     * @return A tuple of the global position, the Jacobian matrix and Hessian matrix
     * @details
     * **Example curve on surface**
     *
     * The local view contains a `parameterSpaceGeometry`, which describes the  composition function
     *
     * \f[g: \begin{cases}\mathbb{R} \rightarrow \mathbb{R}^2 \\ t \mapsto g(t) \end{cases},\f]
     *
     * which describes the curve's local parametrization in terms of the patch geometries parametrization
     * The patch geometry is described by
     *
     * \f[f: \begin{cases}\mathbb{R}^2 \rightarrow \mathbb{R}^3 \\ (u,v) \mapsto * f(u,v) \end{cases},\f]
     *
     * which describes a surface in 3D.
     * Then the curve's realization on the surface in the 3D space can be written by
     *
     * \f[h: \begin{cases}\mathbb{R} \rightarrow \mathbb{R}^3 \\ t \mapsto f(g_1(t),g_2(t)) \end{cases},\f]
     *
     * This function calculates the Jacobian and Hessian matrices for the composition function \f$h(g_1(t), g_2(t))\f$
     * where \f$g(t)\f$ is the original curve and \f$h(u, v)\f$  is the patch geometry function. The computation
     * involves chain and product rules between \f$g(t)\f$ and \f$h(u, v)\f$.
     *     *
     * The Jacobian matrix J is given by
     * \f[
     * J =
     * \frac{\partial f }{\partial g_1}  \frac{\partial g_1}{\partial t} + \frac{\partial f}{\partial g_2}
     * \frac{\partial g_2}{\partial t} \f] and has dimensions \f$3\times 1\f$ or in general
     * \f$\verb+worlddimension+\times \verb+mydimension+\f$.
     *
     * The Hessian matrix H is given by:
     * \f[
     * H =
     * \frac{\partial^2 f}{\partial g_1\partial g_1}  \left(\frac{\partial g_1}{\partial t}\right)^2 +\frac{\partial^2
     * f}{\partial g_2\partial g_2}  \left(\frac{\partial g_2}{\partial t}\right)^2 + 2\frac{\partial^2 f}{\partial
     * g_1\partial g_2}  \frac{\partial g_1}{\partial t} \frac{\partial g_2}{\partial t} + \frac{\partial f}{\partial
     * g_1}\frac{\partial^2 g_1}{\partial t^2}+ \frac{\partial f}{\partial g_2}\frac{\partial^2 g_2}{\partial t^2} \f]
     *
     * and has dimensions \f$1\times 3\f$ or in general \f$\verb|mydimension*(mydimension+1)/2| \times
     * \verb+worlddimension+\f$.
     */
    Hessian hessian(const LocalCoordinate& local) const {
      checkState();
      const auto ouInPatch = globalInParameterSpace(local);

      const auto [p, dfdg, dfdgdg] =
          patchGeometry_->zeroFirstAndSecondDerivativeOfPositionImpl(ouInPatch, nurbsLocalView_, localControlPointNet);

      const JacobianTransposedInParameterSpace dgdt = jacobianTransposedInParameterSpace(local);
      static_assert(JacobianTransposedInParameterSpace::rows == mydimension);
      static_assert(JacobianTransposedInParameterSpace::cols == gridDimension);
      // dgdt (1x2, curve on surface), (2x2 surface in surface), (2x3 surface in 3D), (3x3 volume in 3D)
      FieldMatrix<ctype, patchNumberOfSecondDerivatives, numberOfSecondDerivatives> dgdtSquared(0);
      for (int patchIndex = 0; auto [patchRow, patchCol] : Impl::voigtIndices<gridDimension>()) {
        for (int myIndex = 0; auto [myRow, myCol] : Impl::voigtIndices<mydimension>()) {
          if constexpr (not std::is_same_v<JacobianTransposedInParameterSpace, DiagonalMatrix<ctype, gridDimension>>)
            dgdtSquared[patchIndex][myIndex] =
                (dgdt[myRow][patchRow] * dgdt[myCol][patchCol] + dgdt[myCol][patchRow] * dgdt[myRow][patchCol]) *
                ((patchRow == patchCol) ? 0.5 : 1.0);
          else {
            const auto dgdtmyRowpatchRow = (myRow == patchRow) ? dgdt[myRow][patchRow] : 0;
            const auto dgdtmyColpatchCol = (myCol == patchCol) ? dgdt[myCol][patchCol] : 0;
            const auto dgdtmyColpatchRow = (myCol == patchRow) ? dgdt[myCol][patchRow] : 0;
            const auto dgdtmyRowpatchCol = (myRow == patchCol) ? dgdt[myRow][patchCol] : 0;
            dgdtSquared[patchIndex][myIndex] =
                (dgdtmyRowpatchRow * dgdtmyColpatchCol + dgdtmyColpatchRow * dgdtmyRowpatchCol) *
                ((patchRow == patchCol) ? 0.5 : 1.0);
          }
          ++myIndex;
        }
        ++patchIndex;
      }

      // h = mydimension * (mydimension + 1) / 2 X worlddimension = 1x2 (curve in 2D), 1x3 (curve in 3D), (3x2 surface
      // in 2D), (3x3 surface in 3D), (6x3 volume in 3D)
      //
      // dfdgdg = gridDimension*(gridDimension + 1) / 2 X worlddimension
      // = 3x2 (curve on surface flat),3x3 (curve on surface in 3D), 6x3 (curve in volume), 6x3 (surface in
      // volume),(3x3 surface in surface in 3D)
      //
      // dgdtSquared = gridDimension * (gridDimension + 1) / 2 X mydimension * (mydimension +
      // 1) / 2 = 3x1 (curve on surface), 6x1 (curve in volume), 6x3 (surface in volume), 3x3 ( surface in surface),
      // 6x6 ( volume in volume)
      Hessian h = transposedView(dgdtSquared) * dfdgdg;

      /* if trimming is enabled the parameter space geometry is potentially non-linear,
       * the resutling Hessian has another contribution due to chain-rule, namely the second derivative of g */
      if constexpr (not Trimmer::template isLocalGeometryLinear<codim>) {
        // const auto dgdtdt = parameterSpaceGeometry->hessian(ouInPatch);
        //
        // assert(TrimmerType::template isLocalGeometryLinear<
        //            codim> && "This can not be checked yet. Check if this works with trimming and then remove
        //            assert");
        // h += transpose(transposedView(dfdg) * dgdtdt);
      }
      return h;
    }

    /* returns the midpoint of the domain */
    [[nodiscard]] LocalCoordinate domainMidPoint() const {
      auto dom = domain();
      LocalCoordinate result{};
      for (int i = 0; i < mydimension; ++i)
        result[i] = dom[i].center();

      return result;
    }

    /* @brief returns the domain, i.e. [0,1] */
    [[nodiscard]] std::array<Utilities::Domain<double>, mydimension> domain() const {
      return {};
    }

    [[nodiscard]] bool affine() const {
      // @todo check when this is true
      return false;
    }

    /**
     * @brief Get the polynomial degree of the patch in each dimension.
     * \remark This does return the degree of the underlying patch and not of the mapping of the local view.
     * @return Array of degrees for each dimension.
     */
    [[nodiscard]] std::array<int, mydimension> degree() const {
      return patchGeometry_->degree();
    }

    friend auto referenceElement(const PatchGeometryLocalView& geo) {
      return referenceElement(*geo.parameterSpaceGeometry);
    }

  private:
    void checkState() const {
      assert(parameterSpaceGeometry && "Bind the local view first!");
    }
    GlobalInParameterSpace globalInParameterSpace(const LocalCoordinate& local) const {
      checkState();
      return parameterSpaceGeometry->global(local);
    }

    [[nodiscard]] JacobianTransposedInParameterSpace jacobianTransposedInParameterSpace(
        const LocalCoordinate& u) const {
      checkState();
      return parameterSpaceGeometry->jacobianTransposed(u);
    }

    ControlPointCoordinateNetType localControlPointNet;
    NurbsLocalView nurbsLocalView_;
    std::array<int, gridDimension> spanIndices_;
    const PatchGeometry* patchGeometry_;
    std::shared_ptr<ParameterSpaceGeometry> parameterSpaceGeometry;
  };
} // namespace GeometryKernel
} // namespace Dune::IGANEW

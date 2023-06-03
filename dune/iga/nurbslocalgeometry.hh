// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "dune/iga/nurbsalgorithms.hh"
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune::IGA {
  /** \brief a geometry implementation for NURBS*/
  template <std::integral auto mydim, std::integral auto dimworld, class GridImpl>
  class NURBSLocalGeometry {
   public:
    /** coordinate type */
    using ctype = typename GridImpl::ctype;

    /** \brief Dimension of the cube element */
    static constexpr std::integral auto mydimension = mydim;

    /** \brief Dimension of the world space that the cube element is embedded in*/
    static constexpr std::integral auto coorddimension = GridImpl::dimensionworld;
    static constexpr std::integral auto griddim        = GridImpl::dimension;

    /** \brief Type used for a vector of element coordinates */
    using LocalCoordinate = Dune::FieldVector<ctype, mydimension>;

    /** \brief Type used for a vector of world coordinates */
    using GlobalCoordinate = Dune::FieldVector<ctype, griddim>;

    /** \brief Type for the transposed Jacobian matrix */
    using JacobianTransposed = Dune::FieldMatrix<ctype, mydimension, coorddimension>;

    /** \brief Type for the transposed inverse Jacobian matrix */
    using JacobianInverseTransposed = Dune::FieldMatrix<ctype, coorddimension, mydimension>;

    using ControlPointType    = typename NURBSPatchData<griddim, dimworld, ctype>::ControlPointType;
    using ControlPointNetType = typename NURBSPatchData<griddim, dimworld, ctype>::ControlPointNetType;

   private:
    /* Helper class to compute a matrix pseudo inverse */
    typedef MultiLinearGeometryTraits<ctype>::MatrixHelper MatrixHelper;

   public:
    /** \brief Constructor from NURBSPatchData and an iterator to a specific knot
     *
     *  \param[in] Patchdata shared pointer to an object where the all the data of the NURBSPatch is stored
     *  \param[in] corner Iterator (for each dimension) to the Knot span where the Geometry object is supposed to
     * operate
     */
    explicit NURBSLocalGeometry(int localSubEntityIndex)
        : localIndexInElement_{localSubEntityIndex}

    {
      assert((localIndexInElement_ < 2 && griddim == 1) || (localIndexInElement_ < 4 && griddim == 2)
             || (localIndexInElement_ < 6 && griddim == 3));
    }

    /** \brief Map the center of the element to the geometry */
    GlobalCoordinate center() const {
      LocalCoordinate localcenter(0.5);
      return global(localcenter);
    }

    [[nodiscard]] double volume() const { return 1; }

    /** \brief Return the number of corners of the element */
    [[nodiscard]] int corners() const { return 1 << mydimension; }

    /** \brief Return world coordinates of the k-th corner of the element */
    GlobalCoordinate corner(int k) const {
      LocalCoordinate localcorner;
      for (size_t i = 0; i < mydimension; i++) {
        localcorner[i] = (k & (1 << i)) ? 1 : 0;
      }
      return global(localcorner);
    }

    [[nodiscard]] bool affine() const { return false; }

    /** \brief Type of the element: a hypercube of the correct dimension */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }

    /** \brief Evaluates the mapping from a subKnotSpan to the global one, e.g. from an edge of an surface it construct
     * the 2d coordinates
     *
     * @param  local local coordinates for each dimension in [0,1]^(dim)
     * @return global coordinates in [0,1]^(griddim)
     */
    GlobalCoordinate global(const LocalCoordinate& local) const {
      const double offset = 0;
      if constexpr (mydimension == 0)
        return (localIndexInElement_ == 0) ? GlobalCoordinate(0) : GlobalCoordinate(1 - offset);
      else if constexpr (mydimension == 1)
        switch (localIndexInElement_) {
          case 0:
            return GlobalCoordinate({0, local[0]});
          case 1:
            return GlobalCoordinate({1 - offset, local[0]});
          case 2:
            return GlobalCoordinate({local[0], 0});
          case 3:
            return GlobalCoordinate({local[0], 1 - offset});
          default:
            __builtin_unreachable();
        }
      else if constexpr (mydimension == 3)
        switch (localIndexInElement_) {
          case 0:
            return GlobalCoordinate({0, local[0], local[1]});
          case 1:
            return GlobalCoordinate({1 - offset, local[0], local[1]});
          case 2:
            return GlobalCoordinate({local[0], 0, local[1]});
          case 3:
            return GlobalCoordinate({local[0], 1 - offset, local[1]});
          case 4:
            return GlobalCoordinate({local[0], local[1], 0});
          case 5:
            return GlobalCoordinate({local[0], local[1], 1 - offset});
          default:
            __builtin_unreachable();
        }
      __builtin_unreachable();
    }

    /** \brief Evaluates the mapping from the grid coordinates to the lcaol one, e.g. from an surface to an edge
     * The local coordinates are obtained from global ones by orthogonal projection, therefore if the global coordinates
     * are [0.5,0.7] and the edge points into the second direction the local coordinate [0.7] is returned
     *
     * @param  global global coordinates for each dimension in [0,1]^(griddim)
     * @return local coordinates in [0,1]^(dim)
     */
    LocalCoordinate local(const GlobalCoordinate& global) const {
      if constexpr (mydimension == 0)
        return (localIndexInElement_ == 0) ? LocalCoordinate(0) : LocalCoordinate(1);
      else if constexpr (mydimension == 1)
        switch (localIndexInElement_) {
          case 0:
          case 1:
            return LocalCoordinate({global[1]});
          case 2:
          case 3:
            return LocalCoordinate({global[0]});
          default:
            __builtin_unreachable();
        }
      else if constexpr (mydimension == 3)
        switch (localIndexInElement_) {
          case 0:
          case 1:
            return LocalCoordinate({global[1], global[2]});
          case 2:
          case 3:
            return LocalCoordinate({global[0], global[2]});
          case 4:
          case 5:
            return LocalCoordinate({global[0], global[1]});
          default:
            __builtin_unreachable();
        }
      __builtin_unreachable();
    }

    /** \brief compute the Jacobian transposed matrix
     *
     *  \param[in] local local coordinates for each dimension
     */
    JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
      if constexpr (mydimension == 0)
        return 0;
      else if constexpr (mydimension == 1)
        switch (localIndexInElement_) {
          case 0:
          case 1:
            return JacobianTransposed({FieldVector<ctype, mydimension>({0, 1})});
          case 2:
          case 3:
            return JacobianTransposed({FieldVector<ctype, mydimension>({1, 0})});
          default:
            __builtin_unreachable();
        }
      else if constexpr (mydimension == 3)
        switch (localIndexInElement_) {
          case 0:
          case 1:
            return LocalCoordinate(
                {FieldVector<ctype, mydimension>({0, 1, 0}), FieldVector<ctype, mydimension>({0, 0, 1})});
          case 2:
          case 3:
            return LocalCoordinate(
                {FieldVector<ctype, mydimension>({1, 0, 0}), FieldVector<ctype, mydimension>({0, 0, 1})});
          case 4:
          case 5:
            return LocalCoordinate(
                {FieldVector<ctype, mydimension>({1, 0, 0}), FieldVector<ctype, mydimension>({0, 1, 0})});
          default:
            __builtin_unreachable();
        }
      __builtin_unreachable();
    }

    /** \brief compute the Jacobian determinant of the mapping
     *
     *  \param[in] local local coordinates for each dimension
     */
    ctype integrationElement(const LocalCoordinate& local) const {
      return MatrixHelper::template sqrtDetAAT<mydimension, coorddimension>(jacobianTransposed(local));
    }

    /** \brief compute the Jacobian inverse transposed matrix
     *
     *  \param[in] local local coordinates for each dimension
     */
    JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
      JacobianInverseTransposed jacobianInverseTransposed1;
      MatrixHelper::template rightInvA<mydimension, coorddimension>(jacobianTransposed(local),
                                                                    jacobianInverseTransposed1);
      return jacobianInverseTransposed1;
    }

    /** \brief computes the unit outward normal
     *
     *  \param[in] local local coordinates for each dimension
     */
    GlobalCoordinate unitNormal(const LocalCoordinate& local) const requires(mydimension == 2) {
      auto J = jacobianTransposed(local);
      auto N = cross(J[0], J[1]);
      return N / N.two_norm();
    }

   private : Dune::Geo::ReferenceElements<ctype, mydimension> referenceElement_;
    int localIndexInElement_;
    bool isTrimmed_{false};
  };

  template <std::integral auto mydim, std::integral auto dimworld, class GridImpl>
  auto referenceElement(const NURBSLocalGeometry<mydim, dimworld, GridImpl>& geo) {
    return Dune::ReferenceElements<typename GridImpl::ctype, mydim>::cube();
  };
}  // namespace Dune::IGA

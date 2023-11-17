// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

/** \file
 * \brief The PatchGridGeometry class and its specializations
 */

#include "enums.hh"

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/iga/hierarchicpatch/nurbspatchgeometry.hh>
namespace Dune::IGANEW {

  template <int mydim, int coorddim, class GridImp>  // requires ( coorddim== GridImp::dimensionworld )
  class PatchGridGeometry : public GeometryDefaultImplementation<mydim, coorddim, GridImp, PatchGridGeometry> {
   public:
    static constexpr int mydimension = mydim;
    static constexpr Trimming trim   = GridImp::trim;

    static constexpr std::integral auto worlddimension = coorddim;
    static constexpr std::integral auto griddim        = GridImp::dimension;
    static constexpr std::integral auto codim          = griddim - mydim;
    using ctype                                        = typename GridImp::ctype;
    using PatchGeometry                                = GeometryKernel::NURBSPatchGeometry<GridImp::dimension, worlddimension, ctype>;
    using LocalCoordinateInPatch                       = typename PatchGeometry::LocalCoordinate;
    using LocalCoordinate                              = FieldVector<ctype, mydimension>;
    using GlobalCoordinate                             = FieldVector<ctype, worlddimension>;
    using JacobianTransposed                           = FieldMatrix<ctype, mydimension, worlddimension>;
    using Hessian                   = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, worlddimension>;
    using Jacobian                  = FieldMatrix<ctype, worlddimension, mydimension>;
    using JacobianInverseTransposed = FieldMatrix<ctype, worlddimension, mydimension>;
    using JacobianInverse           = FieldMatrix<ctype, mydimension, worlddimension>;
    using Volume                    = ctype;

    // The codimension of this entitypointer wrt the host grid
    constexpr static int CodimInHostGrid          = GridImp::HostGridType::dimension - mydim;
    constexpr static int DimensionWorldOfHostGrid = GridImp::HostGridType::dimensionworld;

    // select appropriate hostgrid geometry via typeswitch
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridGeometryType;
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridLocalGeometryType;

    using HostGridGeometry = typename std::conditional<coorddim == DimensionWorldOfHostGrid, HostGridGeometryType,
                                                       HostGridLocalGeometryType>::type;

    //! type of the LocalView of the patch geometry
    using GeometryLocalView =
        typename GeometryKernel::NURBSPatchGeometry<GridImp::dimension, worlddimension, ctype>::template GeometryLocalView<codim, trim>;

    /** constructor from host geometry
     */
    PatchGridGeometry(const HostGridGeometry& hostGeometry, GeometryLocalView&& geometryLocalView)
        : hostGeometry_(hostGeometry), geometryLocalView_(std::forward<GeometryLocalView>(geometryLocalView)) {
      geometryLocalView_.bind(hostGeometry_.impl());
    }

    /** \brief Return the element type identifier
     */
    [[nodiscard]] GeometryType type() const { return geometryLocalView_.type(); }

    // return whether we have an affine mapping
    [[nodiscard]] bool affine() const { return geometryLocalView_.affine(); }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    [[nodiscard]] int corners() const { return hostGeometry_.corners(); }

    //! access to coordinates of corners. Index is the number of the corner
    GlobalCoordinate corner(int i) const { return geometryLocalView_.corner(i); }

    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    GlobalCoordinate global(const LocalCoordinate& local) const { return geometryLocalView_.global(local); }

    /** \brief Return the transposed of the Jacobian
     */
    JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
      return geometryLocalView_.jacobianTransposed(local);
    }

    /** \brief Return the transposed of the Jacobian
 */
    Hessian hessian(const LocalCoordinate& local) const {
      return geometryLocalView_.hessian(local);
    }

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    LocalCoordinate local(const GlobalCoordinate& global) const { return geometryLocalView_.local(global); }

    //! Returns true if the point is in the current element
    bool checkInside(const FieldVector<ctype, mydim>& local) const { return geometryLocalView_.checkInside(local); }

    [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const {
      return geometryLocalView_.integrationElement(local);
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const FieldVector<ctype, mydim>& local) const {
      return geometryLocalView_.jacobianInverseTransposed(local);
    }

   private:
    HostGridGeometry hostGeometry_;
    GeometryLocalView geometryLocalView_{};
  };

}  // namespace Dune::IGANEW

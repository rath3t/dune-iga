// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once
/** \file
 * \brief The PatchGridLocalGeometry class and its specializations
 */

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/geometry.hh>

#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>

namespace Dune::IGANEW {

  template <int mydim, int coorddim, class GridImp>
  class PatchGridLocalGeometry
      : public GeometryDefaultImplementation<mydim, coorddim, GridImp, PatchGridLocalGeometry> {
   public:
    static constexpr int mydimension = mydim;
    using TrimmerType                = typename GridImp::TrimmerType;

    static constexpr std::integral auto coorddimension = coorddim;
    static constexpr std::integral auto griddim        = GridImp::dimension;
    static constexpr std::integral auto codim          = griddim - mydim;
    using ctype                                        = typename GridImp::ctype;
    using PatchGeometry             = GeometryKernel::NURBSPatch<GridImp::dimension, coorddimension, ctype>;
    using LocalCoordinateInPatch    = typename PatchGeometry::LocalCoordinate;
    using LocalCoordinate           = FieldVector<ctype, mydimension>;
    using GlobalCoordinate          = FieldVector<ctype, coorddimension>;
    using JacobianTransposed        = FieldMatrix<ctype, mydimension, coorddimension>;
    using Hessian                   = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, coorddimension>;
    using Jacobian                  = FieldMatrix<ctype, coorddimension, mydimension>;
    using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;
    using JacobianInverse           = FieldMatrix<ctype, mydimension, coorddimension>;
    using Volume                    = ctype;

    // The codimension of this entity wrt the host grid
    constexpr static int CodimInHostGrid          = GridImp::ParameterSpaceGrid::dimension - mydim;
    constexpr static int DimensionWorldOfHostGrid = GridImp::ParameterSpaceGrid::dimensionworld;

    // select appropriate hostgrid geometry via typeswitch
    typedef typename GridImp::ParameterSpaceGrid::Traits::template Codim<CodimInHostGrid>::LocalGeometry
        HostGridLocalGeometryType;

    // using LocalNURBSGeometry = GeometryKernel::NURBSPatch<mydim, coorddim, ctype>;
    using LocalHostGeometry =  HostGridLocalGeometryType;
    //! type of the LocalView of the patch geometry
    using GeometryLocalView =
        typename GeometryKernel::NURBSPatch<GridImp::dimension, coorddimension,
                                            ctype>::template GeometryLocalView<codim, TrimmerType>;

    /** constructor from host geometry
     */
    // PatchGridLocalGeometry(const GeometryKernel::NURBSPatch<mydim, coorddim, ctype>& localPatchGeometry)
    //     : localPatchGeometry_(localPatchGeometry) {}

    explicit PatchGridLocalGeometry(const LocalHostGeometry& hostGeometry)
    : hostGeometry_(hostGeometry) {
    }

    // PatchGridLocalGeometry(const HostGridGeometry& hostGeometry) : hostGeometry_(hostGeometry) {}

    /** \brief Return the element type identifier
     */
    [[nodiscard]] GeometryType type() const { return hostGeometry_.type(); }

    // return whether we have an affine mapping
    [[nodiscard]] bool affine() const { return hostGeometry_.affine(); }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    [[nodiscard]] int corners() const { return hostGeometry_.corners(); }

    //! access to coordinates of corners. Index is the number of the corner
    GlobalCoordinate corner(int i) const { return hostGeometry_.corner(i); }

    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    GlobalCoordinate global(const LocalCoordinate& local) const { return hostGeometry_.global(local); }

    /** \brief Return the transposed of the Jacobian
     */
    JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
      return hostGeometry_.jacobianTransposed(local);
    }

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    LocalCoordinate local(const GlobalCoordinate& global) const { return hostGeometry_.local(global); }

    //! Returns true if the point is in the current element
    bool checkInside(const FieldVector<ctype, mydim>& local) const { return hostGeometry_.checkInside(local); }

    [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const {
      return hostGeometry_.integrationElement(local);
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const FieldVector<ctype, mydim>& local) const {
      return hostGeometry_.jacobianInverseTransposed(local);
    }

   private:
    LocalHostGeometry hostGeometry_;
    // GeometryLocalView geometryLocalView_{};
  };

}  // namespace Dune::IGANEW
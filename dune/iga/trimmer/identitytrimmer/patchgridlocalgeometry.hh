// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once
/** \file
 * @brief The PatchGridLocalGeometry class and its specializations
 */

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>

namespace Dune::IGANEW {

/**
 * \brief This class describes entities in the 0..1 space if untrimmed this is in the cube,
 * for the untrimmed case this is simply a wrapper around the Yasp grid local geometry
 * \tparam mydim
 * \tparam coorddim
 * \tparam GridImp
 */
template <int mydim, int coorddim, class GridImp>
class PatchGridLocalGeometry : public GeometryDefaultImplementation<mydim, coorddim, GridImp, PatchGridLocalGeometry>
{
public:
  static constexpr int mydimension = mydim;
  using Trimmer                    = typename GridImp::Trimmer;

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

  using LocalGeometry = typename Trimmer::template Codim<codim>::LocalGeometry;
  // using ParameterSpaceGeometry = typename Trimmer::template LocalGeometry<codim>;

  explicit PatchGridLocalGeometry(const LocalGeometry& localGeometry)
      : localGeometry_(localGeometry) {}

  // PatchGridLocalGeometry(const HostGridGeometry& hostGeometry) : hostGeometry_(hostGeometry) {}

  /** @brief Return the element type identifier
   */
  [[nodiscard]] GeometryType type() const {
    return localGeometry_.type();
  }

  // return whether we have an affine mapping
  [[nodiscard]] bool affine() const {
    return localGeometry_.affine();
  }

  // return the number of corners of this element. Corners are numbered 0...n-1
  [[nodiscard]] int corners() const {
    return localGeometry_.corners();
  }

  // access to coordinates of corners. Index is the number of the corner
  GlobalCoordinate corner(int i) const {
    return localGeometry_.corner(i);
  }

  /** @brief Maps a local coordinate within reference element to
   * global coordinate in element  */
  GlobalCoordinate global(const LocalCoordinate& local) const {
    return localGeometry_.global(local);
  }

  /** @brief Return the transposed of the Jacobian
   */
  JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
    return localGeometry_.jacobianTransposed(local);
  }

  /** @brief Maps a global coordinate within the element to a
   * local coordinate in its reference element */
  LocalCoordinate local(const GlobalCoordinate& global) const {
    return localGeometry_.local(global);
  }

  // Returns true if the point is in the current element
  bool checkInside(const FieldVector<ctype, mydim>& local) const {
    return localGeometry_.checkInside(local);
  }

  [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const {
    return localGeometry_.integrationElement(local);
  }

  // The Jacobian matrix of the mapping from the reference element to this element
  [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const FieldVector<ctype, mydim>& local) const {
    return localGeometry_.jacobianInverseTransposed(local);
  }

private:
  LocalGeometry localGeometry_;
};

} // namespace Dune::IGANEW

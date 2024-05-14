// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

/** \file
 * @brief The PatchGridGeometry class and its specializations
 */


#include <dune/common/fmatrix.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>
namespace Dune::IGA {

template <int mydim, int coorddim, class GridImp>
class PatchGridGeometry : public GeometryDefaultImplementation<mydim, coorddim, GridImp, PatchGridGeometry>
{
public:
  static constexpr int mydimension = mydim;
  // static constexpr Trimming trim   = GridImp::trim;
  using Trimmer = typename GridImp::Trimmer;

  static constexpr std::integral auto worlddimension = coorddim;
  static constexpr std::integral auto griddim        = GridImp::dimension;
  static constexpr std::integral auto codim          = griddim - mydim;
  using ctype                                        = typename GridImp::ctype;
  using PatchGeometry             = GeometryKernel::NURBSPatch<GridImp::dimension, worlddimension, ctype>;
  using LocalCoordinateInPatch    = typename PatchGeometry::LocalCoordinate;
  using LocalCoordinate           = FieldVector<ctype, mydimension>;
  using GlobalCoordinate          = FieldVector<ctype, worlddimension>;
  using JacobianTransposed        = FieldMatrix<ctype, mydimension, worlddimension>;
  using Hessian                   = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, worlddimension>;
  using Jacobian                  = FieldMatrix<ctype, worlddimension, mydimension>;
  using JacobianInverseTransposed = FieldMatrix<ctype, worlddimension, mydimension>;
  using JacobianInverse           = FieldMatrix<ctype, mydimension, worlddimension>;
  using Volume                    = ctype;

  // The geometry in the parameterspace, i.e. in the knotspan domain
  using ParameterSpaceGeometry = typename Trimmer::template Codim<codim>::LocalParameterSpaceGeometry;

  // using LocalGeometryInParameterSpace = typename ReferenceElementType::template Codim<CodimInHostGrid>::Geometry;
  //  type of the LocalView of the patch geometry
  using GeometryLocalView = typename GeometryKernel::NURBSPatch<GridImp::dimension, worlddimension,
                                                                ctype>::template GeometryLocalView<codim, Trimmer>;

  /** constructor from host geometry */
  PatchGridGeometry(const ParameterSpaceGeometry& localGeometry, GeometryLocalView&& geometryLocalView)
      : localGeometry_(localGeometry),
        geometryLocalView_(std::forward<GeometryLocalView>(geometryLocalView)) {
    geometryLocalView_.bind(localGeometry);
  }

  /** @brief Return the element type identifier
   */
  [[nodiscard]] GeometryType type() const {
    return geometryLocalView_.type();
  }

  // return whether we have an affine mapping
  [[nodiscard]] bool affine() const {
    return geometryLocalView_.affine();
  }

  // return the number of corners of this element. Corners are numbered 0...n-1
  [[nodiscard]] int corners() const {
    return geometryLocalView_.corners();
  }

  // access to coordinates of corners. Index is the number of the corner
  [[nodiscard]] GlobalCoordinate corner(int i) const {
    return geometryLocalView_.corner(i);
  }

  /** @brief Maps a local coordinate within reference element to
   * global coordinate in element  */
  [[nodiscard]] GlobalCoordinate global(const LocalCoordinate& local) const {
    return geometryLocalView_.global(local);
  }

  /** @brief Return the transposed of the Jacobian
   */
  [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
    return geometryLocalView_.jacobianTransposed(local);
  }

  /** @brief Return the Hessian */
  [[nodiscard]] Hessian hessian(const LocalCoordinate& local) const {
    return geometryLocalView_.hessian(local);
  }

  /** @brief Maps a global coordinate to a
   * local coordinate in its reference element */
  [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
    return geometryLocalView_.local(global);
  }

  // Returns true if the point is in the current element
  [[nodiscard]] bool checkInside(const FieldVector<ctype, mydim>& local) const {
    return geometryLocalView_.checkInside(local);
  }

  [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const {
    return geometryLocalView_.integrationElement(local);
  }

  // The Jacobian matrix of the mapping from the reference element to this element
  [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const FieldVector<ctype, mydim>& local) const {
    return geometryLocalView_.jacobianInverseTransposed(local);
  }

  auto zeroFirstAndSecondDerivativeOfPosition(const LocalCoordinate& u) const {
    return geometryLocalView_.zeroFirstAndSecondDerivativeOfPosition(u);
  }

  auto secondFundamentalForm(const LocalCoordinate& local) const {
    return geometryLocalView_.secondFundamentalForm(local);
  }

private:
  ParameterSpaceGeometry localGeometry_;
  GeometryLocalView geometryLocalView_{};
};

} // namespace Dune::IGA

// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

namespace Dune::IGANEW::DefaultTrim {

enum class LocalGeometryTag
{
  InParameterSpace,
  InReferenceElement
};

template <int mydim, int coorddim, class GridImp, LocalGeometryTag localGeometryTag>
class TrimmedLocalGeometryImpl
{
public:
  using ctype = typename GridImp::ctype;

  static constexpr int mydimension = mydim;
  using Trimmer                    = typename GridImp::Trimmer;

  static constexpr int coorddimension = coorddim;
  static constexpr int codim          = coorddimension - mydimension;
  using PatchGeometry                 = GeometryKernel::NURBSPatch<mydimension, coorddimension, ctype>;
  using LocalCoordinateInPatch        = typename PatchGeometry::LocalCoordinate;
  using LocalCoordinate               = FieldVector<ctype, mydimension>;
  using GlobalCoordinate              = FieldVector<ctype, coorddimension>;
  using JacobianTransposed            = FieldMatrix<ctype, mydimension, coorddimension>;
  using Hessian                       = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, coorddimension>;
  using Jacobian                      = FieldMatrix<ctype, coorddimension, mydimension>;
  using JacobianInverseTransposed     = FieldMatrix<ctype, coorddimension, mydimension>;
  using JacobianInverse               = FieldMatrix<ctype, mydimension, coorddimension>;
  using Volume                        = ctype;

  //! type of the LocalView of the patch geometry
  using GeometryLocalView = typename GeometryKernel::NURBSPatch<mydimension, coorddimension,
                                                                ctype>::template GeometryLocalView<codim, Trimmer>;

  /** constructor from host geometry  */
  TrimmedLocalGeometryImpl() = default;
  explicit TrimmedLocalGeometryImpl(const GeometryKernel::NURBSPatch<mydimension, coorddimension, ctype>& patchGeometry)
      : patchGeometry{&patchGeometry} {}

  /** @brief Return the element type identifier
   */
  [[nodiscard]] GeometryType type() const { return GeometryTypes::none(mydimension); }

  // return whether we have an affine mapping (true for straight lines)
  [[nodiscard]] bool affine() const {
    // handy if the trims are straight lines for other algorithms
    if constexpr (codim == 0)
      return true;
    else
      return patchGeometry.degree()[0] == 1;
  }

  //! return the number of corners of this element. Corners are numbered 0...n-1
  [[nodiscard]] int corners() const { return patchGeometry.corners(); }

  GlobalCoordinate center() const { return patchGeometry.center(); }

  //! access to coordinates of corners. Index is the number of the corner
  GlobalCoordinate corner(int i) const { return patchGeometry.corner(i); }

  /** @brief Maps a local coordinate within reference element to
   * global coordinate in element  */
  GlobalCoordinate global(const LocalCoordinate& local) const { return patchGeometry.global(local); }

  /** @brief Return the transposed of the Jacobian
   */
  JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
    return patchGeometry.jacobianTransposed(local);
  }

  /** @brief Maps a global coordinate within the element to a
   * local coordinate in its reference element */
  LocalCoordinate local(const GlobalCoordinate& global) const { return patchGeometry.local(global); }

  //! Returns true if the point is in the current element
  // @todo
  bool checkInside(const LocalCoordinate& local) const { return true; }

  [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const { return patchGeometry.volume(); }

  //! The Jacobian matrix of the mapping from the reference element to this element
  // @todo not yet implemented
  [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
    DUNE_THROW(Dune::NotImplemented, "jacobianInverseTransposed() not yet implemented");
    return JacobianInverseTransposed{};
  }

private:
  GeometryKernel::NURBSPatch<mydimension, coorddimension, ctype> patchGeometry;
};
} // namespace Dune::IGANEW::DefaultTrim

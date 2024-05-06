// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/indextransformations.hh>

namespace Dune::IGANEW::DefaultTrim {

enum class LocalGeometryTag
{
  InParameterSpace,
  InReferenceElement
};

template <int mydim, int coorddim, class GridImp, LocalGeometryTag localGeometryTag>
class TrimmedLocalGeometryImpl
{
};

/* elements */
template <int coorddim, class GridImp, LocalGeometryTag localGeometryTag>
class TrimmedLocalGeometryImpl<2, coorddim, GridImp, localGeometryTag>
{
public:
  using ctype = typename GridImp::ctype;

  static constexpr int mydimension = 2;
  using Trimmer                    = typename GridImp::Trimmer;

  static constexpr int coorddimension = coorddim;
  static constexpr int codim          = coorddimension - mydimension;
  using PatchGeometry                 = GeometryKernel::NURBSPatch<mydimension, coorddimension, ctype>;
  using LocalCoordinateInPatch        = typename PatchGeometry::LocalCoordinate;

  using HostGeometry = typename Trimmer::TrimmerTraits::template Codim<0>::HostParameterSpaceGridEntity::Geometry;
  using TrimData     = typename Trimmer::ElementTrimData;

  using LocalCoordinate           = FieldVector<ctype, mydimension>;
  using GlobalCoordinate          = FieldVector<ctype, coorddimension>;
  using JacobianTransposed        = FieldMatrix<ctype, mydimension, coorddimension>;
  using Hessian                   = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, coorddimension>;
  using Jacobian                  = FieldMatrix<ctype, coorddimension, mydimension>;
  using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;
  using JacobianInverse           = FieldMatrix<ctype, mydimension, coorddimension>;
  using Volume                    = ctype;

  // type of the LocalView of the patch geometry
  using GeometryLocalView = typename GeometryKernel::NURBSPatch<mydimension, coorddimension,
                                                                ctype>::template GeometryLocalView<codim, Trimmer>;

  /** constructor from host geometry  */
  TrimmedLocalGeometryImpl() = default;
  explicit TrimmedLocalGeometryImpl(const HostGeometry& hostGeometry, const TrimData& trimData)
      : hostGeometry_{hostGeometry},
        trimData_{trimData} {}
  bool operator==(const TrimmedLocalGeometryImpl& b) const {
    return b.hostGeometry_ == this->hostGeometry_;
  }

  /** @brief Return the element type identifier
   */
  [[nodiscard]] GeometryType type() const {
    return GeometryTypes::none(mydimension);
  }

  // return whether we have an affine mapping (true for straight lines)
  [[nodiscard]] bool affine() const {
    return true;
  }

  // return the number of corners of this element. Corners are numbered 0...n-1
  [[nodiscard]] int corners() const {
    return trimData_.size(2);
  }

  // @todo write a test for that
  GlobalCoordinate center() const {
    // average of corners
    GlobalCoordinate c{0, 0};
    for (auto i : Dune::range(corners()))
      c += corner(i);

    return c / corners();
  }

  // access to coordinates of corners. Index is the number of the corner
  GlobalCoordinate corner(int i) const {
    auto vData = trimData_.vertex(i);
    if (vData.isHost)
      return hostGeometry_.corner(Transformations::mapToDune(2, vData.idx));

    return vData.geometry.value();
  }

  /** @brief Maps a local coordinate within reference element to
   * global coordinate in element  */
  GlobalCoordinate global(const LocalCoordinate& local) const {
    return hostGeometry_.global(local);
  }

  /** @brief Return the transposed of the Jacobian
   */
  JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
    return hostGeometry_.jacobianTransposed(local);
  }

  /** @brief Maps a global coordinate within the element to a
   * local coordinate in its reference element */
  LocalCoordinate local(const GlobalCoordinate& global) const {
    return hostGeometry_.local(global);
  }

  // Returns true if the point is in the current element
  // @todo
  bool checkInside(const LocalCoordinate& local) const {
    return true;
  }

  // @todo
  [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const {
    return hostGeometry_.volume();
  }

  [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
    return hostGeometry_.jacobianInverseTransposed(local);
  }

private:
  HostGeometry hostGeometry_;
  TrimData trimData_;
};

template <int coorddim, class GridImp, LocalGeometryTag localGeometryTag>
class TrimmedLocalGeometryImpl<1, coorddim, GridImp, localGeometryTag>
{
public:
  using ctype = typename GridImp::ctype;

  static constexpr int mydimension = 1;
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

  // type of the LocalView of the patch geometry
  using GeometryLocalView = typename GeometryKernel::NURBSPatch<mydimension, coorddimension,
                                                                ctype>::template GeometryLocalView<codim, Trimmer>;

  /** constructor from host geometry  */
  TrimmedLocalGeometryImpl() = default;
  explicit TrimmedLocalGeometryImpl(const GeometryKernel::NURBSPatch<mydimension, coorddimension, ctype>& patchGeometry)
      : patchGeometry{patchGeometry} {}

  bool operator==(const TrimmedLocalGeometryImpl& b) const {
    return Dune::FloatCmp::eq(b.patchGeometry.corner(0) == this->patchGeometry.corner(0)) and
           Dune::FloatCmp::eq(b.patchGeometry.corner(1) == this->patchGeometry.corner(1));
  };

  /** @brief Return the element type identifier
   */
  [[nodiscard]] GeometryType type() const {
    return GeometryTypes::cube(mydimension);
  }

  // return whether we have an affine mapping (true for straight lines)
  [[nodiscard]] bool affine() const {
    // handy if the trims are straight lines for other algorithms
    if constexpr (codim == 0)
      return true;
    else
      return patchGeometry.degree()[0] == 1;
  }

  // return the number of corners of this element. Corners are numbered 0...n-1
  [[nodiscard]] int corners() const {
    return patchGeometry.corners();
  }

  GlobalCoordinate center() const {
    return patchGeometry.center();
  }

  // access to coordinates of corners. Index is the number of the corner
  GlobalCoordinate corner(int i) const {
    return patchGeometry.corner(i);
  }

  /** @brief Maps a local coordinate within reference element to
   * global coordinate in element  */
  GlobalCoordinate global(const LocalCoordinate& local) const {
    return patchGeometry.global(local);
  }

  /** @brief Return the transposed of the Jacobian
   */
  JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
    return patchGeometry.jacobianTransposed(local);
  }

  /** @brief Maps a global coordinate within the element to a
   * local coordinate in its reference element */
  LocalCoordinate local(const GlobalCoordinate& global) const {
    return patchGeometry.local(global);
  }

  // Returns true if the point is in the current element
  // @todo
  bool checkInside(const LocalCoordinate& local) const {
    return true;
  }

  [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const {
    return patchGeometry.volume();
  }

  // The Jacobian matrix of the mapping from the reference element to this element
  // @todo not yet implemented
  [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
    return patchGeometry.jacobianInverseTransposed(local);
  }

private:
  GeometryKernel::NURBSPatch<mydimension, coorddimension, ctype> patchGeometry;
};

// Template specialization for trimmed vertices
template <int coorddim, class GridImp, LocalGeometryTag localGeometryTag>
class TrimmedLocalGeometryImpl<0, coorddim, GridImp, localGeometryTag>
{
public:
  using ctype = typename GridImp::ctype;

  static constexpr int mydimension = 0;
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

  TrimmedLocalGeometryImpl() = default;
  explicit TrimmedLocalGeometryImpl(const FieldVector<ctype, coorddimension>& pos)
      : pos_{pos} {}
  bool operator==(const TrimmedLocalGeometryImpl& b) const {
    return Dune::FloatCmp::eq(b.pos_, this->pos_);
  }

  // @todo it is unclear to me what the correct bahviour for a vertex is for some of these methods

  /** @brief Return the element type identifier
   */
  [[nodiscard]] GeometryType type() const {
    return GeometryTypes::cube(mydimension);
  }

  // return whether we have an affine mapping (true for vertices??)
  [[nodiscard]] bool affine() const {
    return true;
  }

  // return the number of corners of this element. Corners are numbered 0...n-1
  [[nodiscard]] int corners() const {
    return 1;
  }

  [[nodiscard]] GlobalCoordinate center() const {
    return pos_;
  }

  // access to coordinates of corners. Index is the number of the corner
  [[nodiscard]] GlobalCoordinate corner(int i) const {
    return pos_;
  }

  GlobalCoordinate global(const LocalCoordinate& local) const {
    return pos_;
  }

  /** @brief Return the transposed of the Jacobian
   */
  JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
    return JacobianTransposed{};
  }

  // Yaspgrid returns an empty {} for vertex.local
  LocalCoordinate local(const GlobalCoordinate& global) const {
    return LocalCoordinate{};
  }

  bool checkInside(const LocalCoordinate& local) const {
    return true;
  }

  [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const {
    return 1;
  }

  [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
    return JacobianInverseTransposed{};
  }

private:
  FieldVector<ctype, coorddimension> pos_{};
};

} // namespace Dune::IGANEW::DefaultTrim

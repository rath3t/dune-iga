// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

#include "patchgridentity.hh"

/** \file
 * @brief The PatchGridLeafIntersection and PatchGridLevelIntersection classes
 */

namespace Dune::IGA {

// External forward declarations
template <class Grid>
struct HostGridAccess;

namespace Impl {

  enum class IntersectionType
  {
    Level,
    Leaf
  };

  template <class GridImp, IntersectionType type_>
  struct PatchGridIntersectionImpl
  {
    constexpr static int dim      = GridImp::dimension;
    constexpr static int mydim    = GridImp::dimension - 1;
    constexpr static int dimworld = GridImp::dimensionworld;

    using Geometry      = typename GridImp::template Codim<1>::Geometry;
    using LocalGeometry = typename GridImp::template Codim<1>::LocalGeometry;
    using Entity        = typename GridImp::template Codim<0>::Entity;

    using ctype           = typename GridImp::ctype;
    using NormalVector    = FieldVector<ctype, dimworld>;
    using LocalCoordinate = FieldVector<ctype, mydim>;

    using ParameterSpaceIntersection =
        std::conditional_t<type_ == IntersectionType::Leaf,
                           typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLeafIntersection,
                           typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLevelIntersection>;

    PatchGridIntersectionImpl() = default;

    PatchGridIntersectionImpl(const GridImp* parameterSpaceGrid, const ParameterSpaceIntersection& hostIntersection)
        : patchGrid_(parameterSpaceGrid),
          parameterSpaceIntersection_(hostIntersection) {}

    PatchGridIntersectionImpl(const GridImp* parameterSpaceGrid, ParameterSpaceIntersection&& hostIntersection)
        : patchGrid_(parameterSpaceGrid),
          parameterSpaceIntersection_(std::move(hostIntersection)) {}

    bool equals(const PatchGridIntersectionImpl& other) const {
      return parameterSpaceIntersection_ == other.parameterSpaceIntersection_;
    }

    // return Entity on the inside of this intersection
    // (that is the Entity where we started this Iterator)
    Entity inside() const {
      return PatchGridEntity<0, dim, GridImp>(patchGrid_, parameterSpaceIntersection_.inside());
    }

    // return Entity on the outside of this intersection
    // (that is the neighboring Entity)
    Entity outside() const {
      return PatchGridEntity<0, dim, GridImp>(patchGrid_, parameterSpaceIntersection_.outside());
    }

    // return true if intersection is with boundary.
    bool boundary() const {
      return parameterSpaceIntersection_.boundary();
    }

    /** @brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal() const {
      LocalCoordinate localcenter;
      localcenter = 0.5;
      return this->unitOuterNormal(localcenter);
    }

    // return true if across the edge an neighbor on this level exists
    bool neighbor() const {
      return parameterSpaceIntersection_.neighbor();
    }

    // return the boundary segment index
    size_t boundarySegmentIndex() const {
      return parameterSpaceIntersection_.boundarySegmentIndex();
    }

    // Return true if this is a conforming intersection
    bool conforming() const {
      return parameterSpaceIntersection_.conforming();
    }

    // Geometry type of an intersection
    GeometryType type() const {
      return parameterSpaceIntersection_.type();
    }

    // intersection of codimension 1 of this neighbor with element where
    // iteration started.
    // Here returned element is in LOCAL coordinates of the element
    // where iteration started.
    LocalGeometry geometryInInside() const {
      return LocalGeometry(typename LocalGeometry::Implementation(parameterSpaceIntersection_.geometryInInside()));
    }

    // intersection of codimension 1 of this neighbor with element where iteration started.
    // Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside() const {
      return LocalGeometry(typename LocalGeometry::Implementation(parameterSpaceIntersection_.geometryInOutside()));
    }

    // intersection of codimension 1 of this neighbor with element where iteration started.
    // Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry() const {
      // @todo trim this will be wrong as soon as the intersection geometry has a special geoemtry
      auto geo =
          typename Geometry::Implementation(parameterSpaceIntersection_.geometry(),
                                            patchGridGeometry().template localView<1, typename GridImp::Trimmer>());
      return Geometry(geo);
    }

    // local number of codim 1 entity in self where intersection is contained in
    int indexInInside() const {
      return parameterSpaceIntersection_.indexInInside();
    }

    // local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside() const {
      return parameterSpaceIntersection_.indexInOutside();
    }

    // return outer normal
    FieldVector<ctype, GridImp::dimensionworld> outerNormal(const LocalCoordinate& local) const {
      FieldMatrix<ctype, dimworld, dim> J =
          inside().geometry().jacobianInverseTransposed(geometryInInside().global(local));
      FieldVector<ctype, dimworld> res;
      J.mv(parameterSpaceIntersection_.outerNormal(local), res);
      return res;
    }

    // return outer normal multiplied by the integration element
    FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal(const LocalCoordinate& local) const {
      const ctype detJ                 = this->geometry().integrationElement(local);
      FieldVector<ctype, dimworld> res = unitOuterNormal(local);
      res *= detJ;
      return res;
    }

    // return unit outer normal
    FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal(const LocalCoordinate& local) const {
      auto Jinv = inside().geometry().jacobianInverseTransposed(geometryInInside().global(local));
      FieldVector<ctype, dimworld> res;
      Jinv.mv(parameterSpaceIntersection_.unitOuterNormal(local), res);
      res /= res.two_norm();
      return res;
    }

    bool isTrimmed() const {
      if constexpr (requires { parameterSpaceIntersection_.isTrimmed(); })
        return parameterSpaceIntersection_.isTrimmed();
      return false;
    }

  private:
    auto& patchGridGeometry() const {
      if constexpr (type_ == IntersectionType::Leaf)
        return patchGrid_->patchGeometryAtBack();
      else
        return patchGrid_->patchGeometry(inside().level());
    }
    const GridImp* patchGrid_{};
    ParameterSpaceIntersection parameterSpaceIntersection_{};
  };
} // namespace Impl

/** @brief An intersection with a leaf neighbor element
 * \ingroup PatchGrid
 * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
 * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
 * These neighbors are accessed via a IntersectionIterator. This allows the implement
 * non-matching meshes. The number of neighbors may be different from the number
 * of an element!
 */
template <class GridImp>
class PatchGridLeafIntersection
{
  friend typename GridImp::Traits::LeafIntersectionIterator;

  friend struct HostGridAccess<std::remove_const_t<GridImp>>;

  using Implementation                 = Impl::PatchGridIntersectionImpl<GridImp, Impl::IntersectionType::Leaf>;
  using ParameterSpaceLeafIntersection = typename Implementation::ParameterSpaceIntersection;

public:
  constexpr static int dim      = Implementation::dimension;
  constexpr static int mydim    = Implementation::dimension - 1;
  constexpr static int dimworld = Implementation::dimensionworld;

  using Geometry      = typename Implementation::Geometry;
  using LocalGeometry = typename Implementation::LocalGeometry;
  using Entity        = typename Implementation::Entity;

  using ctype           = typename Implementation::ctype;
  using NormalVector    = typename Implementation::NormalVector;
  using LocalCoordinate = typename Implementation::LocalCoordinate;

  PatchGridLeafIntersection() = default;

  PatchGridLeafIntersection(const GridImp* parameterSpaceGrid, const ParameterSpaceLeafIntersection& hostIntersection)
      : impl_(parameterSpaceGrid, hostIntersection) {}

  PatchGridLeafIntersection(const GridImp* parameterSpaceGrid, ParameterSpaceLeafIntersection&& hostIntersection)
      : impl_(parameterSpaceGrid, std::move(hostIntersection)) {}

  bool equals(const PatchGridLeafIntersection& other) const {
    return impl_.equals(other.impl_);
  }

  Entity inside() const {
    return impl_.inside();
  }

  Entity outside() const {
    return impl_.outside();
  }

  bool boundary() const {
    return impl_.boundary();
  }

  bool neighbor() const {
    return impl_.neighbor();
  }

  size_t boundarySegmentIndex() const {
    return impl_.boundarySegmentIndex();
  }

  bool conforming() const {
    return impl_.conforming();
  }

  GeometryType type() const {
    return impl_.type();
  }

  LocalGeometry geometryInInside() const {
    return impl_.geometryInInside();
  }

  LocalGeometry geometryInOutside() const {
    return impl_.geometryInOutside();
  }

  Geometry geometry() const {
    return impl_.geometry();
  }

  // local number of codim 1 entity in self where intersection is contained in
  int indexInInside() const {
    return impl_.indexInInside();
  }

  int indexInOutside() const {
    return impl_.indexInOutside();
  }

  NormalVector outerNormal(const LocalCoordinate& local) const {
    return impl_.outerNormal(local);
  }

  NormalVector integrationOuterNormal(const LocalCoordinate& local) const {
    return impl_.integrationOuterNormal(local);
  }

  NormalVector centerUnitOuterNormal() const {
    return impl_.centerUnitOuterNormal();
  }

  // return unit outer normal
  NormalVector unitOuterNormal(const LocalCoordinate& local) const {
    return impl_.unitOuterNormal(local);
  }

  bool isTrimmed() const {
    return impl_.isTrimmed();
  }

private:
  Implementation impl_;
};

template <class GridImp>
class PatchGridLevelIntersection
{
  friend typename GridImp::Traits::LevelIntersectionIterator;

  friend struct HostGridAccess<std::remove_const_t<GridImp>>;

  using Implementation                  = Impl::PatchGridIntersectionImpl<GridImp, Impl::IntersectionType::Level>;
  using ParameterSpaceLevelIntersection = typename Implementation::ParameterSpaceIntersection;

public:
  constexpr static int dim      = Implementation::dimension;
  constexpr static int mydim    = Implementation::dimension - 1;
  constexpr static int dimworld = Implementation::dimensionworld;

  using Geometry      = typename Implementation::Geometry;
  using LocalGeometry = typename Implementation::LocalGeometry;
  using Entity        = typename Implementation::Entity;

  using ctype           = typename Implementation::ctype;
  using NormalVector    = typename Implementation::NormalVector;
  using LocalCoordinate = typename Implementation::LocalCoordinate;

  PatchGridLevelIntersection() = default;

  PatchGridLevelIntersection(const GridImp* identityGrid, const ParameterSpaceLevelIntersection& hostIntersection)
      : impl_(identityGrid, hostIntersection) {}

  PatchGridLevelIntersection(const GridImp* identityGrid, ParameterSpaceLevelIntersection&& hostIntersection)
      : impl_(identityGrid, std::move(hostIntersection)) {}

  bool equals(const PatchGridLevelIntersection& other) const {
    return impl_.equals(other.impl_);
  }

  Entity inside() const {
    return impl_.inside();
  }

  Entity outside() const {
    return impl_.outside();
  }

  bool boundary() const {
    return impl_.boundary();
  }

  bool neighbor() const {
    return impl_.neighbor();
  }

  size_t boundarySegmentIndex() const {
    return impl_.boundarySegmentIndex();
  }

  bool conforming() const {
    return impl_.conforming();
  }

  GeometryType type() const {
    return impl_.type();
  }

  LocalGeometry geometryInInside() const {
    return impl_.geometryInInside();
  }

  LocalGeometry geometryInOutside() const {
    return impl_.geometryInOutside();
  }

  Geometry geometry() const {
    return impl_.geometry();
  }

  int indexInInside() const {
    return impl_.indexInInside();
  }

  int indexInOutside() const {
    return impl_.indexInOutside();
  }

  NormalVector outerNormal(const LocalCoordinate& local) const {
    return impl_.outerNormal(local);
  }

  NormalVector integrationOuterNormal(const LocalCoordinate& local) const {
    return impl_.integrationOuterNormal(local);
  }

  NormalVector centerUnitOuterNormal() const {
    return impl_.centerUnitOuterNormal();
  }

  // return unit outer normal
  NormalVector unitOuterNormal(const LocalCoordinate& local) const {
    return impl_.unitOuterNormal(local);
  }

  bool isTrimmed() const {
    return impl_.isTrimmed();
  }

private:
  Implementation impl_;
};

} // namespace Dune::IGA

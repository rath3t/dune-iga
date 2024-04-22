// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

/** \file
 * @brief The TrimmedPatchGridLeafIntersection and PatchGridLevelIntersection classes
 */

namespace Dune::IGANEW::Trim {

// External forward declarations
template <class Grid>
struct HostGridAccess;

/** @brief An intersection with a leaf neighbor element
 * \ingroup PatchGrid
 * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
 * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
 * These neighbors are accessed via a IntersectionIterator. This allows the implement
 * non-matching meshes. The number of neighbors may be different from the number
 * of an element!
 */
template <typename TrimmerType_, class... Implementations>
class IntersectionVariant
{
  using ctype                         = std::common_type_t<typename Implementations::ctype...>;
  using FirstElement                  = std::tuple_element_t<0, std::tuple<Implementations...>>;
  static constexpr int mydimension    = FirstElement::mydimension;
  static constexpr int coorddimension = FirstElement::coorddimension;

  using LocalCoordinate = FieldVector<ctype, mydimension>;

public:
  auto visit(auto&& lambda) const {
    return std::visit(lambda, impl_);
  }

  template <class Implementation>
  IntersectionVariant(const Implementation& impl)
      : impl_(impl) {
  }

  IntersectionVariant()                                 = default;
  IntersectionVariant(const IntersectionVariant& other) = default;
  template <class Implementation>
  requires(!std::is_same_v<Implementation, IntersectionVariant>)
  IntersectionVariant& operator=(const Implementation& impl) {
    impl_ = impl;
    return *this;
  };
  IntersectionVariant(IntersectionVariant&& other) noexcept            = default;
  IntersectionVariant& operator=(const IntersectionVariant& other)     = default;
  IntersectionVariant& operator=(IntersectionVariant&& other) noexcept = default;

  // bool equals(const TrimmedPatchGridLeafIntersection& other) const { return parameterSpaceIntersection ==
  // other.parameterSpaceIntersection; }

  //! return Entity on the inside of this intersection
  //! (that is the Entity where we started this Iterator)
  auto inside() const {
    return visit([](const auto& impl) { return impl.inside(); });
  }

  //! return Entity on the outside of this intersection
  //! (that is the neighboring Entity)
  auto outside() const {
    return visit([](const auto& impl) { return impl.outside(); });
  }

  //! return true if intersection is with boundary.
  [[nodiscard]] bool boundary() const {
    return visit([](const auto& impl) { return impl.boundary(); });
  }

  /** @brief Return unit outer normal (length == 1)
   *
   *   The returned vector is the normal at the center() of the
   *     intersection's geometry.
   *       It is scaled to have unit length. */
  auto centerUnitOuterNormal() const {
    return visit([](const auto& impl) { return impl.centerUnitOuterNormal(); });
  }

  //! return true if across the edge an neighbor on this level exists
  bool neighbor() const {
    return visit([](const auto& impl) { return impl.neighbor(); });
  }

  //! return the boundary segment index
  size_t boundarySegmentIndex() const {
    return visit([](const auto& impl) { return impl.boundarySegmentIndex(); });
  }

  //! Return true if this is a conforming intersection
  bool conforming() const {
    return visit([](const auto& impl) { return impl.conforming(); });
  }

  //! Geometry type of an intersection
  GeometryType type() const {
    return visit([](const auto& impl) { return impl.type(); });
  }

  auto geometryInInside() const {
    return visit([](const auto& impl) { return impl.geometryInInside(); });
  }

  //! intersection of codimension 1 of this neighbor with element where iteration started.
  //! Here returned element is in LOCAL coordinates of neighbor
  auto geometryInOutside() const {
    return visit([](const auto& impl) { return impl.geometryInOutside(); });
  }

  //! intersection of codimension 1 of this neighbor with element where iteration started.
  //! Here returned element is in GLOBAL coordinates of the element where iteration started.
  auto geometry() const {
    return visit([](const auto& impl) { return impl.geometry(); });
  }

  //! local number of codim 1 entity in self where intersection is contained in
  int indexInInside() const {
    return visit([](const auto& impl) { return impl.indexInInside(); });
  }

  //! local number of codim 1 entity in neighbor where intersection is contained
  int indexInOutside() const {
    return visit([](const auto& impl) { return impl.indexInOutside(); });
  }

  //! return outer normal
  auto outerNormal(const LocalCoordinate& local) const {
    return visit([&](const auto& impl) { return impl.outerNormal(local); });
  }

  //! return outer normal multiplied by the integration element
  auto integrationOuterNormal(const LocalCoordinate& local) const {
    return visit([&](const auto& impl) { return impl.integrationOuterNormal(local); });
  }

  //! return unit outer normal
  auto unitOuterNormal(const LocalCoordinate& local) const {
    return visit([&](const auto& impl) { return impl.unitOuterNormal(local); });
  }

private:
  std::variant<Implementations...> impl_;
};
//
// template <class GridImp>
// class PatchGridLevelIntersection {
//   friend  typename GridImp::Traits::LevelIntersectionIterator;
//
//   friend struct HostGridAccess<typename std::remove_const<GridImp>::type>;
//
//   constexpr static int dim   = GridImp::dimension;
//   constexpr static int mydim = GridImp::dimension - 1;
//
//   constexpr static int dimworld = GridImp::dimensionworld;
//
//   using Trimmer = typename GridImp::Trimmer;
//
//   // The type used to store coordinates
//   typedef typename GridImp::ctype ctype;
//
//   typedef typename GridImp::GridFamily::LevelIntersection HostLevelIntersection;
//
//   using LocalCoordinate = FieldVector<ctype, mydim>;
//
//   using MatrixHelper = typename MultiLinearGeometryTraits<double>::MatrixHelper;
//
//  public:
//   typedef typename GridImp::template Codim<1>::Geometry Geometry;
//   typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
//   typedef typename GridImp::template Codim<0>::Entity Entity;
//   typedef FieldVector<ctype, dimworld> NormalVector;
//
//   PatchGridLevelIntersection() = default;
//
//   PatchGridLevelIntersection(const GridImp* identityGrid, const HostLevelIntersection& hostIntersection)
//       : patchGrid_(identityGrid), parameterSpaceIntersection(hostIntersection) {}
//
//   PatchGridLevelIntersection(const GridImp* identityGrid, HostLevelIntersection&& hostIntersection)
//       : patchGrid_(identityGrid), parameterSpaceIntersection(std::move(hostIntersection)) {}
//
//   [[nodiscard]] bool equals(const PatchGridLevelIntersection& other) const {
//     return parameterSpaceIntersection == other.parameterSpaceIntersection;
//   }
//
//   //! return Entity on the inside of this intersection
//   //! (that is the Entity where we started this Iterator)
//   [[nodiscard]] Entity inside() const {
//     return PatchGridEntity<0, dim, GridImp>(patchGrid_, parameterSpaceIntersection.inside());
//   }
//
//   //! return Entity on the outside of this intersection
//   //! (that is the neighboring Entity)
//   [[nodiscard]] Entity outside() const {
//     return PatchGridEntity<0, dim, GridImp>(patchGrid_, parameterSpaceIntersection.outside());
//   }
//
//   /** @brief return true if intersection is with boundary.
//    */
//   [[nodiscard]] bool boundary() const { return parameterSpaceIntersection.boundary(); }
//
//   /** @brief Return unit outer normal (length == 1)
//    *
//    *   The returned vector is the normal at the center() of the
//    *     intersection's geometry.
//    *       It is scaled to have unit length. */
//   NormalVector centerUnitOuterNormal() const {
//     LocalCoordinate localcenter;
//     localcenter = 0.5;
//     return this->unitOuterNormal(localcenter);
//   }
//
//   //! return true if across the edge an neighbor on this level exists
//   [[nodiscard]] bool neighbor() const { return parameterSpaceIntersection.neighbor(); }
//
//   //! return the boundary segment index
//   [[nodiscard]] size_t boundarySegmentIndex() const { return parameterSpaceIntersection.boundarySegmentIndex(); }
//
//   //! Return true if this is a conforming intersection
//   [[nodiscard]] bool conforming() const { return parameterSpaceIntersection.conforming(); }
//
//   //! Geometry type of an intersection
//   [[nodiscard]] GeometryType type() const { return parameterSpaceIntersection.type(); }
//
//   //! intersection of codimension 1 of this neighbor with element where
//   //! iteration started.
//   //! Here returned element is in LOCAL coordinates of the element
//   //! where iteration started.
//   [[nodiscard]] LocalGeometry geometryInInside() const {
//     auto localGeometry = typename
//     LocalGeometry::Implementation::LocalGeometry(parameterSpaceIntersection.geometryInInside()); return
//     LocalGeometry(typename LocalGeometry::Implementation(localGeometry));
//   }
//
//   //! intersection of codimension 1 of this neighbor with element where iteration started.
//   //! Here returned element is in LOCAL coordinates of neighbor
//   [[nodiscard]] LocalGeometry geometryInOutside() const {
//     auto localGeometry = typename
//     LocalGeometry::Implementation::LocalGeometry(parameterSpaceIntersection.geometryInOutside()); return
//     LocalGeometry(typename LocalGeometry::Implementation(localGeometry));
//   }
//
//   //! intersection of codimension 1 of this neighbor with element where iteration started.
//   //! Here returned element is in GLOBAL coordinates of the element where iteration started.
//   [[nodiscard]] Geometry geometry() const {
//     // @todo trim does this make sense?
//     auto geo = typename Geometry::Implementation(
//         parameterSpaceIntersection.geometry(),
//         patchGrid_->patchGeometries_[inside().level()].template localView<1, Trimmer>());
//     return Geometry(geo);
//   }
//
//   //! local number of codim 1 entity in self where intersection is contained in
//   [[nodiscard]] int indexInInside() const { return parameterSpaceIntersection.indexInInside(); }
//
//   //! local number of codim 1 entity in neighbor where intersection is contained
//   [[nodiscard]] int indexInOutside() const { return parameterSpaceIntersection.indexInOutside(); }
//
//   //! return outer normal
//   [[nodiscard]] FieldVector<ctype, dimworld> outerNormal(const LocalCoordinate& local) const {
//     const auto globalInPatch            = geometryInInside().global(local);
//     FieldMatrix<ctype, dimworld, dim> J = inside().geometry().jacobianInverseTransposed(globalInPatch);
//     FieldVector<ctype, dimworld> res;
//     auto refElement                                           = referenceElement(inside().geometry());
//     const int indexInInside                                   = this->indexInInside();
//     const typename LocalGeometry::GlobalCoordinate& refNormal = refElement.integrationOuterNormal(indexInInside);
//     const auto refNormal2                                     = parameterSpaceIntersection.outerNormal(local);
//     // std::cout<<"refNormal"<<refNormal2<<"refNormal"<<refNormal<<std::endl;
//     // std::cout<<"J"<<J<<std::endl;
//     J.mv(refNormal, res);
//     // std::cout<<"res"<<res<<std::endl;
//     return res;
//   }
//
//   //! return outer normal multiplied by the integration element
//   [[nodiscard]] FieldVector<ctype, dimworld> integrationOuterNormal(const LocalCoordinate& local) const {
//     const ctype detJ                 = this->geometry().integrationElement(local);
//     FieldVector<ctype, dimworld> res = unitOuterNormal(local);
//
//     res *= detJ;
//     return res;
//   }
//
//   //! return unit outer normal
//   [[nodiscard]] FieldVector<ctype, dimworld> unitOuterNormal(const LocalCoordinate& local) const {
//     FieldVector<ctype, dimworld> res = outerNormal(local);
//     res /= res.two_norm();
//     return res;
//   }
//
//  private:
//   const GridImp* patchGrid_;
//
//   HostLevelIntersection parameterSpaceIntersection;
// };

} // namespace Dune::IGANEW::Trim

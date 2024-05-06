// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

#include "patchgridintersections.hh"

#include <dune/iga/geometrykernel/nurbspatchtransform.hh>

/** \file
 * @brief The TrimmedPatchGridLeafIntersection and TrimmedLevelIntersection classes
 */

namespace Dune::IGANEW::DefaultTrim {

// External forward declarations
template <class Grid>
struct HostGridAccess;

namespace Impl {

  enum class IntersectionType
  {
    Level,
    Leaf
  };

  namespace IntersectionTraits {

    template <class GridImp, IntersectionType type_>
    using HostIntersection = std::conditional_t<type_ == IntersectionType::Leaf,
                                                typename GridImp::Trimmer::TrimmerTraits::HostLeafIntersection,
                                                typename GridImp::Trimmer::TrimmerTraits::HostLevelIntersection>;

    template <class GridImp>
    using Geometry = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry;

    template <class GridImp>
    using LocalGeometry = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalGeometry;

    template <class GridImp>
    using ParameterSpaceGridEntity =
        typename GridImp::Trimmer::TrimmerTraits::template Codim<0>::ParameterSpaceGridEntity;

  } // namespace IntersectionTraits

  template <class GridImp, IntersectionType type_>
  struct TrimmedIntersectionImpl
  {
    constexpr static int dim      = GridImp::dimension;
    constexpr static int mydim    = GridImp::dimension - 1;
    constexpr static int dimworld = GridImp::dimension;
    using ctype                   = typename GridImp::ctype;
    using IdType                  = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;
    using NormalVector            = FieldVector<ctype, dim>;
    using LocalCoordinate         = FieldVector<ctype, mydim>;
    using EdgeInfo                = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::EntityInfo;

    using IntersectionGeometry     = GeometryKernel::NURBSPatch<mydim, dim, ctype>;
    using HostLeafIntersection     = IntersectionTraits::HostIntersection<GridImp, type_>;
    using ParameterSpaceGridEntity = IntersectionTraits::ParameterSpaceGridEntity<GridImp>;
    using LocalGeometry            = IntersectionTraits::LocalGeometry<GridImp>;
    using Geometry                 = IntersectionTraits::Geometry<GridImp>;
    using TrimmedLocalGeometry = typename GridImp::GridFamily::TrimmerTraits::template Codim<1>::TrimmedLocalGeometry;
    using TrimmedGeometry =
        typename GridImp::GridFamily::TrimmerTraits::template Codim<1>::TrimmedParameterSpaceGeometry;

    TrimmedIntersectionImpl() = default;

    TrimmedIntersectionImpl(const GridImp* patchGrid, const IdType& insideElementId, const EdgeInfo& edgeInfo,
                            int indexInInside)
        : patchGrid_(patchGrid),
          edgeInfo_(edgeInfo),
          insideElementId_(insideElementId),
          geo_(edgeInfo_.intersectionGeometry()),
          indexInInside_(indexInInside) {
      assert(geo_.domain()[0].isUnitDomain());
    }

    ParameterSpaceGridEntity inside() const {
      return patchGrid_->trimmer().entityContainer_.template entity<0>(insideElementId_, level());
    }

    ParameterSpaceGridEntity outside() const {
      DUNE_THROW(GridError, "No Outside Entities for trimmed Intersections");
    }

    bool boundary() const {
      return true;
    }

    bool neighbor() const {
      return false;
    }

    size_t boundarySegmentIndex() const {
      return 0;
    }

    // A boundary has to be conforming
    bool conforming() const {
      return true;
    }

    GeometryType type() const {
      return GeometryTypes::none(mydim);
    }

    LocalGeometry geometryInInside() const {
      auto localGeo = GeometryKernel::transformToSpan(geo_, inside().getHostEntity().geometry());
      return LocalGeometry(TrimmedLocalGeometry(localGeo));
    }

    LocalGeometry geometryInOutside() const {
      DUNE_THROW(GridError, "No Outside Entities for trimmed Intersections");
    }

    Geometry geometry() const {
      return Geometry(TrimmedGeometry(geo_));
    }

    int indexInInside() const {
      return indexInInside_;
    }

    int indexInOutside() const {
      DUNE_THROW(GridError, "No Outside Entities for trimmed Intersections");
    }

    NormalVector outerNormal(const LocalCoordinate& local) const {
      auto jacobian = geo_.jacobianTransposed(local)[0];
      return {jacobian[1], jacobian[0]};
    }

    NormalVector integrationOuterNormal(const LocalCoordinate& local) const {
      return unitOuterNormal(local) * geo_.integrationElement(local);
    }

    NormalVector unitOuterNormal(const LocalCoordinate& local) const {
      auto N = outerNormal(local);
      return N / N.two_norm();
    }

    NormalVector centerUnitOuterNormal() const {
      return unitOuterNormal(0.5);
    }

    bool equals(const TrimmedIntersectionImpl& other) const {
      return geo_ == other.geo_ and indexInInside_ == other.indexInInside_;
    }

  private:
    int level() const {
      if constexpr (type_ == IntersectionType::Leaf)
        return patchGrid_->maxLevel();
      else
        return edgeInfo_.lvl;
    }

    const GridImp* patchGrid_{};
    EdgeInfo edgeInfo_;
    IdType insideElementId_;
    IntersectionGeometry geo_;
    int indexInInside_;
  };

  template <class GridImp, IntersectionType type_>
  struct TrimmedHostIntersectionImpl
  {
    constexpr static int dim      = GridImp::dimension;
    constexpr static int mydim    = GridImp::dimension - 1;
    constexpr static int dimworld = GridImp::dimension;
    using ctype                   = typename GridImp::ctype;
    using IdType                  = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;
    using NormalVector            = FieldVector<ctype, dim>;
    using LocalCoordinate         = FieldVector<ctype, mydim>;
    using EdgeInfo                = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::EntityInfo;

    using IntersectionGeometry     = GeometryKernel::NURBSPatch<mydim, dim, ctype>;
    using HostLeafIntersection     = IntersectionTraits::HostIntersection<GridImp, type_>;
    using ParameterSpaceGridEntity = IntersectionTraits::ParameterSpaceGridEntity<GridImp>;
    using LocalGeometry            = IntersectionTraits::LocalGeometry<GridImp>;
    using Geometry                 = IntersectionTraits::Geometry<GridImp>;
    using TrimmedLocalGeometry = typename GridImp::GridFamily::TrimmerTraits::template Codim<1>::TrimmedLocalGeometry;
    using TrimmedGeometry =
        typename GridImp::GridFamily::TrimmerTraits::template Codim<1>::TrimmedParameterSpaceGeometry;
    TrimmedHostIntersectionImpl() = default;

    TrimmedHostIntersectionImpl(const GridImp* patchGrid, const HostLeafIntersection& hostIntersection,
                                const EdgeInfo& edgeinfo, int indexInInside)
        : patchGrid_(patchGrid),
          hostIntersection_(hostIntersection),
          edgeInfo_(edgeinfo),
          geo_(edgeInfo_.intersectionGeometry()),
          indexInInside_(indexInInside) {
      assert(geo_.domain()[0].isUnitDomain());
    }

    ParameterSpaceGridEntity inside() const {
      return patchGrid_->trimmer().entityContainer_.template entity<0>(insideElementId(), level());
    }

    ParameterSpaceGridEntity outside() const {
      return patchGrid_->trimmer().entityContainer_.template entity<0>(outsideElementId(), level());
    }

    bool boundary() const {
      return hostIntersection_.boundary();
    }

    bool neighbor() const {
      return hostIntersection_.neighbor();
    }

    size_t boundarySegmentIndex() const {
      return 0;
    }

    bool conforming() const {
      return hostIntersection_.conforming();
    }

    // This is up to debate, if its not GeometryTypes::none(mydim);
    GeometryType type() const {
      return hostIntersection_.type();
    }

    LocalGeometry geometryInInside() const {
      auto localGeo = GeometryKernel::transformToSpan(geo_, inside().getHostEntity().geometry());
      return LocalGeometry(TrimmedLocalGeometry(localGeo));
    }

    LocalGeometry geometryInOutside() const {
      auto outsideGeo      = edgeInfo_.intersectionGeometry();
      auto localOutsideGeo = GeometryKernel::transformToSpan(outsideGeo, outside().getHostEntity().geometry());
      return LocalGeometry(TrimmedLocalGeometry(localOutsideGeo));
    }

    Geometry geometry() const {
      return Geometry(TrimmedGeometry(geo_));
    }

    int indexInInside() const {
      return indexInInside_;
    }

    int indexInOutside() const {
      if (patchGrid_->trimmer().entityContainer_.isElementTrimmed(outsideElementId()))
        return patchGrid_->trimmer().entityContainer_.outsideIntersectionIndex(insideElementId(), outsideElementId(),
                                                                               indexInInside());
      return hostIntersection_.indexInOutside();
    }

    // Its fine to use the hostIntersection as we are in the parameter space, and host intersections are always constant
    // straight lines
    NormalVector outerNormal(const LocalCoordinate& local) const {
      return hostIntersection_.outerNormal(local);
    }

    NormalVector integrationOuterNormal(const LocalCoordinate& local) const {
      return unitOuterNormal(local) * geo_.integrationElement(local);
    }

    NormalVector unitOuterNormal(const LocalCoordinate& local) const {
      auto N = outerNormal(local);
      return N / N.two_norm();
    }

    NormalVector centerUnitOuterNormal() const {
      return unitOuterNormal(0.5);
    }

    bool equals(const TrimmedHostIntersectionImpl& other) const {
      return hostIntersection_ == other.hostIntersection_ and indexInInside_ == other.indexInInside_;
    }

  private:
    int level() const {
      if constexpr (type_ == IntersectionType::Leaf)
        return patchGrid_->maxLevel();
      else
        return hostIntersection_.inside().level();
    }
    IdType insideElementId() const {
      auto hostId = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.inside());
      return {.entityIdType = IdType::EntityIdType::host, .id = hostId};
    }
    IdType outsideElementId() const {
      auto hostId = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.outside());
      return {.entityIdType = IdType::EntityIdType::host, .id = hostId};
    }

    const GridImp* patchGrid_{};
    HostLeafIntersection hostIntersection_;
    EdgeInfo edgeInfo_;
    IntersectionGeometry geo_;
    int indexInInside_;
  };

  template <class GridImp, IntersectionType type_>
  struct HostIntersectionImpl
  {
    constexpr static int dim      = GridImp::dimension;
    constexpr static int mydim    = GridImp::dimension - 1;
    constexpr static int dimworld = GridImp::dimension;
    using ctype                   = typename GridImp::ctype;
    using IdType                  = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;
    using NormalVector            = FieldVector<ctype, dim>;
    using LocalCoordinate         = FieldVector<ctype, mydim>;

    using IntersectionGeometry     = GeometryKernel::NURBSPatch<mydim, dim, ctype>;
    using HostLeafIntersection     = IntersectionTraits::HostIntersection<GridImp, type_>;
    using ParameterSpaceGridEntity = IntersectionTraits::ParameterSpaceGridEntity<GridImp>;
    using LocalGeometry            = IntersectionTraits::LocalGeometry<GridImp>;
    using Geometry                 = IntersectionTraits::Geometry<GridImp>;

    HostIntersectionImpl() = default;

    HostIntersectionImpl(const GridImp* patchGrid, const HostLeafIntersection& hostIntersection,
                         std::optional<int> indexInInside = std::nullopt)
        : patchGrid_(patchGrid),
          hostIntersection_(hostIntersection),
          indexInInside_(indexInInside) {}

    ParameterSpaceGridEntity inside() const {
      return patchGrid_->trimmer().entityContainer_.template entity<0>(insideElementId(), level());
    }

    ParameterSpaceGridEntity outside() const {
      return patchGrid_->trimmer().entityContainer_.template entity<0>(outsideElementId(), level());
    }

    bool boundary() const {
      return hostIntersection_.boundary();
    }

    bool neighbor() const {
      return hostIntersection_.neighbor();
    }

    // return the boundary segment index
    size_t boundarySegmentIndex() const {
      return 0;
      // This is not implmented in SubGrid
    }

    // Return true if this is a conforming intersection
    bool conforming() const {
      return hostIntersection_.conforming();
    }

    GeometryType type() const {
      return hostIntersection_.type();
    }

    LocalGeometry geometryInInside() const {
      return hostIntersection_.geometryInInside();
    }

    LocalGeometry geometryInOutside() const {
      return hostIntersection_.geometryInOutside();
    }

    Geometry geometry() const {
      return hostIntersection_.geometry();
    }

    int indexInInside() const {
      return indexInInside_.value_or(hostIntersection_.indexInInside());
    }

    int indexInOutside() const {
      if (patchGrid_->trimmer().entityContainer_.isElementTrimmed(outsideElementId()))
        return patchGrid_->trimmer().entityContainer_.outsideIntersectionIndex(insideElementId(), outsideElementId(),
                                                                               indexInInside());
      return hostIntersection_.indexInOutside();
    }

    NormalVector outerNormal(const LocalCoordinate& local) const {
      return hostIntersection_.outerNormal(local);
    }

    NormalVector integrationOuterNormal(const LocalCoordinate& local) const {
      return hostIntersection_.integrationOuterNormal(local);
    }

    NormalVector unitOuterNormal(const LocalCoordinate& local) const {
      return hostIntersection_.integrationOuterNormal(local);
    }

    NormalVector centerUnitOuterNormal() const {
      return hostIntersection_.centerUnitOuterNormal();
    }

    bool equals(const HostIntersectionImpl& other) const {
      return hostIntersection_ == other.hostIntersection_;
    }

  private:
    int level() const {
      if constexpr (type_ == IntersectionType::Leaf)
        return patchGrid_->maxLevel();
      else
        return hostIntersection_.inside().level();
    }

    IdType insideElementId() const {
      auto hostId = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.inside());
      return {.entityIdType = IdType::EntityIdType::host, .id = hostId};
    }

    IdType outsideElementId() const {
      auto hostId = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.outside());
      return {.entityIdType = IdType::EntityIdType::host, .id = hostId};
    }

    const GridImp* patchGrid_{};
    IntersectionGeometry geo;
    HostLeafIntersection hostIntersection_;
    std::optional<int> indexInInside_;
  };

  template <class GridImp, IntersectionType type_>
  struct IntersectionVariant
  {
    template <class Implementation>
    explicit IntersectionVariant(const Implementation& impl)
        : impl_(impl) {}
    IntersectionVariant() = default;

    auto visit(auto&& lambda) const {
      return std::visit(lambda, impl_);
    }

    auto inside() const {
      return visit([](const auto& impl) { return impl.inside(); });
    }

    auto outside() const {
      return visit([](const auto& impl) { return impl.outside(); });
    }

    bool boundary() const {
      return visit([](const auto& impl) { return impl.boundary(); });
    }

    bool neighbor() const {
      return visit([](const auto& impl) { return impl.neighbor(); });
    }

    size_t boundarySegmentIndex() const {
      return visit([](const auto& impl) { return impl.boundarySegmentIndex(); });
    }

    bool conforming() const {
      return visit([](const auto& impl) { return impl.conforming(); });
    }

    auto type() const {
      return visit([](const auto& impl) { return impl.type(); });
    }

    auto geometryInInside() const {
      return visit([](const auto& impl) { return impl.geometryInInside(); });
    }

    auto geometryInOutside() const {
      return visit([](const auto& impl) { return impl.geometryInOutside(); });
    }

    auto geometry() const {
      return visit([](const auto& impl) { return impl.geometry(); });
    }

    int indexInInside() const {
      return visit([](const auto& impl) { return impl.indexInInside(); });
    }

    int indexInOutside() const {
      return visit([](const auto& impl) { return impl.indexInOutside(); });
    }

    auto outerNormal(const auto& local) const {
      return visit([&](const auto& impl) { return impl.outerNormal(local); });
    }

    auto integrationOuterNormal(const auto& local) const {
      return visit([&](const auto& impl) { return impl.integrationOuterNormal(local); });
    }

    auto unitOuterNormal(const auto& local) const {
      return visit([&](const auto& impl) { return impl.unitOuterNormal(local); });
    }

    auto centerUnitOuterNormal() const {
      return visit([](const auto& impl) { return impl.centerUnitOuterNormal(); });
    }

    bool equals(const IntersectionVariant& other) const {
      return visit([&]<typename T>(const T& impl) { return impl.equals(std::get<T>(other.impl_)); });
    }

  private:
    std::variant<TrimmedIntersectionImpl<GridImp, type_>, TrimmedHostIntersectionImpl<GridImp, type_>,
                 HostIntersectionImpl<GridImp, type_>>
        impl_{};
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
class TrimmedLeafIntersection
{
  friend typename GridImp::Traits::LeafIntersectionIterator;

  friend struct HostGridAccess<std::remove_const_t<GridImp>>;

  constexpr static int dim      = GridImp::dimension;
  constexpr static int mydim    = GridImp::dimension - 1;
  constexpr static int dimworld = GridImp::dimension;

  using Trimmer = typename GridImp::Trimmer;

  using HostLeafIntersection = typename GridImp::Trimmer::TrimmerTraits::HostLeafIntersection;
  using EdgeInfo             = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::EntityInfo;

public:
  // The type used to store coordinates
  typedef typename GridImp::ctype ctype;
  using LocalCoordinate = FieldVector<ctype, mydim>;

  typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry Geometry;
  typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalGeometry LocalGeometry;
  typedef FieldVector<ctype, dim> NormalVector;
  using IdType = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;

  using ParameterSpaceGridEntity =
      typename GridImp::Trimmer::TrimmerTraits::template Codim<0>::ParameterSpaceGridEntity;

  using IntersectionGeometry = GeometryKernel::NURBSPatch<mydim, dim, ctype>;
  using IntersectionImpl     = Impl::IntersectionVariant<GridImp, Impl::IntersectionType::Leaf>;

  TrimmedLeafIntersection() = default;

  TrimmedLeafIntersection(const GridImp* patchGrid, const HostLeafIntersection& hostIntersection,
                          std::optional<int> indexInInside = std::nullopt)
      : underlying_{Impl::HostIntersectionImpl<GridImp, Impl::IntersectionType::Leaf>(patchGrid, hostIntersection,
                                                                                      indexInInside)} {}

  // Trimmed Host Intersection
  TrimmedLeafIntersection(const GridImp* patchGrid, const HostLeafIntersection& hostIntersection,
                          const EdgeInfo& edgeInfo, int indexInInside)
      : underlying_{Impl::TrimmedHostIntersectionImpl<GridImp, Impl::IntersectionType::Leaf>(
            patchGrid, hostIntersection, edgeInfo, indexInInside)} {}

  // Trimmed Intersection
  TrimmedLeafIntersection(const GridImp* patchGrid, const IdType& insideElementId, const EdgeInfo& edgeInfo,
                          int indexInInside)
      : underlying_{Impl::TrimmedIntersectionImpl<GridImp, Impl::IntersectionType::Leaf>(patchGrid, insideElementId,
                                                                                         edgeInfo, indexInInside)} {}

  // TrimmedLeafIntersection(const GridImp* parameterSpaceGrid, HostLeafIntersection&& hostIntersection)
  //     : patchGrid_(parameterSpaceGrid),
  //       hostIntersection_{hostIntersection} {}

  bool operator==(const TrimmedLeafIntersection& other) const {
    return underlying_.equals(other.underlying_);
  }

  /**
   *
   * @return the inside entity
   */
  ParameterSpaceGridEntity inside() const {
    return underlying_.inside();
  }

  /**
   * \brief Entity on the outside of this intersection (that is the neighboring Entity)
   * \return neighbor entity
   */
  ParameterSpaceGridEntity outside() const {
    return underlying_.outside();
  }

  /**
   * \brief True if intersection is part of the domain boundary
   */
  [[nodiscard]] bool boundary() const {
    return underlying_.boundary();
  }

  // return true if across the edge an neighbor on this level exists
  bool neighbor() const {
    return underlying_.neighbor();
  }

  // return the boundary segment index
  /**
   * \brief Boundary segment index, not implemented in Dune::SubGrid
   * @return boundary segment index
   */
  size_t boundarySegmentIndex() const {
    return underlying_.boundarySegmentIndex();
  }

  /**
   * \brief Returns whether the intersection is conforming, i.e., whether it is equivalent to an entire element facet
   * shared by the inside and outside elements
   * @return true if this is a conforming intersection
   */
  bool conforming() const {
    return underlying_.conforming();
  }

  /**
   *\brief type of the reference element for the intersection, i.e., the domain for its parametrization
   * @return  Geometry type of an intersection
   */
  GeometryType type() const {
    return underlying_.type();
  }

  /**
   * @return Geometry of this intersection in local coordinates of the inside entity
   */
  LocalGeometry geometryInInside() const {
    return underlying_.geometryInInside();
  }

  /**
   * @return Geometry of this intersection in local coordinates of the outside entity
   */
  LocalGeometry geometryInOutside() const {
    return underlying_.geometryInOutside();
  }

  /**
   * @return Geometry of the intersection in the grid world space, here parameterspace geometry
   */
  Geometry geometry() const {
    return underlying_.geometry();
  }

  /**
   * \brief Local index of facet in the inside (codim 1) entity that contains the intersection
   * @return local index
   */
  int indexInInside() const {
    return underlying_.indexInInside();
  }

  /**
   * \brief Local index of facet in the outside (codim 1) entity that contains the intersection
   * @return local index
   */
  int indexInOutside() const {
    return underlying_.indexInOutside();
  }

  // return outer normal
  FieldVector<ctype, dim> outerNormal(const LocalCoordinate& local) const {
    return underlying_.outerNormal(local);
  }

  // return outer normal multiplied by the integration element
  FieldVector<ctype, dim> integrationOuterNormal(const LocalCoordinate& local) const {
    return underlying_.integrationOuterNormal(local);
  }

  // return unit outer normal
  FieldVector<ctype, dim> unitOuterNormal(const LocalCoordinate& local) const {
    return underlying_.unitOuterNormal(local);
  }

  /** @brief Return unit outer normal (length == 1)
   *
   *   The returned vector is the normal at the center() of the
   *     intersection's geometry.
   *       It is scaled to have unit length. */
  NormalVector centerUnitOuterNormal() const {
    return underlying_.centerUnitOuterNormal();
  }

private:
  IntersectionImpl underlying_{};
};

template <class GridImp>
class TrimmedLevelIntersection
{
  friend typename GridImp::Traits::LevelIntersectionIterator;

  friend struct HostGridAccess<std::remove_const_t<GridImp>>;

  constexpr static int dim   = GridImp::dimension;
  constexpr static int mydim = GridImp::dimension - 1;

  constexpr static int dimworld = dim;

  using Trimmer = typename GridImp::Trimmer;

  using ParameterSpaceGridEntity =
      typename GridImp::Trimmer::TrimmerTraits::template Codim<0>::ParameterSpaceGridEntity;
  typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry Geometry;
  typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalGeometry LocalGeometry;

  using EdgeInfo = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::EntityInfo;
  using IdType   = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;

  using HostLevelIntersection = typename GridImp::Trimmer::TrimmerTraits::HostLevelIntersection;

  using MatrixHelper = MultiLinearGeometryTraits<double>::MatrixHelper;

public:
  using IntersectionImpl = Impl::IntersectionVariant<GridImp, Impl::IntersectionType::Level>;

  typedef typename GridImp::ctype ctype;
  using LocalCoordinate = FieldVector<ctype, mydim>;

  typedef FieldVector<ctype, dimworld> NormalVector;

  TrimmedLevelIntersection() = default;

  TrimmedLevelIntersection(const GridImp* parameterSpaceGrid, const HostLevelIntersection& hostIntersection,
                           std::optional<int> indexInInside = std::nullopt)
      : underlying_{Impl::HostIntersectionImpl<GridImp, Impl::IntersectionType::Level>(
            parameterSpaceGrid, hostIntersection, indexInInside)} {}

  // Trimmed Host Intersection
  TrimmedLevelIntersection(const GridImp* patchGrid, const HostLevelIntersection& hostIntersection,
                           const EdgeInfo& edgeInfo, int indexInInside)
      : underlying_{Impl::TrimmedHostIntersectionImpl<GridImp, Impl::IntersectionType::Level>(
            patchGrid, hostIntersection, edgeInfo, indexInInside)} {
    assert(edgeInfo.isTrimmed());
  }

  // Trimmed Intersection
  TrimmedLevelIntersection(const GridImp* patchGrid, const IdType& insideElementId, const EdgeInfo& edgeInfo,
                           int indexInInside)
      : underlying_{Impl::TrimmedIntersectionImpl<GridImp, Impl::IntersectionType::Level>(patchGrid, insideElementId,
                                                                                          edgeInfo, indexInInside)} {}

  // TrimmedLevelIntersection(const GridImp* patchGrid, HostLevelIntersection&& hostIntersection)
  //     : patchGrid_(patchGrid),
  //       hostIntersection_{hostIntersection} {}

  bool operator==(const TrimmedLevelIntersection& other) const {
    return underlying_.equals(other.underlying_);
  }

  // returns the inside entity
  ParameterSpaceGridEntity inside() const {
    return underlying_.inside();
  }

  // return Entity on the outside of this intersection
  // (that is the neighboring Entity)
  ParameterSpaceGridEntity outside() const {
    return underlying_.outside();
  }

  // return true if intersection is with boundary.
  [[nodiscard]] bool boundary() const {
    return underlying_.boundary();
  }

  /** @brief Return unit outer normal (length == 1)
   *
   *   The returned vector is the normal at the center() of the
   *     intersection's geometry.
   *       It is scaled to have unit length. */
  NormalVector centerUnitOuterNormal() const {
    return underlying_.centerUnitOuterNormal();
  }

  // return true if across the edge an neighbor on this level exists
  bool neighbor() const {
    return underlying_.neighbor();
  }

  // return the boundary segment index
  size_t boundarySegmentIndex() const {
    return underlying_.boundarySegmentIndex();
  }

  // Return true if this is a conforming intersection
  bool conforming() const {
    return underlying_.conforming();
  }

  // Geometry type of an intersection
  GeometryType type() const {
    return underlying_.type();
  }

  LocalGeometry geometryInInside() const {
    return underlying_.geometryInInside();
  }

  LocalGeometry geometryInOutside() const {
    return underlying_.geometryInOutside();
  }

  Geometry geometry() const {
    return underlying_.geometry();
  }

  // local number of codim 1 entity in self where intersection is contained in
  int indexInInside() const {
    return underlying_.indexInInside();
  }

  // local number of codim 1 entity in neighbor where intersection is contained
  int indexInOutside() const {
    return underlying_.indexInOutside();
  }

  // return outer normal
  FieldVector<ctype, dim> outerNormal(const LocalCoordinate& local) const {
    return underlying_.outerNormal(local);
  }

  // return outer normal multiplied by the integration element
  FieldVector<ctype, dim> integrationOuterNormal(const LocalCoordinate& local) const {
    return underlying_.integrationOuterNormal(local);
  }

  // return unit outer normal
  FieldVector<ctype, dim> unitOuterNormal(const LocalCoordinate& local) const {
    return underlying_.unitOuterNormal(local);
  }

private:
  IntersectionImpl underlying_{};
};

} // namespace Dune::IGANEW::DefaultTrim

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

#include <dune/grid/common/intersection.hh>
#include <dune/iga/hierarchicpatch/patchgridentity.hh>
#include <dune/iga/hierarchicpatch/patchgridintersections.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/indextransformations.hh>

/** \file
 * @brief The PatchGridLeafIntersectionIterator and PatchGridLevelIntersectionIterator classes
 */

namespace Dune::IGANEW::DefaultTrim {

namespace Impl {

  enum class PositionToken
  {
    Begin,
    End
  };

  enum class IntersectionIteratorType
  {
    Level,
    Leaf
  };

  namespace IntersectionIteratorTraits {
    template <class GridImp, IntersectionIteratorType type_>
    using Intersection =
        std::conditional_t<type_ == IntersectionIteratorType::Leaf, typename GridImp::Traits::LeafIntersection,
                           typename GridImp::Traits::LevelIntersection>;

    template <class GridImp, IntersectionIteratorType type_>
    using ParameterSpaceIntersection =
        std::conditional_t<type_ == IntersectionIteratorType::Leaf,
                           typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLeafIntersection,
                           typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLevelIntersection>;

    template <class GridImp, IntersectionIteratorType type_>
    using ParameterSpaceIntersectionIterator =
        std::conditional_t<type_ == IntersectionIteratorType::Leaf,
                           typename GridImp::ParameterSpaceGrid::LeafGridView::IntersectionIterator,
                           typename GridImp::ParameterSpaceGrid::LevelGridView::IntersectionIterator>;

    template <class GridImp>
    using TrimInfo = typename GridImp::Trimmer::ElementTrimData;

    template <class GridImp>
    using EntityInfo = typename GridImp::Trimmer::TrimmerTraits::template Codim<0>::EntityInfo;

    template <class GridImp>
    using LocalEntity = typename GridImp::Trimmer::template Codim<0>::ParameterSpaceGridEntity;

  } // namespace IntersectionIteratorTraits

  template <class GridImp, IntersectionIteratorType type_>
  struct TrimmedIntersectionIterator
  {
    using Intersection               = IntersectionIteratorTraits::Intersection<GridImp, type_>;
    using ParameterSpaceIntersection = IntersectionIteratorTraits::ParameterSpaceIntersection<GridImp, type_>;
    using ParameterSpaceIntersectionIterator =
        IntersectionIteratorTraits::ParameterSpaceIntersectionIterator<GridImp, type_>;

    using TrimInfo    = IntersectionIteratorTraits::TrimInfo<GridImp>;
    using EntityInfo  = IntersectionIteratorTraits::EntityInfo<GridImp>;
    using LocalEntity = IntersectionIteratorTraits::LocalEntity<GridImp>;

    TrimmedIntersectionIterator() = default;

    TrimmedIntersectionIterator(const GridImp* parameterSpaceGrid, const LocalEntity& localEntity, PositionToken pos,
                                const ParameterSpaceIntersectionIterator& hostIterator)
        : patchGrid_(parameterSpaceGrid),
          localEntity_(localEntity),
          maxIterator_(localEntity.trimData().size(1)),
          iterator_(pos == PositionToken::Begin ? 0 : maxIterator_) {
      auto iterator = hostIterator;
      for (const auto i : Dune::range(4)) {
        hostIntersections_.push_back(*iterator);
        ++iterator;
      }
    }

    bool equals(const TrimmedIntersectionIterator& other) const {
      return localEntity_ == other.localEntity_ and iterator_ == other.iterator_;
    }

    void increment() {
      ++iterator_;
    }
    Intersection dereference() const {
      auto& entityInfo   = localEntity_.entityInfo();
      auto& trimData     = localEntity_.trimData();
      auto& edgeTrimInfo = trimData.edge(iterator_);
      auto edgeId   = patchGrid_->trimmer().entityContainer_.globalEdgesIdOfElementsMap_.at(entityInfo.id)[iterator_];
      auto edgeInfo = patchGrid_->trimmer().entityContainer_.idToEdgeInfoMap.at(edgeId);

      if (edgeInfo.isTrimmedHost() or not edgeInfo.isTrimmed()) {
        auto localIntersection = findLocalIntersection(edgeInfo, edgeTrimInfo, edgeInfo.isTrimmed());
        auto realIntersection  = typename Intersection::Implementation(patchGrid_, localIntersection);
        return Intersection(realIntersection);
      }
      /* trimmed Edge */
      auto localIntersection = ParameterSpaceIntersection(patchGrid_, entityInfo.id, edgeInfo, iterator_);
      auto realIntersection  = typename Intersection::Implementation(patchGrid_, localIntersection);
      return Intersection(realIntersection);
    }

  private:
    auto findLocalIntersection(const auto& edgeInfo, const auto& edgeTrimInfo,
                               bool isTrimmed) const -> ParameterSpaceIntersection {
      for (const auto i : Dune::range(4)) {
        auto trimmedIdx = Transformations::mapToTrimmer(1, i);
        if (trimmedIdx == edgeTrimInfo.idx) {
          if (isTrimmed)
            return ParameterSpaceIntersection(patchGrid_, hostIntersections_[i], edgeInfo, iterator_);
          return ParameterSpaceIntersection(patchGrid_, hostIntersections_[i], iterator_);
        }
      }
      DUNE_THROW(GridError, "Couldn't find HostIntersection");
    }

    const GridImp* patchGrid_{};
    std::vector<typename ParameterSpaceIntersectionIterator::Intersection> hostIntersections_{};
    LocalEntity localEntity_;

    unsigned int maxIterator_{};
    unsigned int iterator_{};
  };

  template <class GridImp, IntersectionIteratorType type_>
  struct HostIntersectionIterator
  {
    using Intersection               = IntersectionIteratorTraits::Intersection<GridImp, type_>;
    using ParameterSpaceIntersection = IntersectionIteratorTraits::ParameterSpaceIntersection<GridImp, type_>;
    using ParameterSpaceIntersectionIterator =
        IntersectionIteratorTraits::ParameterSpaceIntersectionIterator<GridImp, type_>;

    HostIntersectionIterator() = default;

    HostIntersectionIterator(const GridImp* parameterSpaceGrid, const ParameterSpaceIntersectionIterator& hostIterator)
        : patchGrid_(parameterSpaceGrid),
          hostIterator_(hostIterator) {}

    bool equals(const HostIntersectionIterator& other) const {
      return hostIterator_ == other.hostIterator_;
    }

    void increment() {
      ++hostIterator_;
    }
    Intersection dereference() const {
      auto parameterspaceIntersection = ParameterSpaceIntersection(patchGrid_, *hostIterator_);
      auto realIntersection           = typename Intersection::Implementation(patchGrid_, parameterspaceIntersection);
      return Intersection(realIntersection);
    }

  private:
    const GridImp* patchGrid_{};
    ParameterSpaceIntersectionIterator hostIterator_{};
  };

  template <class GridImp, IntersectionIteratorType type_>
  struct IntersectionIteratorVariant
  {
    auto visit(auto&& lambda) const {
      return std::visit(lambda, impl_);
    }
    auto visit(auto&& lambda) {
      return std::visit(lambda, impl_);
    }

    template <class Implementation>
    explicit IntersectionIteratorVariant(const Implementation& impl)
        : impl_(impl) {}
    IntersectionIteratorVariant() = default;

    void increment() {
      visit([](auto& impl) { impl.increment(); });
    }
    auto dereference() const {
      return visit([](const auto& impl) { return impl.dereference(); });
    }

    bool equals(const IntersectionIteratorVariant& other) const {
      return visit([&]<typename T>(const T& impl) { return impl.equals(std::get<T>(other.impl_)); });
    }

  private:
    std::variant<TrimmedIntersectionIterator<GridImp, type_>, HostIntersectionIterator<GridImp, type_>> impl_{};
  };
} // namespace Impl

/** @brief Iterator over all element neighbors
 * @ingroup PatchGrid
 * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
 * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
 * These neighbors are accessed via a IntersectionIterator. This allows the implement
 * non-matching meshes. The number of neighbors may be different from the number
 * of an element!
 */
template <class GridImp>
class PatchGridLeafIntersectionIterator
{
  constexpr static int dim      = GridImp::dimension;
  constexpr static int dimworld = GridImp::dimensionworld;
  typedef typename GridImp::ctype ctype;

  using ParameterSpaceIntersectionIterator = typename GridImp::ParameterSpaceGrid::LeafGridView::IntersectionIterator;
  using ParameterSpaceLeafIntersection     = typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLeafIntersection;
  using LeafIntersection                   = typename GridImp::Traits::LeafIntersection;
  using LocalEntity                        = Impl::IntersectionIteratorTraits::LocalEntity<GridImp>;

public:
  using Intersection  = Dune::Intersection<const GridImp, PatchGridLeafIntersection<GridImp>>;
  using IteratorImpl  = Impl::IntersectionIteratorVariant<GridImp, Impl::IntersectionIteratorType::Leaf>;
  using PositionToken = Impl::PositionToken;

  PatchGridLeafIntersectionIterator() = default;

  PatchGridLeafIntersectionIterator(const GridImp* parameterSpaceGrid,
                                    const ParameterSpaceIntersectionIterator& hostIterator)
      : underlying_{Impl::HostIntersectionIterator<GridImp, Impl::IntersectionIteratorType::Leaf>(parameterSpaceGrid,
                                                                                                  hostIterator)} {}

  PatchGridLeafIntersectionIterator(const GridImp* parameterSpaceGrid, const LocalEntity& localEntity,
                                    PositionToken position, const ParameterSpaceIntersectionIterator& hostIterator)
      : underlying_{Impl::TrimmedIntersectionIterator<GridImp, Impl::IntersectionIteratorType::Leaf>(
            parameterSpaceGrid, localEntity, position, hostIterator)} {}

  // equality
  bool equals(const PatchGridLeafIntersectionIterator& other) const {
    return underlying_.equals(other.underlying_);
  }

  // prefix increment
  void increment() {
    underlying_.increment();
  }

  // @brief dereferencing
  LeafIntersection dereference() const {
    return underlying_.dereference();
  }

private:
  IteratorImpl underlying_{};
};

template <class GridImp>
class PatchGridLevelIntersectionIterator
{
  constexpr static int dim      = GridImp::dimension;
  constexpr static int dimworld = GridImp::dimensionworld;
  typedef typename GridImp::ctype ctype;

  using ParameterSpaceIntersectionIterator = typename GridImp::ParameterSpaceGrid::LevelGridView::IntersectionIterator;
  using ParameterSpaceLevelIntersection    = typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLevelIntersection;
  using LevelIntersection                  = typename GridImp::Traits::LevelIntersection;
  using TrimInfo                           = typename GridImp::Trimmer::ElementTrimData;
  using EntityInfo                         = Impl::IntersectionIteratorTraits::EntityInfo<GridImp>;
  using LocalEntity                        = Impl::IntersectionIteratorTraits::LocalEntity<GridImp>;

public:
  using Intersection = Dune::Intersection<const GridImp, PatchGridLevelIntersection<GridImp>>;
  using IteratorImpl = Impl::IntersectionIteratorVariant<GridImp, Impl::IntersectionIteratorType::Level>;

  using PositionToken = Impl::PositionToken;

  PatchGridLevelIntersectionIterator() = default;

  PatchGridLevelIntersectionIterator(const GridImp* parameterSpaceGrid,
                                     const ParameterSpaceIntersectionIterator& hostIterator)
      : underlying_{Impl::HostIntersectionIterator<GridImp, Impl::IntersectionIteratorType::Level>(parameterSpaceGrid,
                                                                                                   hostIterator)} {}

  PatchGridLevelIntersectionIterator(const GridImp* parameterSpaceGrid, const LocalEntity& localEntity,
                                     PositionToken position, const ParameterSpaceIntersectionIterator& hostIterator)
      : underlying_{Impl::TrimmedIntersectionIterator<GridImp, Impl::IntersectionIteratorType::Level>(
            parameterSpaceGrid, localEntity, position, hostIterator)} {}

  // equality
  bool equals(const PatchGridLevelIntersectionIterator& other) const {
    return underlying_.equals(other.underlying_);
  }

  // prefix increment
  void increment() {
    underlying_.increment();
  }

  // @brief dereferencing
  LevelIntersection dereference() const {
    return underlying_.dereference();
  }

private:
  IteratorImpl underlying_{};
};

} // namespace Dune::IGANEW::DefaultTrim

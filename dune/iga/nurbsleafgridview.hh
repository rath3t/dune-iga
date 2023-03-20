// SPDX-FileCopyrightText: 2022 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <ranges>

#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/nurbsgridentity.hh>
#include <dune/iga/nurbsgridindexsets.hh>
#include <dune/iga/nurbsgridleafiterator.hh>
#include <dune/iga/nurbsgridtraits.hh>
#include <dune/iga/nurbspatch.hh>
#include <dune/iga/nurbstrimmedpatch.hh>

namespace Dune::IGA {

  /** \brief Collect several types associated to OneDGrid LeafGridViews */
  template <class GridImp>
  struct NurbsLeafGridViewTraits {
    using Grid        = GridImp;
    using IndexSet    = typename GridImp::Traits::LeafIndexSet;
    using GridViewImp = NURBSLeafGridView<GridImp>;

    typedef typename GridImp ::Traits ::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename GridImp ::Traits ::LeafIntersectionIterator IntersectionIterator;
    typedef typename GridImp ::Traits ::CollectiveCommunication CollectiveCommunication;
    typedef typename GridImp ::Traits ::LeafIntersection Intersection;

    template <int cd>
    struct Codim {
      typedef typename GridImp::Traits ::template Codim<cd>::template Partition<All_Partition>::LeafIterator Iterator;

      typedef typename GridImp::Traits::template Codim<cd>::Entity Entity;
      typedef typename GridImp::Traits::template Codim<cd>::Geometry Geometry;
      typedef typename GridImp::Traits::template Codim<cd>::LocalGeometry LocalGeometry;

      //      /** \brief Define types needed to iterate over entities of a given partition type */
      template <PartitionIteratorType pit>
      struct Partition {
        /** \brief iterator over a given codim and partition type */
        typedef typename GridImp::Traits::template Codim<cd>::template Partition<pit>::LeafIterator Iterator;
      };
    };
    static constexpr bool conforming = true;
  };

  template <int codim, int dim, typename GridImpl>
  class NURBSGridEntity;

  template <typename GridImpl, int griddim, std::integral auto... codim>
  std::tuple<std::vector<typename NurbsLeafGridViewTraits<GridImpl>::template Codim<codim>::Entity>...>
      gridEntityTupleGenerator(std::integer_sequence<std::common_type_t<decltype(codim)...>, codim...>);

  template <typename GridImpl>
  const auto &elements(const NURBSLeafGridView<GridImpl> &gridLeafView);
  template <typename GridImpl>
  auto &elements(NURBSLeafGridView<GridImpl> &gridLeafView);

  /** \brief NURBS leaf grid view, see Dune Book Ch. 5.1 */
  template <typename GridImpl>
  class NURBSLeafGridView {
   public:
    using NurbsGridLinearAlgebraTraits = typename GridImpl::LinearAlgebraTraits;
    using Traits                       = typename GridImpl::Traits;

    using ControlPointNetType = typename GridImpl::ControlPointNetType;

    template <int codim, int dim, typename GridImpl1>
    friend class NURBSGridEntity;

    using LeafIntersectionIterator = typename Traits::LeafIntersectionIterator;
    using IndexSet                 = typename Traits::LeafIndexSet;
    template <int cd>
    using Codim                   = typename Traits::template Codim<cd>;
    using CollectiveCommunication = typename Traits::CollectiveCommunication;
    using Intersection            = typename Traits::LeafIntersection;
    using IntersectionIterator    = LeafIntersectionIterator;

    using ctype                          = typename NurbsGridLinearAlgebraTraits::value_type;
    static constexpr auto dimension      = GridImpl::dimension;
    static constexpr auto dimensionworld = GridImpl::dimensionworld;

    template <int codim, PartitionIteratorType pitype, class GridImp1>
    friend class NURBSGridLeafIterator;

    using Grid = typename Traits::Grid;

    NURBSLeafGridView(
        const std::shared_ptr<std::vector<NURBSPatch<dimension, dimensionworld, NurbsGridLinearAlgebraTraits>>>
            &leafpatches,
        const GridImpl &grid, const std::optional<std::vector<ElementTrimFlag>> &_trimFlags = std::nullopt)
        : leafPatches_(leafpatches),
          grid_{&grid},
          entityVector_{std::make_unique<decltype(gridEntityTupleGenerator<Grid, dimension>(
              std::make_integer_sequence<int, dimension + 1>()))>()},
          trimFlags(_trimFlags) {
      // Make Entities
      createEntities();

      indexSet_ = std::make_unique<NURBSGridLeafIndexSet<GridImpl>>(*this);
    }

    NURBSLeafGridView(const NURBSLeafGridView &other) {
      if (this == &other) return;
      leafPatches_ = other.leafPatches_;
      grid_        = other.grid_;

      trimFlags = other.trimFlags;
      indexSet_ = std::make_unique<NURBSGridLeafIndexSet<GridImpl>>(*this);
      updateGridViewForEntities();
    }

    NURBSLeafGridView &operator=(const NURBSLeafGridView &other) {
      if (this == &other) return *this;
      leafPatches_ = other.leafPatches_;
      grid_        = other.grid_;

      trimFlags = other.trimFlags;
      indexSet_.reset();
      indexSet_ = std::make_unique<NURBSGridLeafIndexSet<GridImpl>>(*this);
      updateGridViewForEntities();
      return *this;
    }

    /** \brief obtain collective communication object, currently this return a "NoComm" object */
    const auto &comm() const { return grid().comm(); }

    template <class Entity>
    bool contains(const Entity &e) const {
      return (e.impl().NURBSGridView_) == this;
    }

    const auto &grid() const { return *grid_; }

    /** \brief obtain how many entities of a give geometry type do live in this grid view */
    [[nodiscard]] int size(const GeometryType &type) const {
      if (type == Dune::GeometryTypes::vertex || type == Dune::GeometryTypes::cube(1)
          || type == Dune::GeometryTypes::cube(2) || type == Dune::GeometryTypes::cube(3))
        return this->size(dimension - type.dim());
      else
        return 0;
    }

    const auto &getPatchData(int patchID = 0) const { return *(leafPatches_->at(patchID).getPatchData()); }
    const auto &getPatch(int patchID = 0) const { return leafPatches_->at(patchID); }

    auto getPreBasis() {
      assert(leafPatches_->size() == 1 && "The basis is only defined for single patch gridview");
      return Dune::Functions::BasisFactory::nurbs<dimension>(this->getPatchData(0));
    }

    template <int cd, Dune::PartitionIteratorType ptype = Dune::All_Partition>
    using LeafIteratorImpl = NURBSGridLeafIterator<cd, ptype, GridImpl>;

    template <int cd>
    typename Codim<cd>::LeafIterator begin() const {
      return LeafIteratorImpl<cd>(std::get<cd>(*entityVector_.get()).begin());
    }

    template <int cd>
    typename Codim<cd>::LeafIterator end() const {
      return LeafIteratorImpl<cd>(std::get<cd>(*entityVector_.get()).end());
    }

    LeafIntersectionIterator ibegin(const typename Codim<0>::Entity &entity) const {
      return entity.impl().ibegin(level_);
    }

    LeafIntersectionIterator iend(const typename Codim<0>::Entity &entity) const { return entity.impl().iend(level_); }

    template <int cd, PartitionIteratorType piType>
    typename Codim<cd>::template Partition<piType>::LeafIterator begin() const {
      if (piType != Ghost_Partition)
        return LeafIteratorImpl<cd, piType>(std::get<cd>(*entityVector_.get()).begin());
      else
        return LeafIteratorImpl<cd, piType>(std::get<cd>(*entityVector_.get()).end());
    }

    template <int cd, PartitionIteratorType piType>
    typename Codim<cd>::template Partition<piType>::LeafIterator end() const {
      return LeafIteratorImpl<cd, piType>(std::get<cd>(*entityVector_.get()).end());
    }

    const IndexSet &indexSet() const { return *indexSet_; }
    [[nodiscard]] int overlapSize(int codim) const { return 0; }
    [[nodiscard]] int ghostSize(int codim) const { return 0; }
    auto size(int codim) const {
      assert(codim <= 3 && codim >= 0);
      if (codim == 0)
        return std::get<0>(*entityVector_.get()).size();
      else if (codim == 1)
        return std::get<1>(*entityVector_.get()).size();
      if constexpr (dimension > 1)
        if (codim == 2) return std::get<2>(*entityVector_.get()).size();
      if constexpr (dimension > 2)
        if (codim == 3) return std::get<3>(*entityVector_.get()).size();
      __builtin_unreachable();
    }

   private:
    void createEntities() {
      for (int currentPatchId = 0; auto &&patch : *leafPatches_.get()) {
        Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<dimension + 1>()), [&](const auto i) {
          std::get<i>(*entityVector_.get()).reserve(patch.size(i));

          // Make Elements (codim = 0)
          if (i == 0) {
            for (unsigned int j = 0; j < patch.size(i); ++j) {
              if (trimFlags.has_value()) {
                std::get<0>(*entityVector_.get())
                    .emplace_back(
                        NURBSGridEntity<0, dimension, GridImpl>(*this, j, currentPatchId, trimFlags.value()[j]));
              } else {
                std::get<0>(*entityVector_.get())
                    .emplace_back(NURBSGridEntity<0, dimension, GridImpl>(*this, j, currentPatchId));
              }
            }
          } else {
            // Codim != 0
            for (unsigned int j = 0; j < patch.size(i); ++j) {
              std::get<i>(*entityVector_.get())
                  .emplace_back(NURBSGridEntity<i, dimension, GridImpl>(*this, j, currentPatchId));
            }
          }
        });
        ++currentPatchId;
      }
    }

    void updateGridViewForEntities() {
      entityVector_ = std::make_unique<decltype(gridEntityTupleGenerator<Grid, dimension>(
          std::make_integer_sequence<int, dimension + 1>()))>();

      createEntities();
    }
    template <int codim>
    typename Codim<codim>::Entity &getEntity(unsigned int directIndex) const {
      if constexpr (codim == 0)  // elements
        return (std::get<0>(*entityVector_.get()).at(directIndex));
      else if constexpr (codim == dimension)  // vertices
        return (std::get<dimension>(*entityVector_.get()).at(directIndex));
      else if constexpr (dimension - codim == 1)  // edges
        return (std::get<dimension - 1>(*entityVector_.get()).at(directIndex));
      else if constexpr (dimension - codim == 2)  // surface
        return (std::get<dimension - 2>(*entityVector_.get()).at(directIndex));
      else
        throw std::logic_error("Your requested entity type does not exist.");
    }

    friend GridImpl;
    friend const auto &elements<GridImpl>(const NURBSLeafGridView<GridImpl> &gridLeafView);
    friend auto &elements<GridImpl>(NURBSLeafGridView<GridImpl> &gridLeafView);
    template <typename GridImp1>
    friend class NURBSintersection;
    std::shared_ptr<std::vector<NURBSPatch<dimension, dimensionworld, NurbsGridLinearAlgebraTraits>>> leafPatches_;
    const GridImpl *grid_;
    using EntityVectorType
        = decltype(gridEntityTupleGenerator<GridImpl, dimension>(std::make_integer_sequence<int, dimension + 1>()));
    std::unique_ptr<EntityVectorType> entityVector_{};
    std::unique_ptr<NURBSGridLeafIndexSet<GridImpl>> indexSet_;
    int level_{};

    std::optional<std::vector<ElementTrimFlag>> trimFlags;
  };

  //  template <typename GridImpl>
  //  const auto &elements(const NURBSLeafGridView<GridImpl> &gridLeafView) {
  //    return std::get<0>(*gridLeafView.entityVector_.get());
  //  }
  //
  //  template <typename GridImpl>
  //  auto &elements(NURBSLeafGridView<GridImpl> &gridLeafView) {
  //    return std::get<0>(*gridLeafView.entityVector_.get());
  //  }

  template <typename GridImpl1, typename ElementEntity>
  auto &intersections(const typename GridImpl1::GridView &gridLeafView, const ElementEntity &e) {
    return *e.intersections_.get();
  }
}  // namespace Dune::IGA

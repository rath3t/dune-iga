// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/nurbsgridentity.hh"
#include "dune/iga/nurbsgridindexsets.hh"
#include "dune/iga/nurbsgridleafiterator.hh"
#include "dune/iga/nurbsgridtraits.hh"
#include "dune/iga/nurbspatch.hh"
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune::IGA {

  /** \brief Collect several types associated to OneDGrid LeafGridViews */
  template <class GridImp>
  struct NurbsLeafGridViewTraits {
    using Grid        = GridImp;
    using IndexSet    = typename GridImp::Traits::LeafIndexSet;
    using GridViewImp = NURBSLeafGridView<const GridImp>;

    typedef typename GridImp ::Traits ::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename GridImp ::Traits ::LeafIntersectionIterator IntersectionIterator;
    typedef typename GridImp ::Traits ::CollectiveCommunication CollectiveCommunication;
    typedef typename GridImp ::Traits ::CollectiveCommunication Communication;
    typedef typename GridImp ::Traits ::LeafIntersection Intersection;

    template <int cd>
    struct Codim {
      typedef typename GridImp::Traits ::template Codim<cd>::template Partition<All_Partition>::LeafIterator Iterator;

      typedef typename GridImp::Traits::template Codim<cd>::Entity Entity;
      typedef typename GridImp::Traits::template Codim<cd>::Geometry Geometry;
      typedef typename GridImp::Traits::template Codim<cd>::LocalGeometry LocalGeometry;

      ///** \brief Define types needed to iterate over entities of a given partition type */
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

  template <typename GridImpl>
  const auto &elements(const NURBSLeafGridView<const GridImpl> &gridLeafView);
  template <typename GridImpl>
  auto &elements(NURBSLeafGridView<const GridImpl> &gridLeafView);

  /** \brief NURBS leaf grid view, see Dune Book Ch. 5.1 */
  template <typename GridImpl>
  class NURBSLeafGridView {
   public:
    using Traits = typename GridImpl::Traits;

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

    using ctype                          = typename GridImpl::ctype;
    static constexpr auto dimension      = GridImpl::dimension;
    static constexpr auto dimensionworld = GridImpl::dimensionworld;

    template <int codim, PartitionIteratorType pitype, class GridImp1>
    friend class NURBSGridLeafIterator;

    using Grid = typename Traits::Grid;

    [[nodiscard]] bool isConforming() const { return true; }

    NURBSLeafGridView(GridImpl &grid, int p_level) : grid_{&grid}, level_{p_level} {
      //      indexSet_ = std::make_unique<NURBSGridLeafIndexSet<GridImpl>>(*this);
    }

    NURBSLeafGridView(const NURBSLeafGridView &other) {
      if (this == &other) return;
      grid_  = other.grid_;
      level_ = other.level_;
    }
    //
    NURBSLeafGridView &operator=(const NURBSLeafGridView &other) {
      if (this == &other) return *this;
      grid_  = other.grid_;
      level_ = other.level_;
      return *this;
    }

    /** \brief obtain collective communication object, currently this return a "NoComm" object */
    const auto &comm() const { return grid().comm(); }

    template <class Entity>
    bool contains(const Entity &e) const {
      return &(e.impl().NURBSGridView_->impl()) == this;
    }

    const auto &grid() const { return *grid_; }

    /** \brief obtain how many entities of a give geometry type do live in this grid view */
    [[nodiscard]] int size(const GeometryType &type) const {
      if constexpr (dimension != 2) {
        if (type == Dune::GeometryTypes::vertex || type == Dune::GeometryTypes::cube(1)
            || type == Dune::GeometryTypes::cube(2) || type == Dune::GeometryTypes::cube(3))
          return this->size(dimension - type.dim());
        else
          return 0;
      } else {
        if (type == Dune::GeometryTypes::vertex || type == Dune::GeometryTypes::cube(1))
          return this->size(dimension - type.dim());
        else if (type == Dune::GeometryTypes::none(2))
          return getPatch().n_trimmedElement;
        else if (type == Dune::GeometryTypes::cube(2))
          return getPatch().n_fullElement;
        else
          return 0;
      }
    }

    const auto &getPatchData(int patchID = 0) const { return *(grid_->leafPatches_->at(patchID).getPatchData()); }
    const auto &lowerOrderPatchData(int patchID = 0) const { return grid_->lowerOrderPatchData(patchID); }
    const auto &getPatch(int patchID = 0) const { return grid_->leafPatches_->at(patchID); }

    template <int cd, Dune::PartitionIteratorType ptype = Dune::All_Partition>
    using LeafIteratorImpl = NURBSGridLeafIterator<cd, ptype, const GridImpl>;

    template <int cd>
    typename Codim<cd>::LeafIterator begin() const {
      return LeafIteratorImpl<cd>(std::get<cd>(*grid_->entityVector.get()).begin());
    }

    template <int cd>
    typename Codim<cd>::LeafIterator end() const {
      return LeafIteratorImpl<cd>(std::get<cd>(*grid_->entityVector.get()).end());
    }

    LeafIntersectionIterator ibegin(const typename Codim<0>::Entity &entity) const {
      return entity.impl().ibegin(level_);
    }

    LeafIntersectionIterator iend(const typename Codim<0>::Entity &entity) const { return entity.impl().iend(level_); }

    template <class DataHandleImp, class DataType>
    void communicate([[maybe_unused]] CommDataHandleIF<DataHandleImp, DataType> &data,
                     [[maybe_unused]] InterfaceType iftype, [[maybe_unused]] CommunicationDirection dir) const {}

    template <int cd, PartitionIteratorType piType>
    typename Codim<cd>::template Partition<piType>::LeafIterator begin() const {
      if (piType != Ghost_Partition)
        return LeafIteratorImpl<cd, piType>(std::get<cd>(*grid_->entityVector.get()).begin());
      else
        return LeafIteratorImpl<cd, piType>(std::get<cd>(*grid_->entityVector.get()).end());
    }

    template <int cd, PartitionIteratorType piType>
    typename Codim<cd>::template Partition<piType>::LeafIterator end() const {
      return LeafIteratorImpl<cd, piType>(std::get<cd>(*grid_->entityVector.get()).end());
    }

    const IndexSet &indexSet() const { return *grid_->indexSet_; }
    [[nodiscard]] int overlapSize(int codim) const { return 0; }
    [[nodiscard]] int ghostSize(int codim) const { return 0; }
    auto size(int codim) const {
      assert(codim <= 3 && codim >= 0);
      if (codim == 0)
        return std::get<0>(*grid_->entityVector.get()).size();
      else if (codim == 1)
        return std::get<1>(*grid_->entityVector.get()).size();
      if constexpr (dimension > 1)
        if (codim == 2) return std::get<2>(*grid_->entityVector.get()).size();
      if constexpr (dimension > 2)
        if (codim == 3) return std::get<3>(*grid_->entityVector.get()).size();
      __builtin_unreachable();
    }

   private:
    template <int codim>
    typename Codim<codim>::Entity &getEntity(unsigned int directIndex) const {
      if constexpr (codim == 0)  // elements
        return (std::get<0>(*grid_->entityVector.get()).at(directIndex));
      else if constexpr (codim == dimension)  // vertices
        return (std::get<dimension>(*grid_->entityVector.get()).at(directIndex));
      else if constexpr (dimension - codim == 1)  // edges
        return (std::get<dimension - 1>(*grid_->entityVector.get()).at(directIndex));
      else if constexpr (dimension - codim == 2)  // surface
        return (std::get<dimension - 2>(*grid_->entityVector.get()).at(directIndex));
      else
        throw std::logic_error("Your requested entity type does not exist.");
    }

    friend GridImpl;

    template <typename GridImp1>
    friend class NURBSintersection;
    //    std::shared_ptr<std::vector<NURBSPatch<dimension, dimensionworld, ctype>>> leafPatches_;
    GridImpl *grid_;
    int level_{};
  };

  //  template <typename GridImpl>
  //  const auto &elements(const NURBSLeafGridView<GridImpl> &gridLeafView) {
  //    return std::get<0>(*gridLeafView.entityVector.get());
  //  }
  //
  //  template <typename GridImpl>
  //  auto &elements(NURBSLeafGridView<GridImpl> &gridLeafView) {
  //    return std::get<0>(*gridLeafView.entityVector.get());
  //  }

  template <typename GridImpl1, typename ElementEntity>
  auto &intersections(const typename GridImpl1::GridView &gridLeafView, const ElementEntity &e) {
    return *e.intersections_.get();
  }
}  // namespace Dune::IGA

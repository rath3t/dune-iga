// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <span>

#include <dune/common/iteratorfacades.hh>
//#include "dune/iga/nurbsleafgridview.hh"
#include "dune/iga/nurbspatch.hh"

namespace Dune::IGA {

  template <int codim, PartitionIteratorType pitype, class GridImp>
  class NURBSGridLeafIterator {
    constexpr static int dim = GridImp::dimension;

   public:
    using Entity                     = typename GridImp::Traits::template Codim<codim>::Entity;
    constexpr static int codimension = codim;

    NURBSGridLeafIterator(typename std::vector<Entity>::const_iterator virtualEntity) : virtualEntity_(virtualEntity) {}

    NURBSGridLeafIterator() = default;
    //! prefix increment
    void increment() { ++virtualEntity_; }

    //! dereferencing
    const Entity& dereference() const { return *virtualEntity_; }

    //! equality
    bool equals(const NURBSGridLeafIterator& other) const { return virtualEntity_ == other.virtualEntity_; }

   private:
    //    /** \brief This increment makes the iterator wander over all entities on all levels */
    //    void globalIncrement() {
    //
    //      // Backup current level because it may not be accessible anymore after
    //      // setting the pointer to the next entity.
    //      const int oldLevel = this->virtualEntity_.level();
    //
    //      // Increment on this level
    //      this->virtualEntity_.impl().setToTarget(this->virtualEntity_.impl().target_->succ_);
    //
    //      // If beyond the end of this level set to first of next level
    //      if (!this->virtualEntity_.impl().target_ && oldLevel < grid_->maxLevel()) {
    //
    //        this->virtualEntity_.impl().setToTarget(const_cast<OneDEntityImp<dim-codim>*>(std::get<1-codim>(grid_->entityImps_[oldLevel+1]).begin()));
    //
    //      }
    //
    //    }

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    //    const GridImp* grid_;

    //! The entity that the iterator is pointing to
    std::vector<Entity>::const_iterator virtualEntity_{nullptr};
  };

  template <class GridImp>
  struct NurbsHierarchicIterator {
    explicit NurbsHierarchicIterator(const typename GridImp::Traits::template Codim<0>::Entity& ent)
        : nurbsEntity{&ent} {}
    using Entity                                           = typename GridImp::Traits::template Codim<0>::Entity;
    auto operator<=>(const NurbsHierarchicIterator&) const = default;
    auto operator*() { return *nurbsEntity; }

    const Entity& dereference() const { return *nurbsEntity; }
    bool equals(const NurbsHierarchicIterator& r) const { return *this == r; }
    void increment() {}
    //    auto operator->() { return nurbsEntity; }
    //    void operator++() {}
    const typename GridImp::Traits::template Codim<0>::Entity* nurbsEntity;
  };

  template <class GridImp>
  class NURBSGridInterSectionIterator : public std::vector<typename GridImp::Traits::LeafIntersection>::const_iterator {
    using Intersection = typename GridImp::Traits::LeafIntersection;
    using Base         = typename std::vector<typename GridImp::Traits::LeafIntersection>::const_iterator;

   public:
    //! copy constructor
    NURBSGridInterSectionIterator(const NURBSGridInterSectionIterator& other) = default;

    const Intersection& dereference() const { return **this; }
    auto operator<=>(const NURBSGridInterSectionIterator&) const = default;
    void increment() { ++(*this); }
    bool equals(const NURBSGridInterSectionIterator& r) const { return *this == r; }
    NURBSGridInterSectionIterator() = default;
    using Reference                 = Intersection;
    explicit NURBSGridInterSectionIterator(typename std::vector<Intersection>::const_iterator spanIter)
        : std::vector<Intersection>::const_iterator(spanIter) {}
  };
}  // namespace Dune::IGA

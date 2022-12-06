// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#pragma once

#include <span>

//#include <dune/common/iteratorfacades.hh>
#include <dune/iga/nurbsleafgridview.hh>
#include <dune/iga/nurbspatch.hh>

namespace Dune::IGA {
  /** \brief NURBS gird leaf iterator */
  template <int codim, PartitionIteratorType pitype, class GridImp>
  class NURBSGridLeafIterator : public std::vector<typename GridImp::Traits::template Codim<codim>::Entity>::const_iterator {
    using Base = typename std::vector<typename GridImp::Traits::template Codim<codim>::Entity>::const_iterator;

  public:
    NURBSGridLeafIterator() = default;
    using Reference         = typename GridImp::Traits::template Codim<codim>::Entity;
    using Entity            = typename GridImp::Traits::template Codim<codim>::Entity;

    void increment() {++(*this); }
    //    using Base::operator*;
    //    using Base::operator->;

    const Entity& dereference() const { return **this; }

    bool equals(const NURBSGridLeafIterator& r) const { return *this == r; }
        explicit NURBSGridLeafIterator(typename std::vector<Entity>::const_iterator spanIter) : std::vector<Entity>::const_iterator(spanIter) {}
  };


    /** \brief Iterator over child elements
     * This is a default implementation since the functionality is not given for the iga grid
     * */
    template<class GridImp>
    struct NurbsHierarchicIterator
    {
      explicit NurbsHierarchicIterator(const  typename GridImp::Traits::template Codim<0>::Entity& ent) : nurbsEntity{&ent}{}
      auto operator<=>(const NurbsHierarchicIterator&) const = default;
      auto operator*(){ return *nurbsEntity;}
      auto operator->(){ return nurbsEntity;}
      void operator++(){}
      const typename GridImp::Traits::template Codim<0>::Entity* nurbsEntity;
    };
    /** \brief Iterator over intersections between elements  */
    template<class GridImp>
    class NURBSGridInterSectionIterator : public std::vector<typename GridImp::Traits::LeafIntersection>::const_iterator
    {
    public:
      NURBSGridInterSectionIterator()=default;
      using Intersection = typename GridImp::Traits::LeafIntersection;
      using Reference =  Intersection;
      explicit NURBSGridInterSectionIterator(typename std::vector<Intersection>::const_iterator spanIter)
          :  std::vector<Intersection>::const_iterator(spanIter)
      {}
    };
  }



  template <class GridImp>
  struct NurbsHierarchicIterator {
    explicit NurbsHierarchicIterator(const typename GridImp::Traits::template Codim<0>::Entity& ent) : nurbsEntity{&ent} {}
    using Entity                                           = typename GridImp::Traits::template Codim<0>::Entity;
    auto operator<=>(const NurbsHierarchicIterator&) const = default;
    auto operator*() { return *nurbsEntity; }

    const Entity& dereference() const { return *nurbsEntity; }
    bool equals(const NurbsHierarchicIterator& r) const { return *this==r; }
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
    NURBSGridInterSectionIterator (const NURBSGridInterSectionIterator& other) = default;

    const Intersection& dereference() const { return **this; }
    auto operator<=>(const NURBSGridInterSectionIterator&) const = default;
    void increment()  { ++(*this); }
    bool equals(const NURBSGridInterSectionIterator& r) const { return *this==r; }
    NURBSGridInterSectionIterator() = default;
    using Reference                 = Intersection;
    explicit NURBSGridInterSectionIterator(typename std::vector<Intersection>::const_iterator spanIter)
        : std::vector<Intersection>::const_iterator(spanIter) {}
  };
}  // namespace Dune::IGA

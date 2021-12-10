// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#pragma once

#include <span>

#include <dune/common/iteratorfacades.hh>
#include <dune/iga/nurbsleafgridview.hh>
#include <dune/iga/nurbspatch.hh>

 namespace Dune::IGA
  {
    /** \brief NURBS gird leaf iterator */
    template<int codim, PartitionIteratorType pitype, class GridImp>
    class NURBSGridLeafIterator : public std::vector<typename GridImp::Traits::template Codim<codim>::Entity>::const_iterator
    {
    public:
      NURBSGridLeafIterator()=default;
      using Reference =  typename GridImp::Traits::template Codim<codim>::Entity;
      using Entity =  typename GridImp::Traits::template Codim<codim>::Entity;

      auto dereference()
      {return *this;}

      explicit NURBSGridLeafIterator(typename std::vector<Entity>::const_iterator spanIter)
          :  std::vector<Entity>::const_iterator(spanIter)
      {}
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


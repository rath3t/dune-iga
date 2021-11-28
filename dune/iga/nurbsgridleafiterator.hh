// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#pragma once

#include <span>

#include <dune/common/iteratorfacades.hh>
#include <dune/iga/nurbsleafgridview.hh>
#include <dune/iga/nurbspatch.hh>
/** \file
 * \brief The NURBSGridIterator class
 */
 namespace Dune::IGA
  {
    /** \brief NURBS gird leaf iterator */
    template<typename NURBSEntity>
  class NURBSGridLeafIterator : public std::vector<NURBSEntity>::const_iterator
    {
    public:
      NURBSGridLeafIterator()=default;
      using Reference =  NURBSEntity;
      using Entity =  NURBSEntity;

      explicit NURBSGridLeafIterator(typename std::vector<NURBSEntity>::const_iterator spanIter)
          :  std::vector<NURBSEntity>::const_iterator(spanIter)
      {}
    };

    template<typename NURBSEntity>
    struct NurbsHierarchicIterator
    {
      explicit NurbsHierarchicIterator(const  NURBSEntity& ent) : nurbsEntity{&ent}{}
      auto operator<=>(const NurbsHierarchicIterator&) const = default;
      auto operator*(){ return *nurbsEntity;}
      auto operator->(){ return nurbsEntity;}
      void operator++(){}
      const NURBSEntity* nurbsEntity;
    };

    template<typename NURBSIntersection>
    class NURBSGridInterSectionIterator : public std::vector<NURBSIntersection>::const_iterator
    {
    public:
      NURBSGridInterSectionIterator()=default;
      using Reference =  NURBSIntersection;
      using Intersection = NURBSIntersection;
      explicit NURBSGridInterSectionIterator(typename std::vector<NURBSIntersection>::const_iterator spanIter)
          :  std::vector<NURBSIntersection>::const_iterator(spanIter)
      {}
    };
  }


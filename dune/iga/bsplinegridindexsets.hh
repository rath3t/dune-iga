// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEGRID_INDEXSETS_HH
#define DUNE_IGA_BSPLINEGRID_INDEXSETS_HH

#include <dune/grid/common/indexidset.hh>
#include <dune/iga/bsplinegridentity.hh>

namespace Dune::IGA
  {
    template<class GridImp>
    class BSplineGridLeafIndexSet
    {
    public:
      typedef int IndexType;
      //! constructor
      explicit BSplineGridLeafIndexSet (const GridImp& g) : grid_(g)
      {}

      //! get index of an entity, need to change the type to a property from GridImp
      template<int codim>
      int index (const BSplineGridEntity<codim, GridImp> & e) const
      {
        return e.getIndex();
      }


    private:
      const GridImp& grid_;
    };
  }

#endif

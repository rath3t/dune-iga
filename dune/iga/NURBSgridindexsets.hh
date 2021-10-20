// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_NURBSGRID_INDEXSETS_HH
#define DUNE_IGA_NURBSGRID_INDEXSETS_HH

#include <dune/grid/common/indexidset.hh>
#include <dune/iga/NURBSgridentity.hh>

namespace Dune::IGA
  {
    template<class GridImp>
    class NURBSGridLeafIndexSet
    {
    public:
      //! constructor
      explicit NURBSGridLeafIndexSet (const GridImp& g) : grid_(g)
      {}

      //! get index of an entity, need to change the type to a property from GridImp
      template<int codim>
      int index (const NURBSGridEntity<codim, GridImp> & e) const
      {
        return e.getIndex();
      }


    private:
      const GridImp& grid_;
    };
  }

#endif

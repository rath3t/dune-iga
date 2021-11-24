// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#pragma once

#include <span>

#include <dune/common/iteratorfacades.hh>
#include <dune/iga/NURBSleafgridview.hh>
#include <dune/iga/NURBSpatch.hh>
/** \file
 * \brief The NURBSGridIterator class
 */

namespace Dune::IGA {
  /** \brief NURBS gird leaf iterator */
  template <typename NURBSEntity>
  class NURBSGridLeafIterator : public std::vector<NURBSEntity>::iterator {
  public:
    NURBSGridLeafIterator() = default;
    using Reference         = NURBSEntity;
    explicit NURBSGridLeafIterator(const typename std::vector<NURBSEntity>::iterator& spanIter)
        : std::vector<NURBSEntity>::iterator(spanIter) {}
  };
}  // namespace Dune::IGA

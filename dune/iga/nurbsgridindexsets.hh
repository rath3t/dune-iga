// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <map>
#include <span>

#include "dune/iga/nurbsgridentity.hh"
#include <dune/grid/common/indexidset.hh>

namespace Dune::IGA {

  template <class GridImpl>
  class NURBSGridLeafIndexSet
      : public IndexSet<GridImpl, NURBSGridLeafIndexSet<GridImpl>, int, std::vector<GeometryType>> {
   public:
    using Types                    = std::vector<GeometryType>;
    using IndexType                = unsigned int;
    using GridView                 = typename GridImpl::Traits::LeafGridView;
    static constexpr auto griddim  = GridImpl::dimension;
    static constexpr auto dimworld = GridImpl::dimensionworld;
    //! Forbid the copy constructor
    NURBSGridLeafIndexSet(const NURBSGridLeafIndexSet&) = delete;
    //! Forbid the assignment operator
    NURBSGridLeafIndexSet& operator=(const NURBSGridLeafIndexSet&) = delete;

    explicit NURBSGridLeafIndexSet(NURBSLeafGridView<GridImpl> const& g) : gridView_(&g) {}

    template <class Entity>
    bool contains(const Entity& e) const {
      return gridView_->contains(e);
    }

    /** \brief get index of an entity */
    template <int codim>
    IndexType index(const typename GridImpl::Traits::template Codim<codim>::Entity& e) const {
      return e.impl().getIndex();
    }

    template <class Entity>
    IndexType index(const Entity& e) const {
      return e.impl().getIndex();
    }

    template <int codimElement>
    IndexType subIndex(const typename GridImpl::Traits::template Codim<codimElement>::Entity& e, int i,
                       unsigned int codim) const {
      if (codimElement == 0 && NURBSGridEntity<codimElement, griddim, GridImpl>::mydimension == codim)
        return gridView_->getPatch(0).getGlobalVertexIndexFromElementIndex(e.impl().getIndex(), i);
      else if (i == 0 && codim == 0 && codimElement == 0)
        return this->index(e);
      else if ((codim == 1 && codimElement == 0 && griddim == 2) || (codim == 2 && codimElement == 0 && griddim == 3))
        return gridView_->getPatch(0).getGlobalEdgeIndexFromElementIndex(e.impl().getIndex(), i);
      else if (codim == 1 && griddim == 3)  // surface case
        return gridView_->getPatch(0).getGlobalSurfaceIndexFromElementIndex(e.impl().getIndex(), i);
      else
        throw std::logic_error("subIndex only defined from element to vertices, edges and surfaces");
    }

    // At construction time the types are not yet known, so we obtain them here
    auto& types(int codim) const {
      for (int id = 0; auto& entTypes : types_)
        entTypes = gridView_->getPatch().typesInCodim(griddim - (id++));
      return types_.at(codim);
    }
    auto size(int codim) const { return gridView_->size(codim); }
    auto size(const GeometryType& gt) const { return gridView_->size(griddim - gt.dim()); }

   private:
    mutable std::array<Types, griddim + 1> types_;
    NURBSLeafGridView<GridImpl> const* gridView_;
  };
}  // namespace Dune::IGA

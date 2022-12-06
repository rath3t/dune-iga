
#pragma once

#include <map>
#include <span>

#include <dune/grid/common/indexidset.hh>
#include <dune/iga/nurbsgridentity.hh>

namespace Dune::IGA {

  template <class GridImpl>
  class NURBSGridLeafIndexSet : public IndexSet< GridImpl, NURBSGridLeafIndexSet< GridImpl >, int ,std::array<GeometryType, 1>> {
  public:
    using Types                    = std::array<GeometryType, 1>;
    using IndexType                = unsigned int;
    using GridView                 = typename GridImpl::Traits::LeafGridView;
    static constexpr auto griddim  = GridImpl::dimension;
    static constexpr auto dimworld = GridImpl::dimensionworld;

    explicit NURBSGridLeafIndexSet(NURBSLeafGridView<GridImpl> const& g) : gridView_(&g) {
      for (int id = 0; auto& entTypes : types_)
        entTypes = {GeometryTypes::cube(griddim - (id++))};
    }

    template <class Entity>
    bool contains(const Entity& e) const {
      return gridView_->contains(e);
    }

    /** \brief get index of an entity */
    template <int codim>
    IndexType index(const typename GridImpl::Traits::template Codim<codim>::Entity& e) const {
      return e.impl().getIndex();
    }

    template<class Entity>
    IndexType index(const  Entity& e) const {
      return e.impl().getIndex();
    }

    template <int codimElement>
    IndexType subIndex(const typename GridImpl::Traits::template Codim<codimElement>::Entity& e, int i, unsigned int codim) const {
      if (codimElement == 0 && NURBSGridEntity<codimElement,griddim, GridImpl>::mydimension == codim)
        return gridView_->getPatch(0).getGlobalVertexIndexFromElementIndex(e.impl().getIndex(), i);
      else if (i == 0 && codim == 0 && codimElement == 0)
        return this->index(e);
      else if ((codim == 1 && codimElement == 0 && griddim == 2 )|| (codim == 2 && codimElement == 0 && griddim == 3))
        return gridView_->getPatch(0).getGlobalEdgeIndexFromElementIndex(e.impl().getIndex(), i);
      else if (codim==1 && griddim == 3) // surface case
        return gridView_->getPatch(0).getGlobalSurfaceIndexFromElementIndex(e.impl().getIndex(), i);
      else
        throw std::logic_error("subIndex only defined from element to vertices, edges and surfaces");
    }

    auto& types(int codim) const { return types_[codim]; }
    auto size(int codim) const { return gridView_->size(codim); }
    auto size(const GeometryType& gt) const { return gridView_->size(griddim - gt.dim()); }

  private:
    std::array<std::array<GeometryType, 1>, griddim + 1> types_;
    NURBSLeafGridView<GridImpl> const* gridView_;
  };
}  // namespace Dune::IGA

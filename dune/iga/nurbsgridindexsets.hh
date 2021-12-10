
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
    using IndexType                =  int;
    using GridView = typename GridImpl::Traits::LeafGridView;
    static constexpr auto griddim  = GridImpl::dimension;
    static constexpr auto dimworld = GridImpl::dimensionworld;

    //! constructor
    explicit NURBSGridLeafIndexSet(const GridView& g) : gridView_(&g) {
      for (int id = 0; auto& entTypes : types_) {
        entTypes = {GeometryTypes::cube(griddim - id)};
        ++id;
      }
    }

    template <class Entity>
    bool contains(const Entity& e) const {
      return gridView_->contains(e);
    }

    //! get index of an entity, need to change the type to a property from GridImp
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
      if (codimElement == 0 && NURBSGridEntity<codimElement,griddim, GridImpl>::mydim == codim)
        return gridView_->impl().NURBSpatch_->getGlobalVertexIndexFromElementIndex(e.impl().getIndex(), i);
      else if (i == 0 && codim == 0 && codimElement == 0)
        return this->index(e);
      else if ((codim == 1 && codimElement == 0 && griddim == 2 )|| (codim == 2 && codimElement == 0 && griddim == 3))
        return gridView_->impl().NURBSpatch_->getGlobalEdgeIndexFromElementIndex(e.impl().getIndex(), i);
      else if (codim==1 && griddim == 3) // surface case
        return gridView_->impl().NURBSpatch_->getGlobalSurfaceIndexFromElementIndex(e.impl().getIndex(), i);
      else
        throw std::logic_error("subIndex only defined from element to vertices, edges and surfaces");
    }

    auto& types(int codim) const { return types_[codim]; }

    auto size(int codim) const { return gridView_->size(codim); }

    auto size(const GeometryType& gt) const {
      const int codim = griddim - gt.dim();
      return gridView_->size(codim);
    }

  private:
    std::array<std::array<GeometryType, 1>, griddim + 1> types_;
    const GridView* gridView_;
  };
}  // namespace Dune::IGA

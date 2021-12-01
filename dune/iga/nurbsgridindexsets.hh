
#pragma once

#include <map>
#include <span>

#include <dune/grid/common/indexidset.hh>
#include <dune/iga/nurbsgridentity.hh>

namespace Dune::IGA {
  template <class GridImpl>
  class NURBSGridLeafIndexSet {
  public:
    using Types                    = std::array<GeometryType, 1>;
    using IndexType                = unsigned int;
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
    int index(const NURBSGridEntity<codim,griddim, GridImpl>& e) const {
      return e.getIndex();
    }

    template <int codimElement>
    int subIndex(const NURBSGridEntity<codimElement,griddim, GridImpl>& e, int i, unsigned int codim) const {
      if (codimElement == 0 && NURBSGridEntity<codimElement,griddim, GridImpl>::mydim == codim)
        return gridView_->NURBSpatch_->getGlobalVertexIndexFromElementIndex(e.getIndex(), i);
      else if (i == 0 && codim == 0 && codimElement == 0)
        return this->index(e);
      else if ((codim == 1 && codimElement == 0 && griddim == 2 )|| (codim == 2 && codimElement == 0 && griddim == 3))
        return gridView_->NURBSpatch_->getGlobalEdgeIndexFromElementIndex(e.getIndex(), i);
      else if (codim==1 && griddim == 3) // surface case
        return gridView_->NURBSpatch_->getGlobalSurfaceIndexFromElementIndex(e.getIndex(), i);
      else
        throw std::logic_error("subIndex only defined from element to vertices and edges");
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

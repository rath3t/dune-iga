//
// Created by lex on 16.11.21.
//

#pragma once

namespace Dune::IGA {
  template <typename GridImpl>
  class IgaIdSet {
  public:
    using IdType = int;
    using GridView                 = typename GridImpl::Traits::LeafGridView;
    explicit IgaIdSet(GridView const& g) : gridView_(&g) {}

    template <class Entity>
    [[nodiscard]] IdType id(const Entity& entity) const {
      return this->offset(Entity::codimension)+gridView_->indexSet().index(entity);
    }

    [[nodiscard]] IdType subId(const typename GridImpl::Traits::template Codim<0>::Entity& entity, int i, unsigned int codim) const {
      return this->offset(codim)+gridView_->indexSet().subIndex(entity, i, codim);
    }

  private:

    [[nodiscard]] IdType offset(const int codim) const {
      IdType sum=0;
    for (int i = 0; i < codim; ++i)
      sum+= gridView_->size(i);
    return sum;
    }
    GridView const* gridView_;
  };
}  // namespace Dune::IGA

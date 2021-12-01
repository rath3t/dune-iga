//
// Created by lex on 16.11.21.
//

#pragma once

namespace Dune::IGA {
  template <typename GridImpl>
  class IgaIdSet {
  public:
    using IdType = int;

    explicit IgaIdSet(const GridImpl& grid) : grid_{&grid} {}

    template <class Entity>
    [[nodiscard]] IdType id(const Entity& entity) const {
      return this->offset(Entity::codimension)+grid_->leafGridView().indexSet().index(entity);
    }

    [[nodiscard]] IdType subId(const typename GridImpl::Traits::template Codim<0>::Entity& entity, int i, unsigned int codim) const {
      return this->offset(codim)+grid_->leafGridView().indexSet().subIndex(entity, i, codim);
    }

  private:

    [[nodiscard]] IdType offset(const int codim) const {
      IdType sum=0;
    for (int i = 0; i < codim; ++i)
      sum+= grid_->size(i);
    return sum;
    }
    const GridImpl* grid_;
  };
}  // namespace Dune::IGA

//
// Created by lex on 16.11.21.
//

#pragma once

namespace Dune::IGA
{
  template<typename IgaGridImpl>
  class IgaIdSet
  {
  public:
    using IdType = std::size_t;

    IgaIdSet(const IgaGridImpl& grid) : grid_{&grid},  indexSet_{grid_->leafGridView().indexSet()}{}

    template<class Entity>
    auto id(const Entity& entity) const
    {
      return indexSet_.index(entity);
    }

    auto subId(const typename IgaGridImpl::Traits::template Codim<0>::Entity& entity,
                 int i,
                 unsigned int codim) const
    {
      return indexSet_.subIndex(entity,i,codim);
//      throw std::logic_error("subId not implemented!");
    }
  private:
    const IgaGridImpl* grid_;
    const typename IgaGridImpl::Traits::IndexSet indexSet_;
  };
}

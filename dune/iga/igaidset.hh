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

    explicit IgaIdSet(const IgaGridImpl& grid) :  gridView_{grid.leafGridView()}{}

    template<class Entity>
    auto id(const Entity& entity) const
    {
      return gridView_.indexSet().index(entity);
    }

    auto subId(const typename IgaGridImpl::Traits::template Codim<0>::Entity& entity,
                 int i,
                 unsigned int codim) const
    {
      return gridView_.indexSet().subIndex(entity,i,codim);
//      throw std::logic_error("subId not implemented!");
    }
  private:
//    const IgaGridImpl* grid_;
    const typename IgaGridImpl::Traits::GridView gridView_;
//    const typename IgaGridImpl::Traits::IndexSet indexSet_;
  };
}

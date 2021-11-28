//
// Created by lex on 16.11.21.
//

#pragma once

namespace Dune::IGA {
  template <typename IgaGridImpl>
  class IgaIdSet {
  public:
    using IdType = int;

    explicit IgaIdSet(const IgaGridImpl& grid) : gridView_{grid.leafGridView()} {}

    template <class Entity>
    int id(const Entity& entity) const {
      return gridView_.indexSet().index(entity);
    }

    auto subId(const typename IgaGridImpl::Traits::template Codim<0>::Entity& entity, int i, unsigned int codim) const {
      return gridView_.indexSet().subIndex(entity, i, codim);
    }

  private:
    const typename IgaGridImpl::Traits::GridView gridView_;
  };
}  // namespace Dune::IGA

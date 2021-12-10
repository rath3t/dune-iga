// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <ranges>

#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/nurbsgridentity.hh>
#include <dune/iga/nurbsgridindexsets.hh>
#include <dune/iga/nurbsgridleafiterator.hh>
#include <dune/iga/nurbspatch.hh>

namespace Dune::IGA {

  /** \brief Collect several types associated to OneDGrid LeafGridViews */
  template <class GridImp>
  struct NurbsLeafGridViewTraits {
    using Grid        = GridImp;
    using IndexSet    = NURBSGridLeafIndexSet<GridImp>;
    using GridViewImp = NURBSLeafGridView<GridImp>;
    //    typedef NurbsLeafGridViewTraits<GridImp> GridViewImp;
    //
    //    /** \brief type of the grid */
    //    typedef typename std::remove_const<GridImp>::type Grid;
    //
    //    /** \brief type of the intersection iterator */
    typedef typename GridImp ::Traits ::LeafIntersectionIterator LeafIntersectionIterator;
    //    typedef typename Grid ::Traits ::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename GridImp ::Traits ::LeafIntersectionIterator IntersectionIterator;
    //    typedef typename Grid ::Traits ::HierarchicIterator HierarchicIterator;
    typedef typename GridImp ::Traits ::CollectiveCommunication CollectiveCommunication;
    typedef typename GridImp ::Traits ::LeafIntersection Intersection;

    template <int cd>
    struct Codim {
      typedef typename GridImp::Traits ::template Codim<cd>::template Partition<All_Partition>::LeafIterator Iterator;

      typedef typename GridImp::Traits::template Codim<cd>::Entity Entity;
      //      typedef typename Grid::Traits::template Codim<cd>::EntitySeed EntitySeed;

      typedef typename GridImp::Traits::template Codim<cd>::Geometry Geometry;
      typedef typename GridImp::Traits::template Codim<cd>::LocalGeometry LocalGeometry;

      //      /** \brief Define types needed to iterate over entities of a given partition type */
      template <PartitionIteratorType pit>
      struct Partition {
        /** \brief iterator over a given codim and partition type */
        typedef typename GridImp::Traits::template Codim<cd>::template Partition<pit>::LeafIterator Iterator;
      };
    };
    //
    static constexpr bool conforming = true;
  };

  template <int codim, int dim, typename GridImpl>
  class NURBSGridEntity;

  template <typename GridImpl, int griddim, std::integral auto... codim>
  std::tuple<std::vector<NURBSGridEntity<codim, griddim, GridImpl>>...> gridEntityTupleGenerator(
      std::integer_sequence<std::common_type_t<decltype(codim)...>, codim...>);

  template <typename GridImpl>
  const auto &elements(const NURBSLeafGridView<GridImpl> &gridLeafView);
  template <typename GridImpl>
  auto &elements(NURBSLeafGridView<GridImpl> &gridLeafView);

  /** \brief NURBS grid manager */
  template <typename GridImpl>
  class NURBSLeafGridView {
  public:
    using NurbsGridLinearAlgebraTraits = typename GridImpl::LinearAlgebraTraits;
    using Traits                       = NurbsLeafGridViewTraits<GridImpl>;


    using ControlPointNetType = typename GridImpl::ControlPointNetType;

    template <int codim, int dim, typename GridImpl1>
    friend class NURBSGridEntity;

    using LeafIntersectionIterator = typename Traits::LeafIntersectionIterator;
    using IndexSet                 = typename Traits::IndexSet;
    template <int cd>
    using Codim = typename Traits::template Codim<cd>;
    //    using IntersectionIterator = LeafIntersectionIterator;
    using CollectiveCommunication = typename Traits::CollectiveCommunication;
    using Intersection            = typename Traits::Intersection;
    using IntersectionIterator    = LeafIntersectionIterator;

    using ctype                          = double;
    static constexpr auto dimension      = GridImpl::dimension;
    static constexpr auto dimensionworld = GridImpl::dimensionworld;

    template <int codim, PartitionIteratorType pitype, class GridImp1>
    friend class NURBSGridLeafIterator;

    using Grid = typename Traits::Grid;

    NURBSLeafGridView(const NURBSPatchData<(size_t)dimension, (size_t)dimensionworld, NurbsGridLinearAlgebraTraits> &patchData,
                      const GridImpl &grid)
        : NURBSLeafGridView(patchData.knotSpans, patchData.controlPoints, patchData.degree, grid) {}

    NURBSLeafGridView(const std::array<std::vector<double>, dimension> &knotSpans, const ControlPointNetType &controlPoints,
                      const std::array<int, dimension> order, const GridImpl &grid)
        : NURBSpatch_(
            std::make_shared<NURBSPatch<dimension, dimensionworld, NurbsGridLinearAlgebraTraits>>(knotSpans, controlPoints, order)),
          grid_{&grid},
          entityVector_{
              std::make_shared<decltype(gridEntityTupleGenerator<Grid, dimension>(std::make_integer_sequence<int, dimension + 1>()))>()},
          indexSet_{*this} {
      Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<dimension + 1>()), [&](const auto i) {
        std::get<i>(*entityVector_.get()).reserve(NURBSpatch_->size(i));
        for (unsigned int j = 0; j < NURBSpatch_->size(i); ++j)
          std::get<i>(*entityVector_.get()).emplace_back(*this, j);
      });
    }

    template <int codim>
    typename Codim<codim>::Entity &getEntity(unsigned int directIndex) const {
      // need to be rewrite for other codims
      if constexpr (codim == 0)  // elements
        return (std::get<0>(*entityVector_.get()).at(directIndex));
      else if constexpr (codim == dimension)  // vertices
        return (std::get<dimension>(*entityVector_.get()).at(directIndex));
      else if constexpr (dimension - codim == 1)  // edges
        return (std::get<dimension - 1>(*entityVector_.get()).at(directIndex));
      else if constexpr (dimension - codim == 2)  // surface
        return (std::get<dimension - 2>(*entityVector_.get()).at(directIndex));
      else
        throw std::logic_error("Your requested entity type does not exist.");
    }

    /** \brief obtain collective communication object */
    const auto &comm() const { return grid().comm(); }

    template <class Entity>
    bool contains(const Entity &e) const {
      return e.NURBSGridView_ == this;
    }

    /** \brief obtain collective communication object */
    const auto &grid() const { return *grid_; }

    [[nodiscard]] int size(const GeometryType &type) const {
      if (type == Dune::GeometryTypes::vertex || type == Dune::GeometryTypes::cube(1) || type == Dune::GeometryTypes::cube(2)
          || type == Dune::GeometryTypes::cube(3))
        return this->size(dimension - type.dim());
      else
        return 0;
    }

    const auto &getPatchData() const { return *NURBSpatch_->getPatchData(); }

    auto getPreBasis() { return Dune::Functions::BasisFactory::nurbs<dimension>(*NURBSpatch_->getPatchData()); }

    template <int cd>
    typename Codim<cd>::Iterator begin() const {
      return typename Codim<cd>::Iterator(std::get<cd>(*entityVector_.get()).cbegin());
    }

    template <int cd>
    typename Codim<cd>::Iterator end() const {
      return typename Codim<cd>::Iterator(std::get<cd>(*entityVector_.get()).cend());
    }

    LeafIntersectionIterator ibegin(const typename Codim<0UL>::Entity &entity) const {
      return entity.ibegin(level_);
    }

    LeafIntersectionIterator iend(const typename Codim<0UL>::Entity &entity) const {
      return entity.iend(level_);
    }

    template <int cd, PartitionIteratorType piType>
    typename Codim<cd>::template Partition<piType>::Iterator begin() const {
      if (piType != Ghost_Partition)
        return typename Codim<cd>::template Partition<piType>::Iterator(this->template begin<cd>());
      else
        return typename Codim<cd>::template Partition<piType>::Iterator(this->template end<cd>());
    }

    template <int cd, PartitionIteratorType piType>
    typename Codim<cd>::template Partition<piType>::Iterator end() const {
      return typename Codim<cd>::template Partition<piType>::Iterator(this->template end<cd>());
    }

    const IndexSet &indexSet() const { return indexSet_; }
    [[nodiscard]] int overlapSize(int codim) const { return 0; }
    [[nodiscard]] int ghostSize(int codim) const { return 0; }
    auto size(int codim) const {
      assert(codim <= 3 && codim >= 0);
      if (codim == 0)
        return std::get<0>(*entityVector_.get()).size();
      else if (codim == 1)
        return std::get<1>(*entityVector_.get()).size();
      if constexpr (dimension > 1)
        if (codim == 2) return std::get<2>(*entityVector_.get()).size();
      if constexpr (dimension > 2)
        if (codim == 3) return std::get<3>(*entityVector_.get()).size();
      __builtin_unreachable();
    }

  private:
    friend class NURBSGridLeafIndexSet<GridImpl>;
    friend const auto &elements<GridImpl>(const NURBSLeafGridView<GridImpl> &gridLeafView);
    friend auto &elements<GridImpl>(NURBSLeafGridView<GridImpl> &gridLeafView);
    template <typename GridImp1>
    friend class NURBSintersection;
    std::shared_ptr<NURBSPatch<dimension, dimensionworld, NurbsGridLinearAlgebraTraits>> NURBSpatch_;
    IndexSet indexSet_;
    const GridImpl *grid_;
    int level_{};
    using EntityVectorType = decltype(gridEntityTupleGenerator<GridImpl, dimension>(std::make_integer_sequence<int, dimension + 1>()));
    std::shared_ptr<EntityVectorType> entityVector_{};
  };

  template <typename GridImpl>
  const auto &elements(const NURBSLeafGridView<GridImpl> &gridLeafView) {
    return std::get<0>(*gridLeafView.entityVector_.get());
  }

  template <typename GridImpl>
  auto &elements(NURBSLeafGridView<GridImpl> &gridLeafView) {
    return std::get<0>(*gridLeafView.entityVector_.get());
  }

  template <typename GridImpl1, typename ElementEntity>
  auto &intersections(const typename GridImpl1::GridView &gridLeafView, const ElementEntity &e) {
    return *e.intersections_.get();
  }
}  // namespace Dune::IGA

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_NURBSGRIDENTITY_HH
#define DUNE_IGA_NURBSGRIDENTITY_HH

#include "nurbsintersection.hh"

#include <array>
#include <numeric>

#include <dune/common/fvector.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/iga/nurbsleafgridview.hh>
#include <dune/iga/nurbspatch.hh>

/** \file
 * \brief The NURBSGridEntity class
 */

namespace Dune::IGA {

  template <int codim,int dim, typename GridImpl>
  class NURBSGridEntity {
  public:

    using NurbsGridLinearAlgebraTraits = typename GridImpl::NurbsGridLinearAlgebraTraits;
    static constexpr auto mydim    = GridImpl::dimension - codim;
    static constexpr int  codimension    = codim;
    static constexpr int  dimension    = GridImpl::dimension;
    static constexpr auto dimworld = GridImpl::dimensionworld;
    using Geometry                 = NURBSGeometry<mydim, dimworld, GridImpl>;
    using EntitySeed = typename GridImpl::Traits::template Codim<codim>::EntitySeed;
    using GridView = typename GridImpl::GridView;
    //! Default Constructor
    NURBSGridEntity() = default;

    NURBSGridEntity(const GridView& gridView, unsigned int directIndex)
        : NURBSGridView_(&gridView),
          directIndex_(directIndex),
          parType_{
              (NURBSGridView_->NURBSpatch_->isBorderElement(directIndex_) ? PartitionType::BorderEntity : PartitionType::InteriorEntity)} {}

    using LocalIntersectionGeometry = typename GridView::Traits::template Codim<1>::LocalGeometry;
    //! Geometry of this entity
    typename GridView::template Codim<codim>::Geometry geometry() const {
      //      std::cerr<< "Error geometry not implemented yet for geometries of codim!=0"<<std::endl;
      return NURBSGridView_->NURBSpatch_->template geometry<codim>(directIndex_);
    }

    [[nodiscard]] unsigned int getIndex() const { return directIndex_; }

    [[nodiscard]] unsigned int subEntities(unsigned int codim1) const {
      return (mydim < codim1 ? 0 : Dune::binomial(static_cast<unsigned int>(mydim), codim1) << codim1);
    }

    [[nodiscard]] auto type() const { return GeometryTypes::cube(GridView::dimension - codim); }
    [[nodiscard]] int level() const { return 0; }
    [[nodiscard]] EntitySeed seed() const { EntitySeed e; e.index= directIndex_; return e; }
    [[nodiscard]] PartitionType partitionType() const { return parType_; }

    LocalIntersectionGeometry localGeometry() const {
      return NURBSGridView_->NURBSpatch_->template geometry<codim,true>(directIndex_);
    }

    auto operator<=>(const NURBSGridEntity&) const = default;

  private:
    friend GridView;
    const GridView* NURBSGridView_{nullptr};
    unsigned int directIndex_{};
    PartitionType parType_{};

  };  // end of OneDGridEntity codim = 0

  /** \brief
   *.
   */
  template < int dim,typename GridImpl>
  class NURBSGridEntity<0, dim,GridImpl> {
  public:
    static constexpr auto mydim    = GridImpl::dimension;
    static constexpr auto dimworld = GridImpl::dimensionworld;
    static constexpr int codimension    = 0;
    static constexpr int dimension    = GridImpl::dimension;
    using Geometry                 = typename GridImpl::Traits::template Codim<0>::Geometry;//NURBSGeometry<mydim, dimworld, dimension, typename GridViewImp::NurbsGridLinearAlgebraTraits>;
    using Intersection             = typename GridImpl::Traits::LeafIntersection;
    using IntersectionIterator     = typename GridImpl::Traits::LeafIntersectionIterator;
    using HierarchicIterator     = typename GridImpl::Traits::HierarchicIterator;
    using EntitySeed = typename GridImpl::Traits::template Codim<0>::EntitySeed;
    using LocalGeometry = typename GridImpl::template Codim<1>::LocalGeometry;
    using GridView = typename GridImpl::Traits::LeafGridView;
    NURBSGridEntity() = default;

    NURBSGridEntity(const GridView& gridView, unsigned int directIndex)
        : NURBSGridView_(&gridView),
          directIndex_(directIndex),
          parType_{
              (NURBSGridView_->NURBSpatch_->isBorderElement(directIndex_) ? PartitionType::BorderEntity : PartitionType::InteriorEntity)}
    {
      intersections_= std::make_shared<std::vector<Intersection>>();
      intersections_->reserve(this->subEntities(1));
      for (int innerLocalIndex = 0,outerLocalIndex=1; innerLocalIndex < this->subEntities(1); ++innerLocalIndex) {

        auto multiIndex = NURBSGridView_->NURBSpatch_->elementNet_->directToMultiIndex(directIndex_);
        multiIndex[static_cast<int>(std::floor(innerLocalIndex/2))]+=((innerLocalIndex%2) ? 1 : -1);
        auto directOuterIndex = NURBSGridView_->NURBSpatch_->elementNet_->index(multiIndex);
        intersections_->emplace_back(innerLocalIndex,outerLocalIndex,directIndex_,directOuterIndex,*NURBSGridView_);
        outerLocalIndex += ((innerLocalIndex-1)%2) ? -1 : 3;
      }

    }

    //! Geometry of this entity
    typename GridImpl::Traits::template Codim<0>::Geometry geometry() const {
      return NURBSGridView_->NURBSpatch_->template geometry<0>(directIndex_);
    }


    [[nodiscard]] unsigned int getIndex() const { return directIndex_; }
    [[nodiscard]] bool hasFather() const { return false; }
    [[nodiscard]] EntitySeed seed() const  { EntitySeed e; e.index= directIndex_; return e; }
    LocalGeometry geometryInFather()const{ }
    [[nodiscard]] NURBSGridEntity father() const { throw std::logic_error("father function not implemented."); }

    [[nodiscard]] unsigned int subEntities(unsigned int codim1) const {
      return (mydim < codim1 ? 0 : Dune::binomial(static_cast<unsigned int>(mydim), codim1) << codim1);
    }
    bool hasBoundaryIntersections()const  { return false;}
    template <int codimSub>
    typename GridImpl::Traits::template Codim<codimSub>::Entity subEntity(int i) const {
      if constexpr (codimSub == 0) {
        assert(i == 0);
        return *this;
      } else if constexpr (codimSub == mydim)  // vertices from elements
      {
        auto globalIndex = NURBSGridView_->NURBSpatch_->getGlobalVertexIndexFromElementIndex(directIndex_, i);
        return NURBSGridView_->template getEntity<codimSub>( globalIndex);
      } else if constexpr (mydim - codimSub == 1)  // edges from elements
      {
        auto globalIndex = NURBSGridView_->NURBSpatch_->getGlobalEdgeIndexFromElementIndex(directIndex_, i);
        return NURBSGridView_->template getEntity<codimSub>( globalIndex);
      }

      throw std::logic_error("The requested subentity codim combination is not supported ");
    }

    bool isNew() const { return false;}
    bool mightVanish() const { return false;}

    IntersectionIterator ibegin([[maybe_unused]] int lvl) const { return IntersectionIterator(intersections_->begin()); }
    HierarchicIterator hbegin([[maybe_unused]] int lvl) const { return  HierarchicIterator(*this); }

    IntersectionIterator iend([[maybe_unused]] int lvl) const { return IntersectionIterator(intersections_->end()); }
    HierarchicIterator hend([[maybe_unused]] int lvl) const {return HierarchicIterator(*this); }

    [[nodiscard]] bool isLeaf() const { return true; }
    [[nodiscard]] auto type() const { return GeometryTypes::cube(mydim); }
    [[nodiscard]] int level() const { return 0; }
    [[nodiscard]] PartitionType partitionType() const { return parType_; }

    auto operator<=>(const NURBSGridEntity&) const = default;
//    auto operator<=>(const NURBSGridEntity&) const = default;

  private:
    friend GridView;
    template <typename GridImpl1, typename ElementEntity>
    friend auto& intersections(const typename GridImpl1::GridView& gridLeafView, const ElementEntity& e);
    const GridView* NURBSGridView_{nullptr};
    std::shared_ptr<std::vector<Intersection>> intersections_;
    unsigned int directIndex_{};
    PartitionType parType_{};

  };  // end of OneDGridEntity codim = 0

  template <std::integral auto codim, std::integral auto dim, typename GridImp>
  auto referenceElement(const NURBSGridEntity<codim,dim, GridImp>& e) {
    return Dune::ReferenceElements<typename GridImp::ctype, NURBSGridEntity<codim, dim,GridImp>::mydim>::cube();
  };

}  // namespace Dune::IGA

#endif  // DUNE_IGA_NURBSGRIDENTITY_HH

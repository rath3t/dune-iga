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

  template <std::integral auto codim, typename GridViewImp>
  class NURBSGridEntity {
  public:
    static constexpr auto mydim    = GridViewImp::dimension - codim;
    static constexpr auto dimworld = GridViewImp::dimensionworld;
    using Geometry                 = NURBSGeometry<mydim, dimworld, mydim, typename GridViewImp::NurbsGridLinearAlgebraTraits>;
    //! Default Constructor
    NURBSGridEntity() = default;

    NURBSGridEntity(const GridViewImp& gridView, unsigned int directIndex)
        : NURBSGridView_(&gridView),
          directIndex_(directIndex),
          parType_{
              (NURBSGridView_->NURBSpatch_->isBorderElement(directIndex_) ? PartitionType::BorderEntity : PartitionType::InteriorEntity)} {}

    using LocalIntersectionGeometry = typename GridViewImp::Traits::template Codim<1>::LocalGeometry;
    //! Geometry of this entity
    typename GridViewImp::template Codim<codim>::Geometry geometry() const {
      //      std::cerr<< "Error geometry not implemented yet for geometries of codim!=0"<<std::endl;
      return NURBSGridView_->NURBSpatch_->template geometry<codim>(directIndex_);
    }

    [[nodiscard]] unsigned int getIndex() const { return directIndex_; }

    [[nodiscard]] unsigned int subEntities(unsigned int codim1) const {
      return (mydim < codim1 ? 0 : Dune::binomial(static_cast<unsigned int>(mydim), codim1) << codim1);
    }

    [[nodiscard]] auto type() const { return GeometryTypes::cube(GridViewImp::dimension - codim); }
    [[nodiscard]] int level() const { return 0; }
    [[nodiscard]] PartitionType partitionType() const { return parType_; }

    LocalIntersectionGeometry localGeometry() const {
      return NURBSGridView_->NURBSpatch_->template geometry<codim,true>(directIndex_);
    }

    auto operator<=>(const NURBSGridEntity&) const = default;

  private:
    friend GridViewImp;
    const GridViewImp* NURBSGridView_{nullptr};
    unsigned int directIndex_{};
    PartitionType parType_{};

  };  // end of OneDGridEntity codim = 0

  /** \brief
   *.
   */
  template <typename GridViewImp>
  class NURBSGridEntity<0, GridViewImp> {
  public:
    static constexpr auto mydim    = GridViewImp::dimension;
    static constexpr auto dimworld = GridViewImp::dimensionworld;
    using Geometry                 = NURBSGeometry<mydim, dimworld, mydim, typename GridViewImp::NurbsGridLinearAlgebraTraits>;
    using Intersection             = NURBSintersection<mydim - 1, GridViewImp>;
    using IntersectionIterator     = typename Intersection::Iterator;


    NURBSGridEntity() = default;

    NURBSGridEntity(const GridViewImp& gridView, unsigned int directIndex)
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
    typename GridViewImp::template Codim<0>::Geometry geometry() const {
      return NURBSGridView_->NURBSpatch_->template geometry<0UL>(directIndex_);
    }


    [[nodiscard]] unsigned int getIndex() const { return directIndex_; }

    [[nodiscard]] unsigned int subEntities(unsigned int codim1) const {
      return (mydim < codim1 ? 0 : Dune::binomial(static_cast<unsigned int>(mydim), codim1) << codim1);
    }

    template <int codimSub>
    typename GridViewImp::template Codim<codimSub>::Entity subEntity(int i) const {
      if constexpr (codimSub == 0) {
        assert(i == 0);
        return *this;
      } else if constexpr (codimSub == mydim)  // vertices from elements
      {
        auto globalIndex = NURBSGridView_->NURBSpatch_->getGlobalVertexIndexFromElementIndex(directIndex_, i);
        return typename GridViewImp::template Codim<codimSub>::Entity(*NURBSGridView_, globalIndex);
      } else if constexpr (mydim - codimSub == 1)  // edges from elements
      {
        auto globalIndex = NURBSGridView_->NURBSpatch_->getGlobalEdgeIndexFromElementIndex(directIndex_, i);
        return typename GridViewImp::template Codim<codimSub>::Entity(*NURBSGridView_, globalIndex);
      }

      throw std::logic_error("The requested subentity codim combination is not supported ");
    }

    IntersectionIterator ibegin([[maybe_unused]] int lvl) const { return IntersectionIterator(intersections_->begin()); }

    IntersectionIterator iend([[maybe_unused]] int lvl) const { return IntersectionIterator(intersections_->end()); }

    [[nodiscard]] bool isLeaf() const { return true; }
    [[nodiscard]] auto type() const { return GeometryTypes::cube(mydim); }
    [[nodiscard]] int level() const { return 0; }
    [[nodiscard]] PartitionType partitionType() const { return parType_; }

    auto operator<=>(const NURBSGridEntity&) const = default;

  private:
    friend GridViewImp;
    template <typename GridImpl, typename ElementEntity>
    friend auto& intersections(const NURBSLeafGridView<GridImpl>& gridLeafView, const ElementEntity& e);
    const GridViewImp* NURBSGridView_{nullptr};
    std::shared_ptr<std::vector<Intersection>> intersections_;
    unsigned int directIndex_{};
    PartitionType parType_{};

  };  // end of OneDGridEntity codim = 0

  template <std::integral auto codim, typename GridViewImp>
  auto referenceElement(const NURBSGridEntity<codim, GridViewImp>& e) {
    return Dune::ReferenceElements<typename GridViewImp::ctype, NURBSGridEntity<codim, GridViewImp>::mydim>::cube();
  };

}  // namespace Dune::IGA

#endif  // DUNE_IGA_NURBSGRIDENTITY_HH

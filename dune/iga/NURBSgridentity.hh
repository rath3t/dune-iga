// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_NURBSGRIDENTITY_HH
#define DUNE_IGA_NURBSGRIDENTITY_HH

#include <array>
#include <numeric>

#include <dune/common/fvector.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/iga/NURBSleafgridview.hh>
#include <dune/iga/NURBSpatch.hh>

/** \file
 * \brief The NURBSGridEntity class
 */

namespace Dune::IGA {

  template <std::integral auto codim, typename GridViewImp>
  class NURBSGridEntity {
  public:
    static constexpr auto mydim = GridViewImp::dimension - codim;
    using Geometry
        = NURBSGeometry<mydim, GridViewImp::dimensionworld, GridViewImp::dimension, typename GridViewImp::NurbsGridLinearAlgebraTraits>;
    //! Default Constructor
    NURBSGridEntity() = default;

    NURBSGridEntity(const GridViewImp& gridView, unsigned int directIndex)
        : NURBSGridView_(&gridView),
          directIndex_(directIndex),
          parType_{
              (NURBSGridView_->NURBSpatch_->isBorderElement(directIndex_) ? PartitionType::BorderEntity : PartitionType::InteriorEntity)}
    // TODO find out if boundary ent

    {}

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
    int level() const { return 0; }
    [[nodiscard]] PartitionType partitionType() const { return parType_; }

    //    bool equals(const BSplinePatchEntity& other) const
    //    {
    //      return target_ == other.target_;
    //    }

    auto operator<=>(const NURBSGridEntity&) const = default;

  private:
    friend GridViewImp;
    const GridViewImp* NURBSGridView_{nullptr};
    unsigned int directIndex_;
    PartitionType parType_;

  };  // end of OneDGridEntity codim = 0

  /** \brief
   *.
   */
  template <typename GridViewImp>
  class NURBSGridEntity<0, GridViewImp> {
  public:
    static constexpr auto mydim = GridViewImp::dimension;
    using Geometry
        = NURBSGeometry<mydim, GridViewImp::dimensionworld, GridViewImp::dimension, typename GridViewImp::NurbsGridLinearAlgebraTraits>;
    //! Default Constructor
    NURBSGridEntity()                       = default;
    NURBSGridEntity(const NURBSGridEntity&) = default;             // NURBSGridEntityopy constructor
    NURBSGridEntity(NURBSGridEntity&&)       noexcept = default;             // Move constructor
    NURBSGridEntity& operator=(const NURBSGridEntity&) = default;  // NURBSGridEntityopy assignment operator
    NURBSGridEntity& operator=(NURBSGridEntity&&)  noexcept = default;       // Move assignment operator
    ~NURBSGridEntity()                          = default;       // Destructor

    NURBSGridEntity(const GridViewImp& gridView, unsigned int directIndex)
        : NURBSGridView_(&gridView),
          directIndex_(directIndex),
          parType_{
              (NURBSGridView_->NURBSpatch_->isBorderElement(directIndex_) ? PartitionType::BorderEntity : PartitionType::InteriorEntity)} {}

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
      } else if constexpr (codimSub == mydim )  // vertices from elements
      {
        auto globalIndex = NURBSGridView_->NURBSpatch_->getGlobalVertexIndexFromElementIndex(directIndex_, i);
        return typename GridViewImp::template Codim<codimSub>::Entity(*NURBSGridView_, globalIndex);
      }
      else if constexpr (mydim-codimSub == 1)  //edges from elements
      {
        auto globalIndex = NURBSGridView_->NURBSpatch_->getGlobalEdgeIndexFromElementIndex(directIndex_, i);
        return typename GridViewImp::template Codim<codimSub>::Entity(*NURBSGridView_, globalIndex);
      }

    throw std::logic_error("The requested subentity codim combination is not supported ");

    }

    bool isLeaf() const { return true;}
    [[nodiscard]] auto type() const { return GeometryTypes::cube(mydim); }
    int level() const { return 0; }
    [[nodiscard]] PartitionType partitionType() const { return parType_; }

    auto operator<=>(const NURBSGridEntity&) const = default;
  private:
//    friend bool operator==(const NURBSGridEntity<0,GridViewImp>& l,const NURBSGridEntity<0,GridViewImp>& r) ;
    friend GridViewImp;
    const GridViewImp* NURBSGridView_{nullptr};
    unsigned int directIndex_;
    PartitionType parType_;

  };  // end of OneDGridEntity codim = 0

  template <std::integral auto codim, typename GridViewImp>
  auto referenceElement(const NURBSGridEntity<codim,GridViewImp>& e) {
    return Dune::ReferenceElements<typename GridViewImp::ctype, NURBSGridEntity<codim,GridViewImp>::mydim>::cube();
  };

}  // namespace Dune::IGA

#endif  // DUNE_IGA_NURBSGRIDENTITY_HH

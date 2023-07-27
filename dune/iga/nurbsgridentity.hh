// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef DUNE_IGA_NURBSGRIDENTITY_HH
#define DUNE_IGA_NURBSGRIDENTITY_HH

#include "dune/iga/nurbsintersection.hh"
#include "dune/iga/nurbsleafgridview.hh"
#include "dune/iga/nurbspatch.hh"
#include <dune/common/fvector.hh>
#include <dune/grid/common/gridenums.hh>

/** \file
 * \brief The NURBSGridEntity class
 */
namespace Dune::IGA {

  template <int codim, int dim, typename GridImpl>
  class NURBSGridEntity : public EntityDefaultImplementation<codim, dim, GridImpl, NURBSGridEntity> {
   public:
    static constexpr auto mydimension = GridImpl::dimension - codim;
    static constexpr int codimension  = codim;
    static constexpr int dimension    = GridImpl::dimension;
    static constexpr auto dimworld    = GridImpl::dimensionworld;
    using Geometry                    = NURBSGeometry<mydimension, dimworld, GridImpl>;
    using EntitySeed                  = typename GridImpl::Traits::template Codim<codim>::EntitySeed;
    using GridView                    = typename GridImpl::GridView;
    NURBSGridEntity()                 = default;
    NURBSGridEntity(const GridView& gridView, unsigned int directIndex, unsigned int patchID)
        : NURBSGridView_(&gridView), directIndex_(directIndex), patchID_{patchID} {}

    using LocalIntersectionGeometry = typename GridView::Traits::template Codim<1>::LocalGeometry;
    //! Geometry of this entity
    typename GridView::template Codim<codim>::Geometry geometry() const {
      return NURBSGridView_->impl().getPatch(patchID_).template geometry<codim>(directIndex_);
    }

    [[nodiscard]] unsigned int getIndex() const { return directIndex_; }

    [[nodiscard]] unsigned int subEntities(unsigned int codim1) const {
      return (mydimension < codim1 ? 0 : Dune::binomial(static_cast<unsigned int>(mydimension), codim1) << codim1);
    }

    [[nodiscard]] auto type() const { return GeometryTypes::cube(GridView::dimension - codim); }
    [[nodiscard]] int level() const { return 0; }
    [[nodiscard]] EntitySeed seed() const {
      EntitySeed e;
      e.impl().index_ = directIndex_;
      e.impl().valid_ = true;
      return e;
    }
    [[nodiscard]] PartitionType partitionType() const { return parType_; }

    auto operator<=>(const NURBSGridEntity&) const = default;

    auto equals(const NURBSGridEntity& r) const {
      return (r.NURBSGridView_ == NURBSGridView_ && directIndex_ == r.directIndex_ && parType_ == r.parType_);
    }

   private:
    friend GridView;
    friend NURBSLeafGridView<GridImpl>;
    void updateGridView(const GridView& other) { NURBSGridView_ = &other; }
    const GridView* NURBSGridView_{nullptr};
    unsigned int directIndex_{};
    unsigned int patchID_{};
    PartitionType parType_{PartitionType::InteriorEntity};
  };

  template <int dim, typename GridImpl>
  class NURBSGridEntity<0, dim, GridImpl> : public EntityDefaultImplementation<0, dim, GridImpl, NURBSGridEntity> {
   public:
    static constexpr auto mydimension = GridImpl::dimension;
    static constexpr auto dimworld    = GridImpl::dimensionworld;
    static constexpr int codimension  = 0;
    static constexpr int dimension    = GridImpl::dimension;
    using Geometry =
        typename GridImpl::Traits::template Codim<0>::Geometry;  // NURBSGeometry<mydimension, dimworld, dimension,
                                                                 // typename GridViewImp::LinearAlgebraTraits>;
    using Intersection         = typename GridImpl::Traits::LeafIntersection;
    using IntersectionIterator = typename GridImpl::Traits::LeafIntersectionIterator;
    using HierarchicIterator   = typename GridImpl::Traits::HierarchicIterator;
    using EntitySeed           = typename GridImpl::Traits::template Codim<0>::EntitySeed;
    using LocalGeometry        = typename GridImpl::template Codim<0>::LocalGeometry;
    using GridView             = typename GridImpl::Traits::LeafGridView;
    NURBSGridEntity()          = default;

    NURBSGridEntity(const GridView& gridView, unsigned int directIndex, unsigned int patchID)
        : NURBSGridView_(&gridView),
          directIndex_(directIndex),
          patchID_{patchID},
          parType_{PartitionType::InteriorEntity},
          trimFlag{NURBSGridView_->impl().getPatch(patchID_).getTrimFlag(directIndex_)},
          intersections_{std::make_shared<std::vector<Intersection>>()}  // FIXME why shared?
    {
      intersections_->reserve(this->subEntities(1));

      for (int innerLocalIndex = 0, outerLocalIndex = 1; innerLocalIndex < this->subEntities(1); ++innerLocalIndex) {
        const auto& eleNet    = NURBSGridView_->impl().getPatch(patchID_).elementNet_;
        auto nurbsDirectIndex = NURBSGridView_->impl().getPatch(patchID_).template getDirectIndex<0>(directIndex_);
        auto multiIndex       = eleNet->directToMultiIndex(nurbsDirectIndex);
        multiIndex[static_cast<int>(std::floor(innerLocalIndex / 2))]
            += ((innerLocalIndex % 2)
                    ? 1
                    : Impl::noNeighbor);  // increase the multiIndex depending on where the outer element should lie
        auto directOuterIndex = (eleNet->isValid(multiIndex)) ? eleNet->index(multiIndex) : Impl::noNeighbor;
        directOuterIndex      = getRealIndexForOuterIndex(directOuterIndex);
        intersections_->emplace_back(NURBSintersection<GridImpl>(innerLocalIndex, outerLocalIndex, directIndex_,
                                                                 directOuterIndex, *NURBSGridView_));
        outerLocalIndex += ((innerLocalIndex - 1) % 2) ? -1 : 3;
      }
    }

    auto elementSubGrid() const { return NURBSGridView_->impl().getPatch(patchID_).getElementSubGrid(directIndex_); }

    void fillQuadratureRule(Dune::QuadratureRule<double, dim>& vector, const std::optional<int>& p_order = std::nullopt,
                            const QuadratureType::Enum qt = QuadratureType::GaussLegendre) const {
      vector.clear();
      int order = p_order.value_or(mydimension
                                   * (*std::ranges::max_element(NURBSGridView_->impl().getPatchData(patchID_).degree)));
      if constexpr (dim == 2)
        if (isTrimmed()) {
          auto elementRepr = NURBSGridView_->impl().getPatch(patchID_).getElementSubGrid(directIndex_);
          fillQuadratureRuleImpl(vector, *elementRepr, order, qt);
          return;
        }
      const auto& rule = Dune::QuadratureRules<double, mydimension>::rule(this->type(), order, qt);
      vector.insert(vector.end(), rule.begin(), rule.end());
    }

    [[nodiscard]] Dune::QuadratureRule<double, dim> getQuadratureRule(const std::optional<int>& p_order = std::nullopt,
                                                                      const QuadratureType::Enum qt
                                                                      = QuadratureType::GaussLegendre) const {
      Dune::QuadratureRule<double, dim> rule;
      fillQuadratureRule(rule, p_order, qt);
      return rule;
    }

    //! Geometry of this entity
    typename GridImpl::Traits::template Codim<0>::Geometry geometry() const {
      return NURBSGridView_->impl().getPatch(patchID_).template geometry<0>(directIndex_);
    }

    [[nodiscard]] unsigned int getIndex() const { return directIndex_; }
    [[nodiscard]] unsigned int getDirectIndexInPatch() const {
      return NURBSGridView_->impl().getPatch(patchID_).template getDirectIndex<0>(directIndex_);
    }
    [[nodiscard]] bool hasFather() const { return false; }
    [[nodiscard]] EntitySeed seed() const {
      EntitySeed e;
      e.impl().index_ = directIndex_;
      e.impl().valid_ = true;
      return e;
    }
    LocalGeometry geometryInFather() const { throw std::logic_error("geometryInFather function not implemented."); }
    [[nodiscard]] NURBSGridEntity father() const { throw std::logic_error("father function not implemented."); }

    [[nodiscard]] unsigned int subEntities(unsigned int codim1) const {
      return (mydimension < codim1 ? 0 : Dune::binomial(static_cast<unsigned int>(mydimension), codim1) << codim1);
    }

    [[nodiscard]] bool hasBoundaryIntersections() const {
      return std::ranges::any_of(*intersections_,
                                 [](const Intersection& intersection) { return intersection.boundary(); });
    }

    [[nodiscard]] ElementTrimFlag getTrimFlag() const { return trimFlag; }
    [[nodiscard]] bool isTrimmed() const { return trimFlag == ElementTrimFlag::trimmed; }

    template <int codimSub>
    typename GridImpl::Traits::template Codim<codimSub>::Entity subEntity(int i) const {
      if constexpr (codimSub == 0) {
        assert(i == 0);
        return *this;
      } else if constexpr (codimSub == mydimension)  // vertices from elements
      {
        auto globalIndex
            = NURBSGridView_->impl().getPatch(patchID_).getGlobalVertexIndexFromElementIndex(directIndex_, i);
        return NURBSGridView_->impl().template getEntity<codimSub>(globalIndex);
      } else if constexpr (mydimension - codimSub == 1)  // edges from elements
      {
        auto globalIndex
            = NURBSGridView_->impl().getPatch(patchID_).getGlobalEdgeIndexFromElementIndex(directIndex_, i);
        return NURBSGridView_->impl().template getEntity<codimSub>(globalIndex);
      } else if constexpr (mydimension - codimSub == 2)  // surfaces from elements
      {
        auto globalIndex
            = NURBSGridView_->impl().getPatch(patchID_).getGlobalSurfaceIndexFromElementIndex(directIndex_, i);
        return NURBSGridView_->impl().template getEntity<codimSub>(globalIndex);
      }
      throw std::logic_error("The requested subentity codim combination is not supported ");
    }

    [[nodiscard]] bool isNew() const { return false; }
    [[nodiscard]] bool mightVanish() const { return false; }

    IntersectionIterator ibegin([[maybe_unused]] int lvl) const {
      return NURBSGridInterSectionIterator<GridImpl>(intersections_->begin());
    }
    HierarchicIterator hbegin([[maybe_unused]] int lvl) const {
      //      NurbsHierarchicIterator<GridImpl>
      return NurbsHierarchicIterator<GridImpl>(*this);
    }

    IntersectionIterator iend([[maybe_unused]] int lvl) const {
      return NURBSGridInterSectionIterator<GridImpl>(intersections_->end());
    }
    HierarchicIterator hend([[maybe_unused]] int lvl) const { return NurbsHierarchicIterator<GridImpl>(*this); }

    [[nodiscard]] bool isLeaf() const { return true; }

    [[nodiscard]] auto type() const {
      if (trimFlag == ElementTrimFlag::full)
        return GeometryTypes::cube(mydimension);
      else
        return GeometryTypes::none(mydimension);
    }
    [[nodiscard]] int level() const { return 0; }
    [[nodiscard]] PartitionType partitionType() const { return parType_; }

    auto operator<=>(const NURBSGridEntity&) const = default;
    auto equals(const NURBSGridEntity& r) const {
      return (r.NURBSGridView_ == NURBSGridView_ && directIndex_ == r.directIndex_ && parType_ == r.parType_
              && intersections_ == r.intersections_);
    }

   private:
    friend GridView;
    friend NURBSLeafGridView<const GridImpl>;

    void updateGridView(const GridView& other) { NURBSGridView_ = &other; }
    template <typename GridImpl1, typename ElementEntity>
    friend auto& intersections(const typename GridImpl1::GridView& gridLeafView, const ElementEntity& e);
    const GridView* NURBSGridView_{nullptr};
    std::shared_ptr<std::vector<Intersection>> intersections_;
    unsigned int directIndex_{};
    PartitionType parType_{PartitionType::InteriorEntity};
    unsigned int patchID_{};

    ElementTrimFlag trimFlag;

    // Helpers
    int getRealIndexForOuterIndex(int outerIndex) { return outerIndex; }

    int getRealIndexForOuterIndex(int outerIndex) requires(dim == 2) {
      if (outerIndex == Impl::noNeighbor) return Impl::noNeighbor;
      return NURBSGridView_->impl().getPatch(patchID_).template getRealIndexOr<0>(outerIndex, Impl::noNeighbor);
    }

  };  // end of Template Specialization for codim = 0

  template <std::integral auto codim, std::integral auto dim, typename GridImp>
  auto referenceElement(const NURBSGridEntity<codim, dim, const GridImp>& e) {
    return Dune::ReferenceElements<typename GridImp::ctype,
                                   NURBSGridEntity<codim, dim, const GridImp>::mydimension>::cube();
  }

  template <std::integral auto codim, std::integral auto dim, typename GridImp>
  auto referenceElement(const NURBSGridEntity<0, dim, const GridImp>& e) {
    if (e.getTrimFlag() == ElementTrimFlag::trimmed)
      return Dune::ReferenceElements<typename GridImp::ctype,
                                     NURBSGridEntity<0, dim, const GridImp>::mydimension>::none();
    else
      return Dune::ReferenceElements<typename GridImp::ctype,
                                     NURBSGridEntity<0, dim, const GridImp>::mydimension>::cube();
  }

}  // namespace Dune::IGA

#endif  // DUNE_IGA_NURBSGRIDENTITY_HH

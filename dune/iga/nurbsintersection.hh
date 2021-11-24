//
// Created by lex on 23.11.21.
//

#pragma once

namespace Dune::IGA {
  namespace Impl {
    static constexpr int noNeighbor = -1;
  }

  template <std::integral auto mydim, typename GridViewImp>
  class NURBSintersection {
  public:
    using Iterator      = NURBSGridInterSectionIterator<NURBSintersection>;
    using Entity        = typename GridViewImp::template Codim<0>::Entity;
    using Geometry      = typename GridViewImp::template Codim<1>::Geometry;
    using LocalGeometry = typename GridViewImp::Traits::LocalGeometryIntersection;

    NURBSintersection(std::array<unsigned int, 2> directIndices,const GridViewImp& gridView )
        :gridView_{gridView},
        interDirectIndex_{directIndices[0]}, outerDirectIndex_{directIndices[1]} {}
    /** coordinate type */
    typedef double ctype;

    /** \brief Dimension of the cube element */
    static constexpr std::integral auto mydimension = mydim;

    /** \brief Dimension of the world space that the cube element is embedded in*/
    static constexpr std::integral auto coorddimension = GridViewImp::dimensionworld;

    /** \brief Type used for a vector of element coordinates */
    using LocalCoordinate = typename GridViewImp::NurbsGridLinearAlgebraTraits::template FixedVectorType<mydim>;

    /** \brief Type used for a vector of world coordinates */
    using GlobalCoordinate = typename GridViewImp::NurbsGridLinearAlgebraTraits::GlobalCoordinateType;
    [[nodiscard]] bool boundary() const { return outerDirectIndex_ != Impl::noNeighbor; }
    [[nodiscard]] bool neighbor() const { return outerDirectIndex_ != Impl::noNeighbor; }
    [[nodiscard]] bool conforming() const { return true; }
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydim); }

    Entity inside() const {return gridView_->template getEntity<0>(interDirectIndex_);}
    Entity outside() const {
      assert(neighbor()&& "Outer Element does not exist.");
      return gridView_->template getEntity<0>(outerDirectIndex_);
    }
    [[nodiscard]] int indexInInside() const {}
    [[nodiscard]] int indexInOutside() const {}
    LocalGeometry geometryInInside() const {}
    LocalGeometry geometryInOutside() const {}
    Geometry geometry() const {}
    GlobalCoordinate outerNormal(const LocalCoordinate& xi) const {}
    GlobalCoordinate unitOuterNormal(const LocalCoordinate& xi) const {}
    GlobalCoordinate integrationOuterNormal(const LocalCoordinate& xi) const {}
    [[nodiscard]] std::size_t boundarySegmentIndex() const {}

  private:
    GridViewImp* gridView_;
    unsigned int interDirectIndex_;
    unsigned int outerDirectIndex_;
  };
}  // namespace Dune::IGA
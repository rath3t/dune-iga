//
// Created by lex on 23.11.21.
//

#pragma once

namespace Dune::IGA {
  template <std::integral auto mydim, std::integral auto dimworld, std::integral auto griddim,
            NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits>
  class nurbsintersection {
  public:
    /** coordinate type */
    typedef double ctype;

    /** \brief Dimension of the cube element */
    static constexpr std::integral auto mydimension = mydim;

    /** \brief Dimension of the world space that the cube element is embedded in*/
    static constexpr std::integral auto coorddimension = dimworld;

    /** \brief Type used for a vector of element coordinates */
    using LocalCoordinate = typename NurbsGridLinearAlgebraTraits::template FixedVectorType<mydim>;

    /** \brief Type used for a vector of world coordinates */
    using GlobalCoordinate = typename NurbsGridLinearAlgebraTraits::GlobalCoordinateType;
    bool boundary() const { return false; }
    bool neighbor() const { return false; }
    bool conforming() const { return false; }
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydim); }
    Entity inside() const {}
    Entity outside() const {}
int indexInInside() const{}
int indexInOutside() const{}
LocalGeometry geometryInInside() const {}
LocalGeometry geometryInOutside() const {}
Geometry geometry() const {}
GlobalCoordinate outerNormal(const LocalCoordinate& xi) const {}
GlobalCoordinate unitOuterNormal(const LocalCoordinate& xi) const {}

  };
}  // namespace Dune::IGA
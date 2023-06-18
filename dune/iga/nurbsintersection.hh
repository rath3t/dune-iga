// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "nurbsleafgridview.hh"
#include "nurbslocalgeometry.hh"
namespace Dune::IGA {
  namespace Impl {
    static constexpr int noNeighbor = -1;
  }

  template <typename NURBSIntersection>
  class NURBSGridInterSectionIterator;

  template <typename GridImp>
  class NURBSintersection {
   public:
    using Iterator      = NURBSGridInterSectionIterator<NURBSintersection>;
    using Entity        = typename GridImp::Traits::template Codim<0>::Entity;
    using Geometry      = typename GridImp::Traits::template Codim<1>::Geometry;
    using LocalGeometry = typename GridImp::Traits::template Codim<1>::LocalGeometry;
    using GridView      = typename GridImp::Traits::LeafGridView;

    using ctype = typename GridImp::ctype;

    static constexpr std::integral auto mydimension = GridImp::dimension - 1;
    static constexpr std::integral auto dimworld    = GridImp::dimensionworld;
    using LocalCoordinate                           = Dune::FieldVector<ctype, mydimension>;
    using GlobalCoordinate                          = Dune::FieldVector<ctype, dimworld>;

    NURBSintersection() = default;

    NURBSintersection(int innerLocalIndex, int outerLocalIndex, int innerDirectIndex, int outerDirectIndex,
                      const GridView& gridView)
        : gridView_{&gridView},
          innerDirectIndex_{innerDirectIndex},
          outerDirectIndex_{outerDirectIndex},
          innerLocalIndex_{innerLocalIndex},
          outerLocalIndex_{outerLocalIndex} {}

    /** \brief Returns true if the intersection is on the boundary */
    [[nodiscard]] bool boundary() const { return outerDirectIndex_ == Impl::noNeighbor; }

    /** \brief Returns true if the intersection has no outer neighbor */
    [[nodiscard]] bool neighbor() const { return outerDirectIndex_ != Impl::noNeighbor; }
    /** \brief Returns true if the intersection is conforming, i.e. inside the patch this is true but between patches
     * not necessarily */
    [[nodiscard]] bool conforming() const { return true; }

    /** \brief Returns the cube type this intersection, i.e. vertex, edge or quadrilateral */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }

    /** \brief Returns the element entity from which this intersection is constructed */
    Entity inside() const { return gridView_->impl().template getEntity<0>(innerDirectIndex_); }

    /** \brief Returns the element entity which intersects with the inside() element */
    Entity outside() const {
      assert(neighbor() && "Outer Element does not exist.");
      return gridView_->impl().template getEntity<0>(outerDirectIndex_);
    }

    /** \brief Returns the index of the inside element */
    [[nodiscard]] int indexInInside() const { return innerLocalIndex_; }

    /** \brief Returns the index of the outside element, this return Impl::noNeighbor if the outside element does not
     * exist */
    [[nodiscard]] int indexInOutside() const {
      assert(neighbor() && "Outer Element does not exist.");
      return outerLocalIndex_;
    }

    /** \brief Returns the local geometry of the intersection in the coordinates of the element inside */
    LocalGeometry geometryInInside() const {
      return LocalGeometry(NURBSLocalGeometry<mydimension, GridImp::dimension, GridImp>(innerLocalIndex_));
    }

    /** \brief Returns the local geometry of the intersection in the coordinates of the element outside */
    LocalGeometry geometryInOutside() const {
      return LocalGeometry(NURBSLocalGeometry<mydimension, GridImp::dimension, GridImp>(outerLocalIndex_));
    }

    /** \brief Returns the global geometry of the intersection, currently this returns the geometry of the intersection
     * as seen from the inside element. This could differ from the outside element */
    Geometry geometry() const { return inside().template subEntity<1>(innerLocalIndex_).geometry(); }

    /** \brief Returns the normal which lies in the tangent plane of the inside element and is perpendicular on the
     * intersection (edge or surface) */
    [[nodiscard]] GlobalCoordinate outerNormal(const LocalCoordinate& xi) const {
      if constexpr (mydimension == 0) {
        const auto xiInElementCoords        = geometryInInside().global(xi);
        const auto insideJacobianTransposed = inside().geometry().jacobianTransposed(xiInElementCoords);
        return (innerLocalIndex_ == 0) ? -insideJacobianTransposed[0] : insideJacobianTransposed[0];
      } else if constexpr (mydimension == 1) {  // edges
        if constexpr (dimworld == 2) {          // edges in R2
          const auto innerJacobianTransposed = this->geometry().jacobianTransposed(xi)[0];
          switch (innerLocalIndex_) {
            case 0:
            case 3:
              return GlobalCoordinate({-innerJacobianTransposed[1], innerJacobianTransposed[0]});
            case 1:
            case 2:
              return GlobalCoordinate({innerJacobianTransposed[1], -innerJacobianTransposed[0]});
            default:
              __builtin_unreachable();
          }
        } else if constexpr (dimworld == 3)  // edges in R3
        {
          const auto xiInElementCoords        = geometryInInside().global(xi);
          const auto insideJacobianTransposed = inside().geometry().jacobianTransposed(xiInElementCoords);
          auto normal                         = cross(insideJacobianTransposed[0], insideJacobianTransposed[1]);
          switch (innerLocalIndex_) {
            case 0:
              return cross(normal, insideJacobianTransposed[1]);
            case 1:
              return cross(insideJacobianTransposed[1], normal);
            case 2:
              return cross(insideJacobianTransposed[0], normal);
            case 3:
              return cross(normal, insideJacobianTransposed[0]);
            default:
              __builtin_unreachable();
          }
        }
      } else if constexpr (mydimension == 2 && dimworld == 3) {  // surfaces in R3
        const auto innerJacobianTransposed = this->geometry().jacobianTransposed(xi);
        switch (innerLocalIndex_) {
          case 0:
          case 3:
          case 4:
            return cross(innerJacobianTransposed[1], innerJacobianTransposed[0]);
          case 1:
          case 2:
          case 5:
            return cross(innerJacobianTransposed[0], innerJacobianTransposed[1]);
          default:
            __builtin_unreachable();
        }
      }
      __builtin_unreachable();
    }

    /** \brief Same as outerNormal() but with unit length */
    [[nodiscard]] GlobalCoordinate unitOuterNormal(const LocalCoordinate& xi) const {
      auto N = this->outerNormal(xi);
      return N / N.two_norm();
    }

    [[nodiscard]] GlobalCoordinate centerUnitOuterNormal() const {
      if constexpr (mydimension == 0)
        return unitOuterNormal({});
      else
        return unitOuterNormal(0.5);
    }

    /** \brief Same as outerNormal() but with the length of the integration element */
    [[nodiscard]] GlobalCoordinate integrationOuterNormal(const LocalCoordinate& xi) const {
      return this->unitOuterNormal(xi) * this->geometry().integrationElement(xi);
    }

    /** \brief Returns the consecutive index if this intersection lies on the boundary */
    [[nodiscard]] std::size_t boundarySegmentIndex() const {
      assert(boundary());
      auto geomEntity = inside().template subEntity<1>(innerLocalIndex_);
      return gridView_->impl().getPatch(0).patchBoundaryIndex(geomEntity.impl().getIndex());
    }

    auto operator<=>(const NURBSintersection&) const = default;
    bool equals(const NURBSintersection& r) const { return *this == r; }

   private:
    const GridView* gridView_;
    int innerDirectIndex_;
    int innerLocalIndex_;
    int outerDirectIndex_;
    int outerLocalIndex_;
  };
}  // namespace Dune::IGA

PathD // SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "nurbstrimboundary.hh"
#include "nurbstrimmer.hh"
#include "subgridhelpers.hh"

#include <clipper2/clipper.core.h>

#include "dune/iga/geometry/geohelper.hh"
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/virtualrefinement.hh>
#include <dune/grid/uggrid.hh>

    namespace Dune::IGA {

  /** @brief representation of the subgrid of trimmed elements in the reference space */
  template <int dim>
  class TrimmedSubGrid
  {
  public:
    using Point   = Dune::FieldVector<double, dim>;
    using Element = MultiLinearGeometry<double, dim, dim>;
    using Index   = std::uint64_t;

    std::vector<Element> elements_{};
    std::vector<Point> vertices_{};
    std::vector<Index> indices_{};

  private:
    std::vector<Boundary> outerBoundaries_;
    std::optional<std::vector<std::vector<Boundary>>> innerBoundaries_;
    TransformToSpan<dim> transformer;

    /// Parameters
    static constexpr int maxOuterBoundaryDivisions{10};
    static constexpr int maxInnerBoundaryDivisions{10};

    static constexpr double outerTargetTolerance{1e-4};
    static constexpr double innerTargetTolerance{1e-4};

  public:
    /// brief: Constructs an trimmed_ elementRepresentation with outer and inner boundaries
    explicit TrimmedSubGrid(Trim::ElementBoundaries& _boundaries,
                            std::pair<std::array<double, dim>, std::array<double, dim>> scalingAndOffset)
        : outerBoundaries_(_boundaries.outerBoundaries),
          innerBoundaries_(_boundaries.innerBoundaries),
          transformer(scalingAndOffset) {
      assert(not outerBoundaries_.empty());

      checkAndDivideSmallLoops();
      constructSubGrid();
    }

  private:
    void constructSubGrid() {
      auto boundaries               = splitBoundaries();
      std::tie(indices_, vertices_) = triangulate<dim>(boundaries, transformer);

      for (auto it = indices_.begin(); it < indices_.end(); it += 3)
        elements_.emplace_back(
            Dune::GeometryTypes::triangle,
            std::vector<FieldVector<double, dim>>{vertices_[*it], vertices_[*(it + 1)], vertices_[*(it + 2)]});
    }

    /// @brief This is necessary for loops with only one or 2 boundaries (e.g. circles)
    void checkAndDivideSmallLoops() {
      auto partsForSize = [](std::size_t size) -> int { return size == 1 ? 4 : 2; };

      auto divideLoop = [](std::vector<Boundary>& loop, int parts) {
        std::vector<Boundary> newLoop;
        for (auto& boundary : loop) {
          auto geometry = boundary.nurbsGeometry;
          auto uVec     = Utilities::linspace(boundary.domain, parts + 1);

          for (auto i : std::views::iota(0u, static_cast<unsigned int>(parts)))
            newLoop.emplace_back(geometry, Utilities::Domain<double>{uVec[i], uVec[i + 1]});
        }
        loop = newLoop;
      };
      if (outerBoundaries_.size() < 3)
        divideLoop(outerBoundaries_, partsForSize(outerBoundaries_.size()));

      if (innerBoundaries_)
        for (auto& innerLoop : innerBoundaries_.value())
          if (innerLoop.size() < 3)
            divideLoop(innerLoop, partsForSize(innerLoop.size()));
    }

  public:
    /// @brief Calculates the area from the simplex elements in the current grid
    [[nodiscard]] double calculateArea() const {
      return std::accumulate(elements_.begin(), elements_.end(), 0.0,
                             [](double rhs, const auto& element) { return rhs + element.volume(); });
    }

  private:
    auto splitBoundaries() const {
      return std::make_pair(splitOuterBoundaries(), splitInnerBoundaryLoops());
    }
    auto splitOuterBoundaries(int maxDivisions       = maxOuterBoundaryDivisions,
                              double targetTolerance = outerTargetTolerance) const {
      return splitBoundariesImpl(outerBoundaries_, maxDivisions, targetTolerance, transformer);
    }

    auto splitInnerBoundaryLoops(int maxDivisions       = maxInnerBoundaryDivisions,
                                 double targetTolerance = innerTargetTolerance) const -> decltype(innerBoundaries_) {
      if (innerBoundaries_) {
        typename decltype(innerBoundaries_)::value_type newInnerBoundaries;
        std::ranges::transform(innerBoundaries_.value(), std::back_inserter(newInnerBoundaries), [&](auto& loop) {
          return splitBoundariesImpl(loop, maxDivisions, targetTolerance, transformer);
        });
        return std::make_optional<std::vector<std::vector<Boundary>>>(newInnerBoundaries);
      } else
        return std::nullopt;
    }

  public:
    auto createRefinedGrid(int refinementSteps) const
        -> std::tuple<decltype(elements_), decltype(vertices_), decltype(indices_)> {
      if (refinementSteps == 0)
        return {elements_, vertices_, indices_};

      decltype(elements_) refElements{};
      decltype(vertices_) refVertices{};
      decltype(indices_) refIndices{};

      int splitSteps  = std::ceil((double)refinementSteps / 2.0);
      int refineSteps = refinementSteps - splitSteps;

      auto boundaries = std::make_pair(splitOuterBoundaries(splitSteps, 0.0), splitInnerBoundaryLoops(splitSteps, 0.0));
      auto [ind, vert] = triangulate<dim>(boundaries, transformer);

      Dune::GridFactory<Dune::UGGrid<dim>> gridFactory;
      for (auto& vertex : vert)
        gridFactory.insertVertex(vertex);

      // The element indices are stored as a flat vector, 3 indices always make 1 triangle
      for (auto it = ind.begin(); it < ind.end(); it += 3)
        gridFactory.insertElement(Dune::GeometryTypes::triangle, std::vector<unsigned int>(it, std::next(it, 3)));

      // Attention: This may break down for other grid implementations, as the insertion Index doesn't have to be the
      // "real" entity index

      auto determineOuterBoundaryIndices = [nBoundaries = (unsigned int)boundaries.first.size()](unsigned int bIndex) {
        return std::vector<unsigned int>{bIndex, (bIndex + 1) % nBoundaries};
      };
      auto determineInnerBoundaryIndices = [](unsigned int bIndex, unsigned int nBoundariesInLoop,
                                              unsigned int innerIndex) {
        auto increment = ((innerIndex + 1) % nBoundariesInLoop) - innerIndex;
        return std::vector<unsigned int>{bIndex, bIndex + increment};
      };

      unsigned int bIndex = 0;
      for (auto& boundary : boundaries.first)
        gridFactory.insertBoundarySegment(determineOuterBoundaryIndices(bIndex++),
                                          std::make_shared<GridBoundarySegment<dim>>(boundary, transformer));

      if (boundaries.second)
        for (auto& innerLoop : boundaries.second.value()) {
          unsigned int innerIndex = 0;
          for (auto& boundary : innerLoop)
            gridFactory.insertBoundarySegment(determineInnerBoundaryIndices(bIndex++, innerLoop.size(), innerIndex++),
                                              std::make_shared<GridBoundarySegment<dim>>(boundary, transformer));
        }

      // Create Grid
      auto grid = gridFactory.createGrid();
      grid->globalRefine(refineSteps);

      auto gridView  = grid->leafGridView();
      auto& indexSet = gridView.indexSet();

      // Collect vertices and indices and elements
      auto eleInd = std::vector<unsigned int>();

      for (auto& v : vertices(gridView))
        refVertices.push_back(v.geometry().center());

      for (auto& ele : elements(gridView)) {
        eleInd.clear();
        for (auto i : std::views::iota(0u, ele.subEntities(dim)))
          eleInd.push_back(indexSet.template index<dim>(ele.template subEntity<dim>(i)));

        refElements.emplace_back(Dune::GeometryTypes::triangle,
                                 std::vector<FieldVector<double, dim>>{refVertices[eleInd[0]], refVertices[eleInd[1]],
                                                                       refVertices[eleInd[2]]});
        std::ranges::copy(eleInd, std::back_inserter(refIndices));
      }

      return {refElements, refVertices, refIndices};
    }
  };
} // namespace Dune::IGA

// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "nurbstrimboundary.hh"
#include "nurbstrimmer.hh"
#include "subgridhelpers.hh"

#include <clipper2/clipper.core.h>
#include <mapbox/earcut.hpp>

#include "dune/iga/geometry/geohelper.hh"
#include <dune/grid/uggrid.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/virtualrefinement.hh>



namespace Dune::IGA {

  /** \brief representation of the trimmed_ element in the parameter space */
  template <int dim>
  class TrimmedElementRepresentation {
   public:
    using Point    = Dune::FieldVector<double, dim>;
    using Element = MultiLinearGeometry<double, dim, dim>;
    using Index = std::uint64_t;

    std::vector<Element> elements_{};
    std::vector<Point> vertices_{};
    std::vector<Index> indices_{};


   private:
    std::vector<Boundary> outerBoundaries_;
    std::optional<std::vector<std::vector<Boundary>>> innerBoundaries_;
    std::array<double, dim> scaling_{};
    std::array<double, dim> offset_{};

    /// Parameters
    static constexpr int maxPreSamplesOuterBoundaries{4};
    int innerLoopPreSample{3};
    double targetTolerance{1e-3};

    double targetArea{};

   public:
    /// brief: Constructs an trimmed_ elementRepresentation with outer and inner boundaries
    explicit TrimmedElementRepresentation(Trim::ElementBoundaries& _boundaries,
                                          std::pair<std::array<double, dim>, std::array<double, dim>> scalingAndOffset)
        : outerBoundaries_(_boundaries.outerBoundaries),
          innerBoundaries_(_boundaries.innerBoundaries),
          scaling_{scalingAndOffset.first},
          offset_{scalingAndOffset.second},
          targetArea{calculateTargetArea(200)} {
      assert(not outerBoundaries_.empty());
      reconstructTrimmedElement();
    }

   private:
    // FIXME Helper function
    template <typename VecType>
    [[nodiscard]] auto toLocal(const VecType& cp) const -> VecType {
      VecType local;
      for (int i = 0; i < dim; ++i) {
        local[i] = (cp[i] - offset_[i]) / scaling_[i];
        local[i] = std::clamp(local[i], 0.0, 1.0);
      }
      return local;
    }

    // FIXME Helper function
    [[nodiscard]] auto toLocal(const Clipper2Lib::PointD& p) const -> Clipper2Lib::PointD {
      std::array<double, 2> cp{p.x, p.y};
      auto cpLoc = toLocal(cp);
      return {cpLoc[0], cpLoc[1]};
    }

    // FIXME explanation TODO
    /// \brief
    void reconstructTrimmedElement() {
      if (innerBoundaries_) prepareInnerBoundaries();

      auto boundaries = splitBoundaries();
      std::tie(indices_, vertices_) = triangulate(boundaries);

      for (auto it = indices_.begin(); it < indices_.end(); it += 3)
        elements_.emplace_back(Dune::GeometryTypes::triangle, std::vector<FieldVector<double, dim>>{vertices_[*it], vertices_[*(it + 1)], vertices_[*(it + 2)]});
    }

    [[nodiscard]] auto triangulate(auto& boundaries) const  -> std::pair<std::vector<Index>, std::vector<Point>>{
      // Construct mesh with Earcut
      // C.f. https://github.com/mapbox/earcut.hpp

      auto [splitOuter, splitInner] = boundaries;

      // Create array of points
      auto vertices = std::vector<Point>();

      vertices.reserve(splitOuter.size());
      for (auto& boundary : splitOuter)
        vertices.push_back(toLocal(boundary.endPoints.front()));

      std::vector<std::vector<Point>> polygonInput;
      polygonInput.push_back(vertices);

      if (splitInner)
        for (auto& innerLoop : splitInner.value()) {
          assert(innerLoop.size() > 1);
          std::vector<Point> holeInput;
          for (auto& boundary : innerLoop) {
            auto v = toLocal(boundary.endPoints.front());
            vertices.push_back(v);
            holeInput.push_back(v);
          }
          polygonInput.push_back(holeInput);
        }

      return {mapbox::earcut<Index>(polygonInput), vertices};
    }

    // FIXME better name preparing for what?
    /// \brief
    void prepareInnerBoundaries() {
      assert(innerBoundaries_.has_value());
      for (auto& innerLoop : innerBoundaries_.value()) {
        if (innerLoop.size() == 1) {
          // Divide into 3 parts
          auto boundary = innerLoop.front();
          assert(Dune::FloatCmp::eq(boundary.endPoints.front(), boundary.endPoints.back()));
          auto geometry = boundary.nurbsGeometry;

          innerLoop.clear();
          auto u = Utilities::linspace(boundary.domain, 4);

          innerLoop.emplace_back(geometry, Utilities::Domain<double>{u[0], u[1]});
          innerLoop.emplace_back(geometry, Utilities::Domain<double>{u[1], u[2]});
          innerLoop.emplace_back(geometry, Utilities::Domain<double>{u[2], u[3]});
        }
      }
    }



   public:
    /// \brief Calculates the area from the actual trim paths in [0, 1] domain
    [[nodiscard]] double calculateTargetArea(unsigned int div = 200) const {
      Clipper2Lib::PathD polygon;
      polygon.reserve(div * outerBoundaries_.size());

      for (auto& boundary : outerBoundaries_) {
        auto path = boundary.path(div);
        for (Clipper2Lib::PointD point : path)
          polygon.emplace_back(toLocal(point));
      }
      // ClipperLib returns signed area
      return std::fabs(Clipper2Lib::Area(polygon));
    }

    /// \brief Calculates the area from the simplex elements in the current grid
    [[nodiscard]] double calculateArea() const {
      return std::accumulate(elements_.begin(), elements_.end(), 0.0,
                             [](double rhs, const auto& element) { return rhs + element.volume(); });
    }

   private:
    static constexpr double tolerance = double(16) * std::numeric_limits<double>::epsilon();
    static bool approxSamePoint(const Point& a, const Point& b) { return Dune::FloatCmp::eq(a, b, tolerance); };

    auto splitBoundaries(int maxSplit = maxPreSamplesOuterBoundaries) const {
      return std::make_pair(splitOuterBoundaries(maxSplit), splitInnerBoundaryLoops());
    }
    auto splitOuterBoundaries(int maxSplit) const  {
      return splitBoundariesImpl(outerBoundaries_, maxSplit);
    }

    auto splitInnerBoundaryLoops() const -> decltype(innerBoundaries_) {
      if (innerBoundaries_) {
        typename decltype(innerBoundaries_)::value_type newInnerBoundaries;
        std::ranges::transform(innerBoundaries_.value(), std::back_inserter(newInnerBoundaries), [&](auto& boundaries) {
          return splitBoundariesImpl<false>(boundaries, innerLoopPreSample); });
        return std::make_optional<std::vector<std::vector<Boundary>>>(newInnerBoundaries);
      } else
        return std::nullopt;
    }

    using DomainType = Utilities::Domain<double>;
    struct DomainInformation {
      DomainInformation(const DomainType& d, int i) : domain{d}, localIndex{i} {}
      DomainType domain{};
      int localIndex{};
    };

    template <bool checkArea = true>
    auto splitBoundariesImpl(auto& boundaries, int subSample) const {
      auto refineMap = determineCurvedBoundaries(boundaries);

      // Get all domains
      std::vector<DomainInformation> domains;
      domains.reserve((1 << subSample) * boundaries.size());

      for (auto i : std::views::iota(0u, boundaries.size()))
        domains.emplace_back(boundaries[i].domain, i);

      std::vector<DomainInformation> tempDomainInfos;
      tempDomainInfos.reserve((1 << subSample) * boundaries.size());

      for (auto _ : std::views::iota(0, subSample)) {
        tempDomainInfos.clear();
        tempDomainInfos.insert(tempDomainInfos.end(), domains.begin(), domains.end());
        domains.clear();
        for (auto& domainInfo : tempDomainInfos) {
          bool refine = refineMap[domainInfo.localIndex];

          if (refine) {
            std::array<DomainType, 2> newDomains = Utilities::splitDomainInHalf(domainInfo.domain);
            domains.emplace_back(newDomains[0], domainInfo.localIndex);
            domains.emplace_back(newDomains[1], domainInfo.localIndex);
          } else
            domains.push_back(domainInfo);
        }
        // Check convergence
        if constexpr (checkArea) {
          auto areaFromDomainInformations = [&](const std::vector<DomainInformation>& domainInfos) {
            Clipper2Lib::PathD path;
            for (const auto& domainInfo : domainInfos) {
              auto u = toLocal(outerBoundaries_[domainInfo.localIndex].nurbsGeometry(domainInfo.domain.left()));
              path.emplace_back(u[0], u[1]);
            }
            return std::fabs(Clipper2Lib::Area(path));
          };
          auto actArea = areaFromDomainInformations(domains);
          if (std::fabs(actArea - targetArea) < targetTolerance) break;
        }
      }

      // Create split outerBoundaries from DomainInformation
      std::vector<Boundary> newBoundaries;
      std::ranges::transform(domains, std::back_inserter(newBoundaries), [&boundaries](const auto& domain_info) -> Boundary {
        return {boundaries[domain_info.localIndex].nurbsGeometry, domain_info.domain};
      });
      return newBoundaries;
    }


    [[nodiscard]] auto getVerticesIndices(std::vector<Point>& vertices, Boundary& boundary) const
        -> std::vector<unsigned int> {
      auto it1 = std::ranges::find_if(
          vertices, [&](const Point& v) { return approxSamePoint(v, toLocal(boundary.endPoints.front())); });

      auto it2 = std::ranges::find_if(
          vertices, [&](const Point& v) { return approxSamePoint(v, toLocal(boundary.endPoints.back())); });

      assert(it1 != vertices.end());
      assert(it2 != vertices.end());

      const unsigned int idx1 = std::distance(vertices.begin(), it1);
      const unsigned int idx2 = std::distance(vertices.begin(), it2);

      assert(idx1 != idx2);
      return {idx1, idx2};
    }

    [[nodiscard]] static auto determineCurvedBoundaries(const std::vector<Boundary>& boundaries) -> std::vector<bool> {
      std::vector<bool> result;
      result.reserve(boundaries.size());
      std::ranges::transform(boundaries, std::back_inserter(result),
                             [](auto& boundary) { return boundary.degree() > 1; });

      return result;
    }
   public:
    auto createRefinedGrid(int refinementSteps) const -> std::tuple<decltype(elements_), decltype(vertices_), decltype(indices_)> {

      if (refinementSteps == 0)
        return {elements_, vertices_, indices_};

      decltype(elements_) refElements{};
      decltype(vertices_) refVertices{};
      decltype(indices_) refIndices{};


      int splitSteps = std::ceil((double) refinementSteps / 2.0);
      int refineSteps = refinementSteps - splitSteps;

      auto boundaries = splitBoundaries(splitSteps);
      auto [ind, vert] = triangulate(boundaries);

      Dune::GridFactory<Dune::UGGrid<dim>> gridFactory;
      for (auto& vertex : vert)
        gridFactory.insertVertex(vertex);

      // The element indices are stored as a flat vector, 3 indices always make 1 triangle
      for (auto it = ind.begin(); it < ind.end(); it += 3)
        gridFactory.insertElement(Dune::GeometryTypes::triangle, std::vector<unsigned int>(it, std::next(it, 3)));

      auto toLocalLambda = [&](const auto& cp) { return toLocal(cp); };
      // TODO Without search
      for (auto& boundary : boundaries.first)
        gridFactory.insertBoundarySegment(getVerticesIndices(vert, boundary),
                                          std::make_shared<GridBoundarySegment<dim>>(boundary, toLocalLambda));

      if (boundaries.second)
        for (auto& innerLoop : boundaries.second.value())
          for (auto& boundary : innerLoop)
          gridFactory.insertBoundarySegment(getVerticesIndices(vert, boundary),
                                            std::make_shared<GridBoundarySegment<dim>>(boundary, toLocalLambda));

      // Create Grid
      auto grid = gridFactory.createGrid();
      grid->globalRefine(refineSteps);

      auto gridView = grid->leafGridView();
      auto& indexSet = gridView.indexSet();

      // Collect vertices and indices and elements
      auto eleInd = std::vector<unsigned int>();

      for (auto& v : vertices(gridView))
        refVertices.push_back(v.geometry().center());

      for (auto& ele : elements(gridView)) {
        eleInd.clear();
        for (auto i : std::views::iota(0u, ele.subEntities(dim)))
          eleInd.push_back(indexSet.template index<dim>(ele.template subEntity<dim>(i)));

        refElements.emplace_back(Dune::GeometryTypes::triangle, std::vector<FieldVector<double, dim>>{refVertices[eleInd[0]], refVertices[eleInd[1]],
                                                                                                      refVertices[eleInd[2]]});
        std::ranges::copy(eleInd, std::back_inserter(refIndices));
      }

      return {refElements, refVertices, refIndices};
    }

  };
}  // namespace Dune::IGA

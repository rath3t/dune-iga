// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "nurbstrimboundary.hh"
#include "nurbstrimmer.hh"

#include <clipper2/clipper.core.h>
#include <mapbox/earcut.hpp>

#include "dune/iga/geometry/geohelper.hh"
#include <dune/alugrid/grid.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/virtualrefinement.hh>

// Add support for Dune::FieldVector in Earcut
namespace mapbox::util {

  template <typename T>
  struct nth<0, Dune::FieldVector<T, 2>> {
    inline static auto get(const Dune::FieldVector<T, 2>& t) { return t[0]; };
  };

  template <typename T>
  struct nth<1, Dune::FieldVector<T, 2>> {
    inline static auto get(const Dune::FieldVector<T, 2>& t) { return t[1]; };
  };
}  // namespace mapbox::util

namespace Dune::IGA {

  template <int dim>
  struct GridBoundarySegment : Dune::BoundarySegment<dim, dim, double> {
    explicit GridBoundarySegment(Boundary& _boundary, auto _trFct) : boundary(_boundary), transformFct(_trFct) {}

    Dune::FieldVector<double, dim> operator()(const Dune::FieldVector<double, 1>& localI) const override {
      // u has to be mapped on the domain of 0 to 1
      const auto local = std::clamp(localI[0], 0.0, 1.0);
      double u         = Utilities::mapToRange(local, Utilities::Domain<double>{}, boundary.domain);
      return transformFct(boundary.nurbsGeometry(u));
    }
    std::function<Dune::FieldVector<double, dim>(const Dune::FieldVector<double, dim>&)> transformFct;
    Boundary boundary;
  };

  /** \brief representation of the trimmed element in the parameter space */
  template <int dim>
  class TrimmedElementRepresentation {
   public:
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming, Dune::ALUGridNoComm>;
    using Point    = Dune::FieldVector<double, dim>;

    using Element = MultiLinearGeometry<double, dim, dim>;
    std::vector<Element> subElements{};
    std::vector<Point> subVertices{};
    std::vector<unsigned int> indices{};

   private:

    std::vector<Boundary> outerBoundaries;
    std::optional<std::vector<std::vector<Boundary>>> innerBoundaries;
    bool verbose{};
    bool trimmed{};
    std::array<double, dim> scaling{};
    std::array<double, dim> offset{};

    /// Parameters
    int maxPreSamplesOuterBoundaries{4};
    int innerLoopPreSample{3};
    int preGlobalRefine{0};
    int edgeRefinements{0};
    double targetTolerance{1e-3};
    double targetArea{};

   public:
    /// brief: Constructs an trimmed elementRepresentation with outer and inner boundaries
    explicit TrimmedElementRepresentation(Trim::ElementBoundaries& _boundaries,
                                          std::pair<std::array<double, dim>, std::array<double, dim>> scalingAndOffset)
        : outerBoundaries(_boundaries.outerBoundaries),
          innerBoundaries(_boundaries.innerBoundaries),
          trimmed(true),
          scaling{scalingAndOffset.first},
          offset{scalingAndOffset.second},
          targetArea{calculateTargetArea(200)} {
      reconstructTrimmedElement();
    }
    /// brief: Constructs an untrimmed elementRepresentation gets lazily constructed when called upon the GridView
    explicit TrimmedElementRepresentation(std::pair<std::array<double, dim>, std::array<double, dim>> scalingAndOffset)
        : trimmed(false), scaling{scalingAndOffset.first}, offset{scalingAndOffset.second} {
      constructFromCorners();
    }

    [[nodiscard]] bool isTrimmed() const { return trimmed; }

    void refineAndConstructGrid(unsigned int refinement = 0) {
      if (refinement == 0)
        return;
      if (not trimmed) {
        refineUntrimmed(refinement);
      }
    }

    auto geometryType() {
      return subElements.front().type();
    }

    std::size_t vertexSubIndex(std::uint64_t eleIdx, std::size_t subIndex) {
      if (trimmed)
        return indices[eleIdx * 3 + subIndex];
      else
        return indices[eleIdx * 4 + subIndex];
    }

   private:
    // FIXME Helper function
    template <typename VecType>
    [[nodiscard]] auto toLocal(const VecType& cp) const -> VecType {
      VecType local;
      for (int i = 0; i < dim; ++i) {
        local[i] = (cp[i] - offset[i]) / scaling[i];
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
      if (innerBoundaries) prepareInnerBoundaries();

      splitBoundaries();
      triangulate();

      for (auto it = indices.begin(); it < indices.end(); it += 3)
        subElements.emplace_back(Dune::GeometryTypes::triangle, std::vector<FieldVector<double, dim>>{subVertices[*it], subVertices[*(it + 1)], subVertices[*(it + 2)]});


      /*
      Dune::GridFactory<Grid> gridFactory;
      for (auto& vertex : vertices)
        gridFactory.insertVertex(vertex);

      // The element indices are stored as a flat vector, 3 indices always make 1 triangle
      for (auto it = indices.begin(); it < indices.end(); it += 3)
        gridFactory.insertElement(Dune::GeometryTypes::triangle, std::vector<unsigned int>(it, std::next(it, 3)));

      auto toLocalLambda = [&](const auto& cp) { return toLocal(cp); };
      // TODO Without search
      for (auto& boundary : outerBoundaries)
        gridFactory.insertBoundarySegment(getVerticesIndices(vertices, boundary),
                                          std::make_shared<GridBoundarySegment<dim>>(boundary, toLocalLambda));

      if (innerBoundaries)
        for (auto& innerLoop : innerBoundaries.value())
          for (auto& boundary : innerLoop)
            gridFactory.insertBoundarySegment(getVerticesIndices(vertices, boundary),
                                              std::make_shared<GridBoundarySegment<dim>>(boundary, toLocalLambda));

      // Create Grid
      grid = gridFactory.createGrid();

      if (preGlobalRefine > 0) grid->globalRefine(preGlobalRefine);

      refineGridOnEdges(edgeRefinements);

      GridView gridViewRefined = grid->leafGridView();

      if (verbose)
        std::cout << "Reconstructed Grid with " << gridViewRefined.size(0)
                  << " elements. Area of elements: " << calculateArea()
                  << ". Approx area of trimmed element: " << calculateTargetArea() << std::endl;
      */
    }

    [[nodiscard]] auto triangulate() {
      // Construct mesh with Earcut
      // C.f. https://github.com/mapbox/earcut.hpp

      // Create array of points
      subVertices.reserve(outerBoundaries.size());
      for (auto& boundary : outerBoundaries)
        subVertices.push_back(toLocal(boundary.endPoints.front()));

      std::vector<std::vector<Point>> polygonInput;
      polygonInput.push_back(subVertices);

      if (innerBoundaries)
        for (auto& innerLoop : innerBoundaries.value()) {
          assert(innerLoop.size() > 1);
          std::vector<Point> holeInput;
          for (auto& boundary : innerLoop) {
            auto v = toLocal(boundary.endPoints.front());
            subVertices.push_back(v);
            holeInput.push_back(v);
          }
          polygonInput.push_back(holeInput);
        }

      indices = mapbox::earcut<unsigned int>(polygonInput);
    }

    // FIXME better name preparing for what?
    /// \brief
    void prepareInnerBoundaries() {
      assert(innerBoundaries.has_value());
      for (auto& innerLoop : innerBoundaries.value()) {
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

    void constructFromCorners() {
      subVertices = std::vector<FieldVector<double, dim>>{{0, 0}, {1, 0}, {0, 1}, {1, 1}};
      indices = {0, 1, 2, 3};
      subElements.emplace_back(Dune::GeometryTypes::quadrilateral, subVertices);
    }

   public:
    /// \brief Calculates the area from the actual trim paths in [0, 1] domain
    [[nodiscard]] double calculateTargetArea(unsigned int div = 200) const {
      if (not isTrimmed()) throw std::runtime_error("calculateTargetArea only defined for trimmed elements");

      Clipper2Lib::PathD polygon;
      polygon.reserve(div * outerBoundaries.size());

      for (auto& boundary : outerBoundaries) {
        auto path = boundary.path(div);
        for (Clipper2Lib::PointD point : path)
          polygon.emplace_back(toLocal(point));
      }
      // ClipperLib returns signed area
      return std::fabs(Clipper2Lib::Area(polygon));
    }

    /// \brief Calculates the area from the simplex elements in the current grid
    [[nodiscard]] double calculateArea() const {
      return std::accumulate(subElements.begin(), subElements.end(), 0.0,
                             [](double rhs, const auto& element) { return rhs + element.volume(); });
    }

   private:
    static constexpr double tolerance = double(16) * std::numeric_limits<double>::epsilon();
    static bool approxSamePoint(const Point& a, const Point& b) { return Dune::FloatCmp::eq(a, b, tolerance); };

    void splitBoundaries() {
      splitOuterBoundaries();
      splitInnerBoundaryLoops();
    }

    using DomainType = Utilities::Domain<double>;
    struct DomainInformation {
      DomainInformation(const DomainType& d, int i) : domain{d}, localIndex{i} {}
      DomainType domain{};
      int localIndex{};
    };

    template <bool checkArea = true>
    void splitBoundariesImpl(auto& boundaries, int subSample) {
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
              auto u = toLocal(outerBoundaries[domainInfo.localIndex].nurbsGeometry(domainInfo.domain.left()));
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
      std::ranges::transform(domains, std::back_inserter(newBoundaries), [&](const auto& domain_info) -> Boundary {
        return {boundaries[domain_info.localIndex].nurbsGeometry, domain_info.domain};
      });
      boundaries = std::move(newBoundaries);
    }
    void splitOuterBoundaries() { splitBoundariesImpl(outerBoundaries, maxPreSamplesOuterBoundaries); }

    void splitInnerBoundaryLoops() {
      if (innerBoundaries)
        std::ranges::for_each(innerBoundaries.value(),
                              [&](auto& boundaries) { splitBoundariesImpl<false>(boundaries, innerLoopPreSample); });
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

    void refineUntrimmed(int refinementSteps) {
        assert(geometryType() == Dune::GeometryTypes::quadrilateral && not trimmed);

        subVertices.clear();
        subElements.clear();
        indices.clear();

        Dune::RefinementIntervals tag(refinementSteps + 1);
        Dune::VirtualRefinement<dim, double>& refinement = Dune::buildRefinement<dim, double>(Dune::GeometryTypes::quadrilateral, Dune::GeometryTypes::quadrilateral);

        auto eSubEnd = refinement.eEnd(tag);
        auto eSubIt = refinement.eBegin(tag);

        auto vSubEnd = refinement.vEnd(tag);
        auto vSubIt = refinement.vBegin(tag);

        subElements.reserve(refinement.nElements(tag));
        subVertices.reserve(refinement.nVertices(tag));
        indices.reserve(refinement.nElements(tag) * 3);

        for (; vSubIt != vSubEnd; ++vSubIt)
          subVertices.push_back(vSubIt.coords());

        std::vector<Point> eleCoords;
        eleCoords.reserve(4);

        for (; eSubIt != eSubEnd; ++eSubIt) {
          eleCoords.clear();
          std::ranges::copy(eSubIt.vertexIndices(), std::back_inserter(indices));

          for (auto idx : eSubIt.vertexIndices())
            eleCoords.push_back(subVertices[idx]);

          subElements.emplace_back(Dune::GeometryTypes::quadrilateral, eleCoords);
        }


    }

  };
}  // namespace Dune::IGA

// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/virtualrefinement.hh>
#include <dune/iga/trimmer/defaulttrimmer/integrationrules/simplexintegrationrulegenerator.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

namespace Dune::IGA {

template <class GridView>
class IGARefinedGeometries
{
public:
  static constexpr int dim = GridView::dimension;

  using Element = MultiLinearGeometry<double, dim, dim>;
  using Point   = Dune::FieldVector<double, dim>;
  using Index   = std::uint64_t;
  using IDType  = typename GridView::Grid::GlobalIdSet::IdType;

private:
  struct ElementData
  {
    std::vector<Element> elements{};
    std::vector<Point> vertices{};
    std::vector<Index> indices{};
  };

  std::unordered_map<IDType, ElementData> trimmedElementData_;
  ElementData cubeData{};

public:
  IGARefinedGeometries(const GridView& gridView, const int subSampleFull, const int subSampleTrimmed) {
    assert(subSampleFull >= 0 and subSampleTrimmed >= 0 && "subSamples have to be zero or positive");

    using Trimmer = typename GridView::GridViewImp::TrimmerType;
    static_assert(
        std::is_same_v<Trimmer, IGANEW::DefaultTrim::TrimmerImpl<Trimmer::mydimension, Trimmer::dimensionworld,
                                                                 typename Trimmer::ctype>>);

    createCubeRefinement(subSampleFull);

    const auto& idSet = gridView.grid().globalIdSet();

    for (const auto& element : elements(gridView)) {
      if (element.impl().isTrimmed()) {
        auto [ele, vert, ind] =
            IGANEW::DefaultTrim::SimplexIntegrationRuleGenerator<typename GridView::Grid>::createSimplicies(element);
        trimmedElementData_.emplace(idSet.id(element), ElementData{ele, vert, ind});
      }
    }
  }

  [[nodiscard]] const std::vector<Element>& getElements(IDType eIndex) const {
    if (isTrimmed(eIndex))
      return trimmedElementData_.at(eIndex).elements;
    else
      return cubeData.elements;
  }

  [[nodiscard]] const std::vector<Point>& getVertices(IDType eIndex) const {
    if (isTrimmed(eIndex))
      return trimmedElementData_.at(eIndex).vertices;
    else
      return cubeData.vertices;
  }

  [[nodiscard]] const std::vector<Index>& getIndices(IDType eIndex) const {
    if (isTrimmed(eIndex))
      return trimmedElementData_.at(eIndex).indices;
    else
      return cubeData.indices;
  }

  [[nodiscard]] Index vertexSubIndex(IDType eIndex, Index subEleIndex, Index subEntityIndex) const {
    // As in a trimmed element only triangles are present, the index offset is 3, for untrimmed elements, the subgrid
    // + is made up of quadrilaterals, therefore the offset is 4
    int offset = isTrimmed(eIndex) ? 3 : 4;
    if (isTrimmed(eIndex))
      return trimmedElementData_.at(eIndex).indices[subEleIndex * offset + subEntityIndex];
    else
      return cubeData.indices.at(subEleIndex * offset + subEntityIndex);
  }

  [[nodiscard]] std::size_t nElements(IDType eIndex) const {
    if (isTrimmed(eIndex))
      return trimmedElementData_.at(eIndex).elements.size();
    else
      return cubeData.elements.size();
  }

  [[nodiscard]] std::size_t nVertices(IDType eIndex) const {
    if (isTrimmed(eIndex))
      return trimmedElementData_.at(eIndex).vertices.size();
    else
      return cubeData.vertices.size();
  }

  [[nodiscard]] GeometryType geometryType(IDType eIndex) const {
    if (isTrimmed(eIndex))
      return GeometryTypes::simplex(dim);
    else
      return GeometryTypes::cube(dim);
  }

private:
  void createCubeRefinement(const int subSample) {
    Dune::RefinementIntervals tag{subSample + 1};
    Dune::VirtualRefinement<dim, double>& refinement =
        Dune::buildRefinement<dim, double>(Dune::GeometryTypes::cube(dim), Dune::GeometryTypes::cube(dim));

    auto eSubEnd = refinement.eEnd(tag);
    auto eSubIt  = refinement.eBegin(tag);

    auto vSubEnd = refinement.vEnd(tag);
    auto vSubIt  = refinement.vBegin(tag);

    cubeData.elements.reserve(refinement.nElements(tag));
    cubeData.vertices.reserve(refinement.nVertices(tag));
    cubeData.indices.reserve(refinement.nElements(tag) * 4);

    for (; vSubIt != vSubEnd; ++vSubIt) {
      cubeData.vertices.push_back(vSubIt.coords());
    }
    std::vector<Point> eleCoords;
    eleCoords.reserve(4);

    for (; eSubIt != eSubEnd; ++eSubIt) {
      eleCoords.clear();
      std::ranges::copy(eSubIt.vertexIndices(), std::back_inserter(cubeData.indices));

      for (auto idx : eSubIt.vertexIndices())
        eleCoords.push_back(cubeData.vertices[idx]);

      cubeData.elements.emplace_back(Dune::GeometryTypes::cube(dim), eleCoords);
    }
  }

  [[nodiscard]] bool isTrimmed(IDType eIndex) const {
    return trimmedElementData_.contains(eIndex);
  }
};

} // namespace Dune::IGA
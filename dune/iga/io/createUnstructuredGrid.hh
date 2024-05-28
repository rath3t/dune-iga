// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/io/vtk/vtkrefinedgeometries.hh>

namespace Dune::IGA {

namespace Impl {
  struct Comparator
  {
    bool operator()(const Dune::FieldVector<double, 2>& p1, const Dune::FieldVector<double, 2>& p2) const {
      if (p1[0] < p2[0])
        return true;
      if (p1[0] > p2[0])
        return false;
      return p1[1] < p2[1];
    };
  };
} // namespace Impl

template <typename UnstructuredGrid, typename PatchGrid>
std::unique_ptr<UnstructuredGrid> createUnstructuredGridImpl(PatchGrid* patchGrid) {
  auto gridView = patchGrid->leafGridView();
  auto& idSet   = patchGrid->globalIdSet();

  IGARefinedGeometries geometries(gridView, 0, 0);

  std::set<FieldVector<double, 2>, Impl::Comparator> vertices;

  for (const auto& element : elements(gridView)) {
    auto id       = idSet.id(element);
    auto geometry = element.geometry();

    std::ranges::transform(geometries.getVertices(id), std::inserter(vertices, vertices.begin()),
                           [&](const auto& v) { return geometry.global(v); });
  }

  auto gridFactory = Dune::GridFactory<UnstructuredGrid>();

  for (const auto& v : vertices) {
    gridFactory.insertVertex(v);
  }

  // Reconstruct grid
  for (const auto& element : elements(gridView)) {
    auto geometry = element.geometry();
    auto id       = idSet.id(element);

    auto gt            = geometries.geometryType(id);
    unsigned int nSubI = gt == GeometryTypes::simplex(2) ? 3 : 4;

    auto& eleVertices = geometries.getVertices(id);

    for (auto subEleIdx : std::views::iota(0u, geometries.nElements(id))) {
      std::vector<unsigned int> elementVertices;

      for (auto subEntityIndex : std::views::iota(0u, nSubI)) {
        auto localVertexIdx                 = geometries.vertexSubIndex(id, subEleIdx, subEntityIndex);
        Dune::FieldVector<double, 2> vertex = geometry.global(eleVertices[localVertexIdx]);

        // Find Idx
        auto it = vertices.find(vertex);
        assert(it != vertices.end());

        elementVertices.push_back(std::distance(vertices.begin(), it));
      }
      gridFactory.insertElement(gt, elementVertices);
    }
  }

  return gridFactory.createGrid();
}
} // namespace Dune::IGA

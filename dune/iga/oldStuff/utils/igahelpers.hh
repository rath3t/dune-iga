// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/bitsetvector.hh>
#include <dune/geometry/referenceelements.hh>

template <class GridView, int ncomp = 1>
class BoundaryPatchEnclosingVerticesPropertyTrimmed
{
  typedef typename GridView::IndexSet IndexSet;
  static const int dim = GridView::dimension;

public:
  typedef typename GridView::Intersection Intersection;

  /** @brief Create property from a marker vector
   */
  BoundaryPatchEnclosingVerticesPropertyTrimmed(const GridView& gridView, const Dune::BitSetVector<ncomp>& vertices)
      : indexSet_(gridView.indexSet()),
        vertices_(vertices) {}

  /** @brief Check if intersection is enclosed by vertices in the vector
   */
  bool operator()(const Intersection& i) const {
    const auto inside  = i.inside();
    int localFaceIndex = i.indexInInside();

    if (inside.impl().isTrimmed())
      return false;

    auto refElement = Dune::ReferenceElements<double, dim>::general(inside.type());

    // Get global node ids
    int n = refElement.size(localFaceIndex, 1, dim);

    // Using ReferenceElement::subEntity is OK here, because we loop
    // over all subEntities (i.e. sub-vertices) and just return false
    // if _any_ of them is not marked.
    for (int ii = 0; ii < n; ii++) {
      int localVertexIndex = refElement.subEntity(localFaceIndex, 1, ii, dim);
      if (not(vertices_[indexSet_.subIndex(inside, localVertexIndex, dim)].any()))
        return false;
    }
    return true;
  }

private:
  const IndexSet& indexSet_;
  const Dune::BitSetVector<ncomp>& vertices_;
};

namespace Dune::Functions {

template <class Basis, class F>
void forEachUntrimmedBoundaryDOF(const Basis& basis, F&& f) {
  auto localView       = basis.localView();
  auto seDOFs          = subEntityDOFs(basis);
  const auto& gridView = basis.gridView();
  for (auto&& element : elements(gridView))
    if (element.hasBoundaryIntersections() and not element.impl().isTrimmed()) {
      localView.bind(element);
      for (const auto& intersection : intersections(gridView, element))
        if (intersection.boundary())
          for (auto localIndex : seDOFs.bind(localView, intersection))
            f(localIndex, localView, intersection);
    }
}
} // namespace Dune::Functions

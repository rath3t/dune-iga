// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#if HAVE_DUNE_FUFEM
  #include "dune/iga/utils/igahelpers.hh"
  #include <dune/fufem/boundarypatch.hh>
  #include <dune/python/pybind11/eigen.h>
  #include <dune/python/pybind11/functional.h>
  #include <dune/python/pybind11/pybind11.h>
  #include <dune/python/pybind11/stl.h>

namespace Dune::IGA::Python {

// Python wrapper for the FVAssembler C++ class
template <class BoundaryPatch, class... options>
void registerBoundaryPatch(pybind11::handle scope, pybind11::class_<BoundaryPatch, options...> cls) {
  using pybind11::operator""_a;

  using GridView = typename BoundaryPatch::GridView;

  cls.def(pybind11::init([](const GridView& gv, Eigen::Ref<Eigen::VectorX<bool>> vec) {
            Dune::BitSetVector<1> bitSetVector;
            bitSetVector.resize(vec.size());
            for (size_t i = 0; i < vec.size(); ++i)
              bitSetVector[i] = vec[i];

            BoundaryPatch* neumannBoundary = new BoundaryPatch(gv);

            // create property for insertion by vertex vector
            BoundaryPatchEnclosingVerticesPropertyTrimmed<GridView, 1> prop(gv, bitSetVector);
            neumannBoundary->insertFacesByProperty(prop);
            return neumannBoundary;
          }),
          pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());

  cls.def("gridView", &BoundaryPatch::gridView);
}

} // namespace Dune::IGA::Python
#endif

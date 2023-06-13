// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include "dune/python/iga/gridenums.hh"
#include <dune/iga/io/ibra/ibrareader.hh>
#include <dune/iga/nurbsgrid.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/functions/globalbasis.hh>
#include <dune/python/grid/capabilities.hh>
#include <dune/python/grid/enums.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/pybind11/pybind11.h>

namespace Dune::Python::IGA {

  template <class NURBSGrid, class... options>
  void registerGlobalBasis(pybind11::module handle, pybind11::class_<NURBSGrid, options...> cls) {
    Dune::Python::registerGlobalBasis(handle, cls);
  }
}  // namespace Dune::Python::IGA

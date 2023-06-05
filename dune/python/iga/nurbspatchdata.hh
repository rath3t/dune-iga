// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/fufem/boundarypatch.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/functions/globalbasis.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include "dune/iga/nurbspatchdata.hh"

namespace Dune::Python {

template <class NURBSPatchData, class... options>
void registerNurbsPatchData(pybind11::handle scope, pybind11::class_<NURBSPatchData, options...> cls) {
  using pybind11::operator""_a;

  using GlobalCoordinateType    = typename NURBSPatchData::GlobalCoordinateType;
   constexpr int dimWorld    =  NURBSPatchData::GlobalCoordinateType::dimension;
  using ControlPointType      = typename NURBSPatchData::ControlPointType;
  using ControlPointNetType      = typename NURBSPatchData::ControlPointNetType;
  using GridView       = typename GlobalBasis::GridView;
  using Element        = typename LinearElastic::Element;
  using Traits         = typename LinearElastic::Traits;
  using FErequirements = typename LinearElastic::FERequirementType;

  cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu) {
            return new LinearElastic(basis, element, emod, nu);
          }),
          pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());

  using LoadFunction = std::function<Eigen::Vector<double, Traits::worlddim>(Eigen::Vector<double, Traits::worlddim>,
                                                                             const double&)>;
  cls.def(pybind11::init(
              [](const GlobalBasis& basis, const Element& element, double emod, double nu,
                 const LoadFunction volumeLoad) { return new LinearElastic(basis, element, emod, nu, volumeLoad); }),
          );

}

}  // namespace Ikarus::Python

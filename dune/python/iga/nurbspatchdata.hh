// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/nurbsbasis.hh"
#include "dune/iga/splines/nurbsalgorithms.hh"
#include <dune/python/functions/globalbasis.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/iga/utils/mdnet.hh>
#include <dune/iga/geometrykernel/controlpoint.hh>


namespace Dune::Python {

template <class ControlPoint, class... options>
void registerControlPoint(pybind11::handle scope, pybind11::class_<ControlPoint, options...> cls) {
  using pybind11::operator""_a;

  using VectorType = typename ControlPoint::VectorType;
  using ScalarType = typename VectorType::value_type;

  cls.def(pybind11::init());
  cls.def(pybind11::self + pybind11::self)
      .def(pybind11::self - pybind11::self)
      .def(pybind11::self += pybind11::self)
      .def(-pybind11::self)
      .def(pybind11::self *= double())
      .def(double() * pybind11::self)
      .def(pybind11::self * double());
  cls.def(pybind11::init([](const VectorType& vec, const ScalarType& weight = 1) {
            return new ControlPoint({.p = vec, .w = weight});
          }),
          "Coordinates"_a, "Weight"_a);
  cls.def_readwrite("coords", &ControlPoint::p).def_readwrite("weight", &ControlPoint::w);
}

template <class MultiDimensionNet, class... options>
void registerMultiDimensionalNet(pybind11::handle scope, pybind11::class_<MultiDimensionNet, options...> cls) {
  using pybind11::operator""_a;

  using ValueType              = typename MultiDimensionNet::value_type;
  constexpr std::size_t netDim = MultiDimensionNet::netDim;

  cls.def(pybind11::init());
  if constexpr (netDim == 2)
    cls.def(pybind11::init(
        [](const std::vector<std::vector<ValueType>>& values) { return new MultiDimensionNet(values); }));
  else if constexpr (netDim == 1)
    cls.def(pybind11::init([](const std::vector<ValueType>& values) { return new MultiDimensionNet(values); }));
  else if constexpr (netDim == 3)
    cls.def(pybind11::init(
        [](const std::vector<std::vector<std::vector<ValueType>>>& values) { return new MultiDimensionNet(values); }));

  cls.def("__len__", &MultiDimensionNet::size);
  cls.def("strideSizes", &MultiDimensionNet::strideSizes);
  cls.def("set", &MultiDimensionNet::set);
  cls.def("directGet", [](MultiDimensionNet& self, int directIndex) { return self.directGet(directIndex); });
  cls.def("get", [](MultiDimensionNet& self, const std::array<int, netDim>& multIndex) { return self.get(multIndex); });
  cls.def_property_readonly_static("netDim", [](pybind11::object /* self */) { return MultiDimensionNet::netDim; });
  cls.def_property_readonly_static("valueType",
                                   [](pybind11::object /* self */) { return Dune::className<ValueType>(); });
}

template <class NURBSPatchData, class... options>
void registerNurbsPatchData(pybind11::handle scope, pybind11::class_<NURBSPatchData, options...> cls) {
  using pybind11::operator""_a;

  using GlobalCoordinateType     = typename NURBSPatchData::GlobalCoordinateType;
  constexpr int dimWorld         = NURBSPatchData::GlobalCoordinateType::dimension;
  using ControlPointType         = typename NURBSPatchData::ControlPointType;
  using ControlPointNetType      = typename NURBSPatchData::ControlPointNetType;
  constexpr std::size_t patchDim = NURBSPatchData::patchDim;
  constexpr std::size_t dimworld = NURBSPatchData::dimworld;

  cls.def(pybind11::init());
  cls.def(pybind11::init([](const std::array<std::vector<double>, patchDim>& knotSpansI,
                            const ControlPointNetType& controlPointsI, const std::array<int, patchDim>& degreeInput) {
    return new NURBSPatchData(knotSpansI, controlPointsI, degreeInput);
  }));
  cls.def_readwrite("knotSpans", &NURBSPatchData::knotSpans)
      .def_readwrite("controlPoints", &NURBSPatchData::controlPoints)
      .def_readwrite("degree", &NURBSPatchData::degree);

  cls.def_property_readonly_static("patchDim", [](pybind11::object /* self */) { return patchDim; });
  cls.def_property_readonly_static("dimworld", [](pybind11::object /* self */) { return dimworld; });

  cls.def("asBasis", [](NURBSPatchData& self) { return Dune::Functions::BasisFactory::nurbs(self); });
  cls.def("degreeElevate", [](NURBSPatchData& self, const int refinementDirection, const int elevationFactor) {
    return Dune::IGA::Splines::degreeElevate(self, refinementDirection, elevationFactor);
  });
  cls.def("knotRefinement",
          [](NURBSPatchData& self, const std::vector<double>& newKnots, const int refinementDirection) {
            return Dune::IGA::Splines::knotRefinement(self, newKnots, refinementDirection);
          });
}

} // namespace Dune::Python

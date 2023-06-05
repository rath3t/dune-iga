// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

//#include <dune/fufem/boundarypatch.hh>
//#include <dune/functions/functionspacebases/lagrangebasis.hh>
//#include <dune/functions/functionspacebases/powerbasis.hh>
//#include <dune/grid/yaspgrid.hh>
//#include <dune/python/common/typeregistry.hh>
#include <dune/python/functions/globalbasis.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include "dune/iga/nurbspatchdata.hh"
#include "dune/iga/nurbsbasis.hh"
#include "dune/iga/controlpoint.hh"

namespace Dune::Python {

template <class ControlPoint, class... options>
void registerControlPoint(pybind11::handle scope, pybind11::class_<ControlPoint, options...> cls) {
  using pybind11::operator""_a;

  using VectorType    = typename ControlPoint::VectorType;
  using ScalarType    = typename VectorType::value_type;

  cls.def(pybind11::init());
  cls.def(pybind11::self + pybind11::self)
      .def(pybind11::self - pybind11::self)
      .def(pybind11::self += pybind11::self)
      .def(-pybind11::self)
      .def(pybind11::self *= double())
      .def(double() * pybind11::self)
      .def(pybind11::self * double());
  cls.def(pybind11::init([](const VectorType& vec, const ScalarType& weight=1) {
            return new ControlPoint({.p=vec,.w=weight});
          }),"Coordinates"_a,"Weight"_a
  );
  cls.def_readwrite("coords", &ControlPoint::p)
  .def_readwrite("weight", &ControlPoint::w);

}


template <class MultiDimensionNet, class... options>
void registerMultiDimensionNet(pybind11::handle scope, pybind11::class_<MultiDimensionNet, options...> cls) {
  using pybind11::operator""_a;

  using ValueType    = typename MultiDimensionNet::value_type;
  constexpr std::size_t netDim =      MultiDimensionNet::netDim;

  cls.def(pybind11::init());
  cls.def(pybind11::init([](const std::vector<std::vector<ValueType>>& values) {
            return new MultiDimensionNet(values);
          })
  );

  cls.def("set", &MultiDimensionNet::set);
  cls.def("get", [](MultiDimensionNet& self,const std::array<int, netDim>& multIndex  ) { return self.get(multIndex); });
  cls.def_property_readonly_static("netDim", [](pybind11::object /* self */) { return MultiDimensionNet::netDim; });
//  cls.def(pybind11::init([](std::array<int, netDim> dimSize, const std::vector<std::vector<ValueType>> values) {
//            return new MultiDimensionNet(dimSize,values);
//          })
//  );
//
//  cls.def(pybind11::init([](std::array<int, netDim> dimSize, const std::vector<std::vector<std::vector<ValueType>>>& values) {
//            return new MultiDimensionNet(dimSize,values);
//          })
//  );

}

template <class NURBSPatchData, class... options>
void registerNurbsPatchData(pybind11::handle scope, pybind11::class_<NURBSPatchData, options...> cls) {
  using pybind11::operator""_a;

  using GlobalCoordinateType    = typename NURBSPatchData::GlobalCoordinateType;
   constexpr int dimWorld    =  NURBSPatchData::GlobalCoordinateType::dimension;
  using ControlPointType      = typename NURBSPatchData::ControlPointType;
  using ControlPointNetType      = typename NURBSPatchData::ControlPointNetType;
  constexpr std::size_t netDim = ControlPointNetType::netdim;

  cls.def(pybind11::init());
  cls.def(pybind11::init([](const std::array<std::vector<double>, netDim>& knotSpansI, const ControlPointNetType& controlPointsI,
                            const std::array<int, netDim>& degreeInput) {
            return new NURBSPatchData(netDim,controlPointsI,degreeInput);
          })
          );
  cls.def_readwrite("knotSpans", &NURBSPatchData::knotSpans)
      .def_readwrite("controlPoints", &NURBSPatchData::controlPoints)
      .def_readwrite("degree", &NURBSPatchData::degree);


  cls.def("asNurbs", [](NURBSPatchData& self) {
    return Dune::Functions::BasisFactory::nurbs(self);
  },pybind11::keep_alive<0, 1>());
}

}  // namespace Ikarus::Python

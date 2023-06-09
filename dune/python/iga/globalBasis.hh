// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include <dune/python/grid/enums.hh>
#include <dune/iga/io/ibra/ibrareader.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/grid/capabilities.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/common/typeregistry.hh>
#include "dune/python/iga/gridenums.hh"
#include <dune/iga/nurbsgrid.hh>
#include <dune/python/functions/globalbasis.hh>




namespace Dune::Python::IGA
{

template <class NURBSGrid, class... options>
void registerIGAGlobalBasis(pybind11::module handle, pybind11::class_<NURBSGrid, options...> cls) {
  using pybind11::operator""_a;


  Dune::Python::registerGlobalBasis(handle,cls);

//static constexpr std::integral auto dimension      = NURBSGrid::dimension;
//static constexpr std::integral auto dimensionworld = NURBSGrid::dimensionworld;
//using ctype                                        = typename NURBSGrid::ctype;

//  module.def( "reader", [] ( const pybind11::dict &args_ ) { return Dune::Python::IGA::reader< NURBSGrid >( args_ ); } );
//
//    Dune::Python::registerHierarchicalGrid (module, cls);
//
//  auto clsLeafView = insertClass< typename NURBSGrid::LeafGridView >( module, "LeafGrid", GenerateTypeName( cls, "LeafGridView" ) );
//  if( clsLeafView.second )
//  registerGridView( module, clsLeafView.first );
//
//  clsLeafView.first.def("preBasis",[](const typename NURBSGrid::LeafGridView& self){return self.impl().preBasis();});
//
//
//using ControlPointNetType    = typename NURBSGrid::ControlPointNetType;
//  using NURBSPatchDataType    = typename NURBSGrid::NURBSPatchDataType;
//
//  cls.def(pybind11::init([](const std::array<std::vector<double>, dimension>& knotSpans, const ControlPointNetType& controlPoints,
//                            const std::array<int, dimension>& order){return new NURBSGrid(knotSpans,controlPoints,order);}));
//
//  cls.def(pybind11::init([](const NURBSPatchDataType& nurbsPatchData){return new NURBSGrid(nurbsPatchData);}));
//  cls.def("globalRefineInDirection",[]( NURBSGrid& self,const int dir, const int refinementLevel, bool omitTrim = false){self.globalRefineInDirection(dir,refinementLevel,omitTrim);});
//  cls.def("patchData",[](const NURBSGrid& self,int i = 0){return self.patchData(i);});




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
}

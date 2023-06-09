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


namespace Dune::Python {
// we have to ahead of https://gitlab.dune-project.org/core/dune-grid/-/blob/releases/2.9/dune/python/grid/hierarchical.hh?ref_type=heads#L233 thus we use requires(IsSpecializationTwoNonTypesAndType<Dune::IGA::NURBSGrid,Grid>::value)
// to be sure this overload is used
//make sure this type is used if an iga grid is passed this is function is enabled and used by adl
template<template<auto, auto, typename> class Type, typename>
struct IsSpecializationTwoNonTypesAndType : std::false_type {};

template<template<auto, auto, typename> class Type, auto T, auto T2, typename S>
struct IsSpecializationTwoNonTypesAndType<Type, Type<T, T2, S>> : std::true_type {};

template<int dim, int dimworld, typename ScalarType>
struct Capabilities::HasGridFactory<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>> : public std::integral_constant<
    bool,
    false> {
};


}
namespace Dune::Python::IGA
{
template<class Grid>
requires(IsSpecializationTwoNonTypesAndType<Dune::IGA::NURBSGrid, Grid>::value)
inline static std::shared_ptr<Grid> reader(const pybind11::dict &dict) {

  std::string file_path;
  Dune::Python::IGA::Reader reader = IGA::Reader::json;
  bool trim = true;
  std::array<int, 2> elevateDegree = {0, 0};
  std::array<int, 2> preKnotRefine = {0, 0};
  std::array<int, 2> postKnotRefine = {0, 0};
  if (dict.contains("reader"))
    reader = dict["reader"].cast<Dune::Python::IGA::Reader>();

  switch (reader) {
    case IGA::Reader::json:

      if (dict.contains("file_path"))
        file_path = dict["file_path"].cast<std::string>();
      else
        DUNE_THROW(Dune::IOError, "No field in dict with name file_path. Unable to read grid");
      if (dict.contains("trim"))
        trim = dict["trim"].cast<bool>();
      if (dict.contains("elevate_degree"))
        elevateDegree = dict["elevate_degree"].cast<std::array<int, 2>>();
      if (dict.contains("pre_knot_refine"))
        preKnotRefine = dict["pre_knot_refine"].cast<std::array<int, 2>>();
      if (dict.contains("post_knot_refine"))
        postKnotRefine = dict["post_knot_refine"].cast<std::array<int, 2>>();

      static constexpr std::integral auto dim = Grid::dimension;
      static constexpr std::integral auto dimworld = Grid::dimensionworld;
      using ScalarType = typename Grid::ctype;
      return Dune::IGA::IbraReader < dim, dimworld, ScalarType
          > ::read(file_path, trim, elevateDegree, preKnotRefine, postKnotRefine);
    default:DUNE_THROW(Dune::NotImplemented, "Your requested reader is not implemeneted");
  }
}

template <class NURBSGrid, class... options>  requires(IsSpecializationTwoNonTypesAndType<Dune::IGA::NURBSGrid,NURBSGrid>::value)
void registerHierarchicalGrid(pybind11::module module, pybind11::class_<NURBSGrid, options...> cls) {
  using pybind11::operator""_a;

static constexpr std::integral auto dimension      = NURBSGrid::dimension;
static constexpr std::integral auto dimensionworld = NURBSGrid::dimensionworld;
using ctype                                        = typename NURBSGrid::ctype;

  module.def( "reader", [] ( const pybind11::dict &args_ ) { return Dune::Python::IGA::reader< NURBSGrid >( args_ ); } );

    Dune::Python::registerHierarchicalGrid (module, cls);

  auto clsLeafView = insertClass< typename NURBSGrid::LeafGridView >( module, "LeafGrid", GenerateTypeName( cls, "LeafGridView" ) );
  if( clsLeafView.second )
  registerGridView( module, clsLeafView.first );

  clsLeafView.first.def("preBasis",[](const typename NURBSGrid::LeafGridView& self){return self.impl().preBasis();});


using ControlPointNetType    = typename NURBSGrid::ControlPointNetType;
  using NURBSPatchDataType    = typename NURBSGrid::NURBSPatchDataType;

  cls.def(pybind11::init([](const std::array<std::vector<double>, dimension>& knotSpans, const ControlPointNetType& controlPoints,
                            const std::array<int, dimension>& order){return new NURBSGrid(knotSpans,controlPoints,order);}));

  cls.def(pybind11::init([](const NURBSPatchDataType& nurbsPatchData){return new NURBSGrid(nurbsPatchData);}));
//  cls.def("hierarchicalGrid",[](const self& nurbsPatchData){return new NURBSGrid(nurbsPatchData);}));




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

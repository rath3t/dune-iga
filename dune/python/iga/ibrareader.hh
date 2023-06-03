// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/python/grid/hierarchical.hh>

#include "dune/iga/io/ibra/ibrareader.hh"
#include "dune/python/grid/range.hh"
#include "dune/iga/nurbsleafgridview.hh"
#include < dune/python/grid/enums.hh>
//#include "modregisterGridView.hh"


namespace Dune::Python {

template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
inline static std::shared_ptr< Dune::IGA::NURBSGrid<dim, dimworld, ScalarType> > reader ( const std::tuple< Reader, std::string > &args )
{
  switch( std::get< 0 >( args ) )
  {
    case Reader::dgf:
      std::cout<<"Reader::dgf"<<std::endl;
      return readDGF< Grid >( std::get< 1 >( args ) );

    case Reader::dgfString:
    {
      std::istringstream input( std::get< 1 >( args ) );
      std::cout<<"Reader::dgfString"<<std::endl;
      return readDGF< Grid >( input );
    }

    case Reader::gmsh:
      return readGmsh< Grid >( std::get< 1 >( args ) );

    default:
      return nullptr;
  }
}



// Python wrapper for the FVAssembler C++ class
template <class IbraReader, class... options>
void registerIbraReader(pybind11::handle scope, pybind11::class_<IbraReader, options...> cls) {
  namespace py = pybind11;
  using pybind11::operator""_a;
  pybind11::module::import( "dune.geometry" );
  pybind11::module::import( "dune.grid" );
  auto module = pybind11::module::import( "dune.iga" );
  using Grid                = typename IbraReader::Grid;
   constexpr  auto dimension      = Grid::dimension;
   constexpr  auto dimensionworld = Grid::dimensionworld;
  using GridView                = typename Grid::GridView;

  cls.def(pybind11::init());

  auto includes = Dune::Python::IncludeFiles{"dune/iga/nurbsgrid.hh"};
  auto clsName = Dune::className<IbraReader>();
  auto clsGrid = insertClass< Grid >( module, "HierarchicalGrid", Dune::Python::GenerateTypeName(clsName),
                                                     includes );
  if( clsGrid.second )
    registerHierarchicalGrid ( module, clsGrid.first );


  cls.def("read",[](const IbraReader& self, const std::string& fileName, const bool trim = true,
                            std::array<int, 2> elevateDegree={0,0} ,
                            std::array<int, 2> preKnotRefine={0,0} ) {
    return self.read(fileName,trim,elevateDegree,preKnotRefine);
  },"fileName"_a, "trim"_a = true, "elevateDegree"_a = std::array<int, 2>{0,0}, "preKnotRefine"_a = std::array<int, 2>{0,0} );



}

}  // namespace Ikarus::Python

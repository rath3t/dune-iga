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
#include "modregisterGridView.hh"


namespace Dune::Python {




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
//  using GVTraits                = Dune::IGA::NurbsLeafGridViewTraits<Grid>;
//  auto typeName2 = GenerateTypeName( "NURBSGrid", dimension ,dimensionworld);
//  auto traits2 = insertClass< Grid >( scope, "NURBSGrid", typeName2, IncludeFiles{ "dune/iga/nurbsgrid.hh" } );
//  auto typeName = GenerateTypeName( "NurbsLeafGridViewTraits", MetaType< Grid >() );
//  auto traits = insertClass< GVTraits >( scope, "NurbsLeafGridViewTraits", typeName, IncludeFiles{ "dune/iga/nurbsleafgridview.hh" } );
//
//
//
////  if( iter.second )
////    Dune::Python::registerPyGridViewIterator< GridView, codim >(module);
//  typedef PyGridViewIterator< GridView, 0 > Iterator;
//  std::cout<<"Found in reg "<< Dune::Python::detail::findInTypeRegistry<Iterator>().second<<std::endl;
//
//

//  typedef PyGridViewIterator< GridView, 0 > Iterator;
//  std::cout<<"Found in reg "<< Dune::Python::detail::findInTypeRegistry<Iterator>().second<<std::endl;
//  typedef PyGridViewIterator< GridView, 1 > Iterator1;
//  std::cout<<"Found in reg "<< Dune::Python::detail::findInTypeRegistry<Iterator1>().second<<std::endl;




  cls.def(pybind11::init());
////
////  auto clsLeafView = insertClass< typename Grid::LeafGridView >( scope, "LeafGrid", GenerateTypeName( cls, "LeafGridView" ) );
////  if( clsLeafView.second )
////    registerGridView( scope, clsLeafView.first );
//
  auto includes = Dune::Python::IncludeFiles{"dune/iga/nurbsgrid.hh","dune/iga/io/ibra/ibrareader.hh","dune/iga/nurbsleafgridview.hh","dune/iga/nurbsgridleafiterator.hh"};
  auto clsName = Dune::className<IbraReader>();
  auto clsGrid = insertClass< Grid >( module, "NURBSGrid", Dune::Python::GenerateTypeName(clsName),
                                                     includes );
  if( clsGrid.second )
    registerHierarchicalGrid ( module, clsGrid.first );



//  auto clsLeafView = insertClass< typename Grid::LeafGridView >( scope, "LeafGrid", GenerateTypeName( cls, "LeafGridView" ) );
//  if( clsLeafView.second )
//    registerGridViewALT( scope, clsLeafView.first );
//
//  auto lv       = Dune::Python::insertClass<Grid>(
//      scope, "NURBSGrid",
//      Dune::Python::GenerateTypeName("Dune::IGA::NURBSGrid"),
//      includes)
//      .first;
//  lv.def("size",[](const Grid& self,int i) { return self.size(i); });
//  lv.def("globalRefine",[]( Grid& self,int i) { return self.globalRefine(i); });
//  lv.def_property_readonly( "leafView", pybind11::cpp_function( [] ( const Grid &self ) {
//                               return self.leafGridView();
//                             }, pybind11::keep_alive< 0, 1 >() ),
//                             R"doc(
//          Obtain leaf view of the grid
//          Returns:  leaf grid view
//        )doc" );

//  auto includes2 = Dune::Python::IncludeFiles{"dune/iga/nurbsgrid.hh","dune/iga/io/ibra/ibrareader.hh"};
//  Dune::Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [&] ( auto codim ) {
//
//    using IterType= PyGridViewIterator< GridView, codim >;
//    auto typeName = GenerateTypeName( "PyGridViewIterator", MetaType< GridView >(), codim );
//    auto iter = insertClass< IterType >( scope, "EntityIterator", typeName, IncludeFiles{ "dune/python/grid/range.hh" } );
//
//    if( iter.second )
//      Dune::Python::registerPyGridViewIterator< GridView, codim >(module);
//    typedef PyGridViewIterator< GridView, 0 > Iterator;
//    std::cout<<"Found in reg "<< Dune::Python::detail::findInTypeRegistry<Iterator>().second<<std::endl;
//
//  } );


  cls.def("read",[](const IbraReader& self, const std::string& fileName, const bool trim = true,
                            std::array<int, 2> elevateDegree={0,0} ,
                            std::array<int, 2> preKnotRefine={0,0} ) {
    return self.read(fileName,trim,elevateDegree,preKnotRefine).release();
  },"fileName"_a, "trim"_a = true, "elevateDegree"_a = std::array<int, 2>{0,0}, "preKnotRefine"_a = std::array<int, 2>{0,0} );



}

}  // namespace Ikarus::Python

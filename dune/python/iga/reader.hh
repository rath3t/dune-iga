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



template< class GridView, class... options >
inline static void registerGridViewImpl ( pybind11::handle scope, pybind11::class_< GridView, options... > cls )
{
  typedef typename GridView::Grid Grid;
  typedef PyGridViewIterator< GridView, 0 > PyElementIterator;

  using pybind11::operator""_a;

  pybind11::options opts;
  opts.disable_function_signatures();

  detail::registerGridViewConstructorFromGrid( cls, PriorityTag< 42 >() );

  cls.attr( "dimGrid" ) = pybind11::int_( static_cast< int >( GridView::dimension ) );
  cls.attr( "dimWorld" ) = pybind11::int_( static_cast< int >( GridView::dimensionworld ) );

//  registerGridEntities< GridView >( cls );
//  registerGridIntersection< GridView >( cls );
//
//  cls.def( "_mapper", [] ( GridView &self, pybind11::object layout ) {
//             return makeMultipleCodimMultipleGeomTypeMapper( self, layout );
//           }, pybind11::keep_alive< 0, 1 >(), "layout"_a,
//           R"doc(
//          Set up a mapper to attach data to the grid. The layout argument defines how many
//          degrees of freedom to assign to each subentity of a geometry type.
//
//          Args:
//              layout:     function, dict, tuple, or list defining the number of indices to reserve
//                          for each geometry type.
//
//          If layout is a dict, is must map geometry types to integers. All types not mentioned in
//          the dictionary are assumed to be zero.
//
//          If layout is a tuple or a list, it must contain exactly dimension+1 integers, one for
//          each codimension in the grid.
//
//          if layout is a function it maps a geometry type to the number of degrees of freedom to
//          reserve. Here a return value of 0 or `False` indicates that no data is to be attach, `True` can be used instead of 1.
//
//          Returns:   the mapper
//        )doc" );

//  // register iterators
//  Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [] ( auto codim ) {
//    registerPyGridViewIterator< GridView, codim >();
//  } );
//
//  if( Capabilities::canIterate< Grid, 0 >::value )
//    cls.def_property_readonly( "elements", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, 0 >( self ); },
//                               R"doc(
//            sequence of all elements (i.e., entities of codimension 0)
//          )doc" );
//  if( Capabilities::canIterate< Grid, 1 >::value )
//    cls.def_property_readonly( "facets", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, 1 >( self ); },
//                               R"doc(
//            range of all facets (i.e., entities of codimension 1)
//          )doc" );
//  if( Capabilities::canIterate< Grid, GridView::dimension-1 >::value )
//    cls.def_property_readonly( "edges", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, GridView::dimension-1 >( self ); },
//                               R"doc(
//            range of all edges (i.e., entities of dimension 1)
//          )doc" );
//  if( Capabilities::canIterate< Grid, GridView::dimension >::value )
//    cls.def_property_readonly( "vertices", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, GridView::dimension >( self ); },
//                               R"doc(
//            range of all vertices (i.e., entities of dimension 0)
//          )doc" );
//
//  std::array< pybind11::object (*) ( pybind11::object ), GridView::dimension+1 > makePyGridViewIterators;
//  Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &makePyGridViewIterators ] ( auto codim ) {
//    makePyGridViewIterators[ codim ] = makePyGridViewIterator< GridView, codim >;
//  } );
//  cls.def( "entities", [ makePyGridViewIterators ] ( pybind11::object self, int codim ) {
//             if( (codim < 0) || (codim > GridView::dimension) )
//               throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
//             return makePyGridViewIterators[ codim ]( self );
//           }, "codim"_a,
//           R"doc(
//          get range of entities for a codimension
//
//          Args:
//              codim:    Codimension to obtain range of entities for
//        )doc" );

//  registerPyIntersectionIterator< GridView >();
//  cls.def( "intersections", [] ( const GridView &self, const typename GridView::template Codim< 0 >::Entity &e ) {
//             return PyIntersectionIterator< GridView >( self.ibegin( e ), self.iend( e ) );
//           }, pybind11::keep_alive< 0, 1 >(), "element"_a,
//           R"doc(
//          get range of all codim-1 intersections for an element
//
//          Args:
//              element:    Element of obtain intersections for
//        )doc" );
//
//  registerPyBoundaryIntersectionIterator< GridView, PyElementIterator >();
//  cls.def_property_readonly( "boundaryIntersections", [] ( const GridView &self ) {
//                               return PyBoundaryIntersectionIterator< GridView, PyElementIterator >( self, PyElementIterator( self.template begin< 0 >(), self.template end< 0 >() ) );
//                             }, pybind11::keep_alive< 0, 1 >(),
//                             R"doc(
//          range of all codim-1 boundary intersections of the grid
//        )doc" );

  // register partitions

//  registerGridViewPartition< GridView, Interior_Partition >();
//  registerGridViewPartition< GridView, InteriorBorder_Partition >();
//  registerGridViewPartition< GridView, Overlap_Partition >();
//  registerGridViewPartition< GridView, OverlapFront_Partition >();
//  registerGridViewPartition< GridView, All_Partition >();
//  registerGridViewPartition< GridView, Ghost_Partition >();

//  cls.def_property_readonly( "interiorPartition", [] ( pybind11::object self ) {
//    return GridViewPartition< GridView, Interior_Partition >( self );
//  }, R"doc(
//          Partition containing only interior entities.
//        )doc" );
//  cls.def_property_readonly( "interiorBorderPartition", [] ( pybind11::object self ) {
//    return GridViewPartition< GridView, InteriorBorder_Partition >( self );
//  }, R"doc(
//          Partition containing only interior and border entities.
//        )doc" );
//  cls.def_property_readonly( "overlapPartition", [] ( pybind11::object self ) {
//    return GridViewPartition< GridView, Overlap_Partition >( self );
//  }, R"doc(
//          Partition containing only interior, border and overlap entities.
//        )doc" );
//  cls.def_property_readonly( "overlapFrontPartition", [] ( pybind11::object self ) {
//    return GridViewPartition< GridView, OverlapFront_Partition >( self );
//  }, R"doc(
//          Partition containing only interior, border, overlap, and front entities.
//        )doc" );
//  cls.def_property_readonly( "allPartition", [] ( pybind11::object self ) {
//    return GridViewPartition< GridView, All_Partition >( self );
//  }, R"doc(
//          Partition containing all entities.
//        )doc" );
//  cls.def_property_readonly( "ghostPartition", [] ( pybind11::object self ) {
//    return GridViewPartition< GridView, Ghost_Partition >( self );
//  }, R"doc(
//          Partition containing only ghost entities.
//        )doc" );

//  cls.def("__repr__",
//          [] (const GridView &gridView) -> std::string {
//            return "LeafGrid with " + std::to_string(gridView.indexSet().size(0)) + " elements";
//          });
//
//  cls.def_property_readonly( "hierarchicalGrid", [] ( const GridView &self ) -> const Grid & { return self.grid(); },
//                             R"doc(
//          associated hierarchical grid
//        )doc" );
//
  cls.def_property_readonly_static( "dimension", [] ( pybind11::object ) { return int(GridView::dimension); } );
  cls.def_property_readonly_static( "dimensionworld", [] ( pybind11::object ) { return int(GridView::dimensionworld); } );
//
//  cls.def( "size", [] ( const GridView &self, int codim ) {
//             if( (codim < 0) || (codim > GridView::dimension) )
//               throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
//             return self.size( codim );
//           }, "codim"_a,
//           R"doc(
//          Args:
//              codim:     required codimension
//          Returns:       number of subentities of codimension `codim`
//        )doc" );
//  cls.def( "size", [] ( const GridView &self, Dune::GeometryType gt ) {
//             if( (gt.dim() < 0) || (gt.dim() > GridView::dimension) )
//               throw pybind11::value_error( "Invalid geometry type (dimension must be in [0, " + std::to_string( GridView::dimension ) + "])." );
//             return self.size( gt );
//           }, "gt"_a,
//           R"doc(
//          Args:
//              gt:        a geometry type
//          Returns:       number of subentities of the given geometry type
//        )doc" );
//
//  registerVTKWriter< GridView >( cls );
//  cls.def( "vtkWriter", [] ( const GridView &self, const bool nonconforming = false ) {
//    const VTK::DataMode dm = nonconforming ? VTK::nonconforming : VTK::conforming;
//    return new VTKWriter< GridView >( self, dm );
//  }, pybind11::keep_alive< 0, 1 >() );
//  cls.def( "vtkWriter", [] ( const GridView &self, int subsampling ) {
//    return new SubsamplingVTKWriter< GridView >( self,
//                                                 Dune::refinementIntervals(1<<subsampling) );
//  }, pybind11::keep_alive< 0, 1 >(), "subsampling"_a );

//  cls.def( "overlapSize", [] ( const GridView &self, int codim ) {
//    if( (codim < 0) || (codim > GridView::dimension) )
//      throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
//    return self.overlapSize( codim );
//  }, "codim"_a );
//  cls.def( "ghostSize", [] ( const GridView &self, int codim ) {
//    if( (codim < 0) || (codim > GridView::dimension) )
//      throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
//    return self.ghostSize( codim );
//  }, "codim"_a );

//  cls.def_property_readonly( "_indexSet", [] ( const GridView &self ) -> const typename GridView::IndexSet & {
//                               return self.indexSet();
//                             }, pybind11::return_value_policy::reference, pybind11::keep_alive< 0, 1 >(),
//                             R"doc(
//          index set for the grid
//        )doc" );

//  cls.def_property_readonly( "comm", [] ( const GridView &gridView ) -> const typename Grid::CollectiveCommunication & {
//                               return gridView.grid().comm();
//                             }, pybind11::return_value_policy::reference, pybind11::keep_alive< 0, 1 >(),
//                             R"doc(
//          collective communication for the grid
//
//          Note: For collective (or global) operations, all processes in this
//                collective communication must call the corresponding method.
//        )doc" );
//
//  typedef NumPyCommDataHandle< MultipleCodimMultipleGeomTypeMapper< GridView >, double, std::function< double ( double, double ) > > CommDataHandle;
//  cls.def( "communicate", [] ( const GridView &gridView, CommDataHandle &dataHandle, InterfaceType iftype, CommunicationDirection dir ) {
//    gridView.communicate( dataHandle, iftype, dir );
//  } );
//  cls.def( "communicate", [] ( const GridView &gridView, pybind11::object dataHandle, InterfaceType iftype, CommunicationDirection dir ) {
//    ProxyDataHandle proxyDataHandle( std::move( dataHandle ) );
//    gridView.communicate( proxyDataHandle, iftype, dir );
//  });

  // export grid capabilities

//  cls.def_property_readonly( "conforming", [] ( pybind11::object ) { return static_cast< bool >( GridView::conforming ); } );
//
//  if( Capabilities::hasSingleGeometryType< Grid >::v )
//    cls.def_property_readonly_static( "type", [] ( pybind11::object ) {
//      return GeometryType( Capabilities::hasSingleGeometryType< Grid >::topologyId, Grid::dimension );
//    } );
//
//  cls.def_property_readonly_static( "isCartesian", [] ( pybind11::object ) { return Capabilities::isCartesian< Grid >::v; } );
//  cls.def_property_readonly_static( "canCommunicate", [] ( pybind11::object ) {
//    pybind11::tuple canCommunicate( Grid::dimension+1 );
//    Hybrid::forEach( std::make_integer_sequence< int, Grid::dimension+1 >(), [ &canCommunicate ] ( auto codim ) {
//      canCommunicate[ codim ] = pybind11::cast( bool( Capabilities::canCommunicate< Grid, codim >::v ) );
//    } );
//    return canCommunicate;
//  } );

//  cls.def_property_readonly_static( "threadSafe", [] ( pybind11::object ) { return Capabilities::viewThreadSafe< Grid >::v; } );

  // export utility methods

//  cls.def( "coordinates", [] ( const GridView &self ) { return coordinates( self ); },
//           R"doc(
//          Returns: `numpy` array with the coordinates of all vertices in the grid in
//                   the format `[ [x_1,y_1], [x_2,y_2], ..., [x_N,y_N] ]` for example
//                   in 2d.
//        )doc" );
//  cls.def( "tessellate", [] ( const GridView &self, int level ) { return tessellate( self, level ); }, "level"_a = 0,
//           R"doc(
//          Generated a possibly refined tessellation using only simplices.
//
//          Args:
//              level: virtual refinement level to use to generate the tessellation
//
//          Returns: (coordinates,simplices) where coordinates is a `numpy` array
//                   of the vertex coordinates
//                   (e.g. in 2d `[ [x_1,y_1], [x_2,y_2], ..., [x_N,y_N] ]` )
//                   and simplices is a `numpy` array of the vertices of the simplices
//                   (e.g. in 2d `[s_11,s_12,s_13], [s_21,s_22,s_23], ..., [s_N1,s_N2,s_N3] ]` )
//
//        )doc" );
//  cls.def( "polygons", [] ( const GridView &self ) { return polygons( self ); },
//           R"doc(
//          Store the grid in numpy arrays.
//
//          Returns: coordinate array storing the vertex coordinate of each polygon
//                   in the grid.
//        )doc" );
//
//  cls.def( "contains", [] ( GridView &self, pybind11::object entity ) {
//             bool found = false, contained = false;
//             Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &self, entity, &found, &contained ] ( auto codim ) {
//               typedef typename GridView::template Codim< decltype( codim )::value >::Entity Entity;
//               if( pybind11::isinstance< Entity >( entity ) )
//               {
//                 found = true;
//                 contained = self.contains( pybind11::cast< const Entity & >( entity ) );
//               }
//             } );
//             if( found )
//               return contained;
//             else
//               throw pybind11::value_error( "Argument 'entity' is not a valid entity for this grid." );
//           }, "entity"_a,
//           R"doc(
//          Check whether an entity is contained in the grid instance
//
//          Args:
//              entity:   entity to check
//
//          Note:
//          - The entity must be contained in the corresponding hierarchical grid.
//        )doc" );
//
//#if HAVE_DUNE_VTK
//  using VirtualizedGF = Dune::Vtk::Function<GridView>;
//      auto vgfClass = Python::insertClass<VirtualizedGF>(scope,"VtkFunction",
//          Python::GenerateTypeName("Dune::Vtk::Function", MetaType<GridView>()),
//          Python::IncludeFiles{"dune/vtk/function.hh"});
//      if( vgfClass.second )
//      {
//        vgfClass.first.def("name",[](VirtualizedGF &self) { return self.name(); });
//      }
//#endif
  auto addAttr = pybind11::module::import( "dune.grid.grid_generator" ).attr("addAttr");
  addAttr(scope, cls);
}



template <class Grid, class... options>  requires(IsSpecializationTwoNonTypesAndType<Dune::IGA::NURBSGrid,Grid>::value)
void registerHierarchicalGridImpl(pybind11::module module, pybind11::class_<Grid, options...> cls) {
  pybind11::module::import( "dune.geometry" );
  pybind11::module::import( "dune.grid" );

  using pybind11::operator""_a;

  pybind11::options opts;
  opts.disable_function_signatures();

  auto clsLeafView = insertClass< typename Grid::LeafGridView >( module, "LeafGrid", GenerateTypeName( cls, "LeafGridView" ) );
  if( clsLeafView.second )
    registerGridViewImpl( module, clsLeafView.first );

//  module.def( "reader", [] ( const std::tuple< Reader, std::string > &args ) { return reader< Grid >( args ); } );
//  module.def( "reader", [] ( const std::string &args ) { return reader< Grid >( std::make_tuple( Reader::dgf,args ) ); } );
//  module.def( "reader", [] ( const StructuredReader<Grid> &args ) { return reader< Grid >( args ); } );
//  module.def( "reader", [] ( const pybind11::dict &args ) { return reader< Grid >( args ); } );

//  registerHierarchicalGridPicklingSupport( cls );
//
//  registerHierarchicalGridIdSets( cls );

  cls.def_property_readonly( "leafView", pybind11::cpp_function( [] ( const Grid &self ) {
                               return self.leafGridView();
                             }, pybind11::keep_alive< 0, 1 >() ),
                             R"doc(
          Obtain leaf view of the grid

          Returns:  leaf grid view
        )doc" );
//  cls.def( "_levelView", [] ( const Grid &self, int level ) {
//             return self.levelGridView( level );
//           }, pybind11::keep_alive< 0, 1 >(), "level"_a,
//           R"doc(
//          Obtain level view of the grid
//
//          Args:
//              level:    level to obtain view for
//
//          Returns:  level grid view
//        )doc" );

//  typedef typename Grid::template Codim< 0 >::Entity Element;
//  cls.def( "mark", [] ( Grid &self, const Element &element, Marker marker ) {
//             self.mark( static_cast< int >( marker ), element );
//           }, "element"_a, "marker"_a,
//           R"doc()doc" );
//  cls.def( "mark", [] ( Grid &self, const std::function< Marker( const Element &e ) > &marking ) {
//             std::pair< int, int > marked;
//             for( const Element &element : elements( self.leafGridView() ) )
//             {
//               Marker marker = marking( element );
//               marked.first += static_cast< int >( marker == Marker::Refine );
//               marked.second += static_cast< int >( marker == Marker::Coarsen );
//               self.mark( static_cast< int >( marker ), element );
//             }
//             return marked;
//           }, "marking"_a,
//           R"doc(
//          Set the grid's adaptation markers
//
//          Args:
//              marking:    callback returning a dune.grid.Marker for each leaf
//                          element in the grid
//        )doc" );
//
//  cls.def( "adapt", [] ( Grid &self ) {
//             const auto &range = detail::gridModificationListenersRange(self);
//             for( const auto &listener : range )
//               listener->preModification( self );
//             self.preAdapt();
//             self.adapt();
//             self.postAdapt();
//             for( const auto &listener : range )
//               listener->postModification( self );
//           },
//           R"doc(
//          Refine or coarsen the hierarchical grid to match the current marking
//
//          All elements marked for refinement will be refined by this operation.
//          However, due to closure rules, additional elements might be refined.
//          Similarly, not all elements marked for coarsening are necessarily
//          coarsened.
//
//          Note:
//          - This is a collective operation.
//          - The grid implementation defines the rule by which elements are split.
//        )doc" );
//
//  cls.def( "adapt", [] ( Grid &self, const std::function< Marker( const Element &e ) > &marking ) {
//             std::pair< int, int > marked;
//             for( const Element &element : elements( self.leafGridView() ) )
//             {
//               Marker marker = marking( element );
//               marked.first += static_cast< int >( marker == Marker::Refine );
//               marked.second += static_cast< int >( marker == Marker::Coarsen );
//               self.mark( static_cast< int >( marker ), element );
//             }
//             if (marked.first + marked.second)
//             {
//               const auto &range = detail::gridModificationListenersRange(self);
//               for( const auto &listener : range )
//                 listener->preModification( self );
//               self.preAdapt();
//               self.adapt();
//               self.postAdapt();
//               for( const auto &listener : range )
//                 listener->postModification( self );
//             }
//             return marked;
//           },
//           R"doc(
//          Refine or coarsen the hierarchical grid to match the provided marking function.
//
//          Args:
//              marking:    callback returning a dune.grid.Marker for each leaf
//                          element in the grid
//
//          All elements for which are marked for refinement by the callback function
//          will be refined by this operation.
//          However, due to closure rules, additional elements might be refined.
//          Similarly, not all elements marked for coarsening are necessarily
//          coarsened.
//
//          Note:
//          - This is a collective operation.
//          - The grid implementation defines the rule by which elements are split.
//        )doc" );

  cls.def( "globalRefine", [] ( Grid &self, int level ) {
             const auto &range = detail::gridModificationListenersRange(self);
             for( const auto &listener : range )
               listener->preModification( self );
             self.globalRefine( level );
             for( const auto &listener : range )
               listener->postModification( self );
           }, "iterations"_a = 1,
           R"doc(
          Refine each leaf element of the grid.

          Args:
              iterations:   Number of global refinement iterations to perform (defaults to 1)

          Note:
          - This is a collective operation.
          - The grid implementation defines the rule by which elements are split.
        )doc" );

//  cls.def( "loadBalance", [] ( Grid &self ) {
//             const auto &range = detail::gridModificationListenersRange(self);
//             for( const auto &listener : range )
//               listener->preModification( self );
//             self.loadBalance();
//             for( const auto &listener : range )
//               listener->postModification( self );
//           },
//           R"doc(
//          Redistribute the grid to equilibrate the work load on each process.
//
//          Note:
//          - This is a collective operation.
//          - The redistribution strategy is chosen by the grid implementation.
//        )doc" );
//
//  cls.def_property_readonly( "maxLevel", [] ( const Grid &self ) -> int { return self.maxLevel(); } );
  cls.def_property_readonly_static( "dimension", [] ( pybind11::object ) { return int(Grid::dimension); } );
  cls.def_property_readonly_static( "dimensionworld", [] ( pybind11::object ) { return int(Grid::dimensionworld); } );
//  // evaluate this information at grid creation time since this changes for
//  // ALUGrid conform when used as simplex grid.
//  const int refineStepsForHalf = DGFGridInfo< Grid >::refineStepsForHalf();
//  cls.def_property_readonly_static( "refineStepsForHalf", [refineStepsForHalf] ( pybind11::object ) { return refineStepsForHalf; } );

  // export grid capabilities

//  if( Capabilities::hasSingleGeometryType< Grid >::v )
//  {
//    cls.def_property_readonly_static( "type", [] ( pybind11::object ) {
//                                        return GeometryType( Capabilities::hasSingleGeometryType< Grid >::topologyId, Grid::dimension );
//                                      },
//                                      R"doc(
//            "All elements in this grid have this geometry type"
//          )doc" );
//  }

//  cls.def_property_readonly_static( "isCartesian", [] ( pybind11::object ) { return Capabilities::isCartesian< Grid >::v; } );
//  cls.def_property_readonly_static( "canCommunicate", [] ( pybind11::object ) {
//    pybind11::tuple canCommunicate( Grid::dimension+1 );
//    Hybrid::forEach( std::make_integer_sequence< int, Grid::dimension+1 >(), [ &canCommunicate ] ( auto codim ) {
//      canCommunicate[ codim ] = pybind11::cast( bool( Capabilities::canCommunicate< Grid, codim >::v ) );
//    } );
//    return canCommunicate;
//  } );

//  cls.def_property_readonly_static( "threadSafe", [] ( pybind11::object ) { return Capabilities::threadSafe< Grid >::v; } );
//  cls.def_property_readonly_static( "viewThreadSafe", [] ( pybind11::object ) { return Capabilities::viewThreadSafe< Grid >::v; } );
//
//  auto [ clsComm, notRegistered ] = insertClass< typename Grid::CollectiveCommunication >( cls, "CollectiveCommunication", GenerateTypeName( cls, "CollectiveCommunication" ) );
//  if( notRegistered )
//    registerCommunication( clsComm );
//
//  cls.def_property_readonly( "comm", [] ( const Grid &grid ) -> const typename Grid::CollectiveCommunication & {
//                               return grid.comm();
//                             }, pybind11::return_value_policy::reference, pybind11::keep_alive< 0, 1 >(),
//                             R"doc(
//          collective communication for the grid
//
//          Note: For collective (or global) operations, all processes in this
//                collective communication must call the corresponding method.
//        )doc" );

  auto addHAttr = pybind11::module::import( "dune.grid.grid_generator" ).attr("addHAttr");
  addHAttr(module);

}

template <class NURBSGrid, class... options>  requires(IsSpecializationTwoNonTypesAndType<Dune::IGA::NURBSGrid,NURBSGrid>::value)
void registerHierarchicalGrid(pybind11::module module, pybind11::class_<NURBSGrid, options...> cls) {
  using pybind11::operator""_a;

static constexpr std::integral auto dimension      = NURBSGrid::dimension;
static constexpr std::integral auto dimensionworld = NURBSGrid::dimensionworld;
using ctype                                        = typename NURBSGrid::ctype;

  module.def( "reader", [] ( const pybind11::dict &args_ ) { return Dune::Python::IGA::reader< NURBSGrid >( args_ ); } );

    registerHierarchicalGridImpl (module, cls);

//  auto clsLeafView = insertClass< typename NURBSGrid::LeafGridView >( module, "LeafGrid", GenerateTypeName( cls, "LeafGridView" ) );
//  if( clsLeafView.second )
//  registerGridView( module, clsLeafView.first );

//  clsLeafView.first.def("preBasis",[](const typename NURBSGrid::LeafGridView& self){return self.impl().preBasis();});


using ControlPointNetType    = typename NURBSGrid::ControlPointNetType;
  using NURBSPatchDataType    = typename NURBSGrid::NURBSPatchDataType;

  cls.def(pybind11::init([](const std::array<std::vector<double>, dimension>& knotSpans, const ControlPointNetType& controlPoints,
                            const std::array<int, dimension>& order){return new NURBSGrid(knotSpans,controlPoints,order);}));

  cls.def(pybind11::init([](const NURBSPatchDataType& nurbsPatchData){return new NURBSGrid(nurbsPatchData);}));
  cls.def("globalRefineInDirection",[]( NURBSGrid& self,const int dir, const int refinementLevel, bool omitTrim = false){self.globalRefineInDirection(dir,refinementLevel,omitTrim);});
  cls.def("patchData",[](const NURBSGrid& self,int i = 0){return self.patchData(i);});




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

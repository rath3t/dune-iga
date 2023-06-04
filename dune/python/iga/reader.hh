// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include <dune/python/grid/enums.hh>
#include <dune/iga/io/ibra/ibrareader.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/grid/capabilities.hh>

namespace Dune::Python::IGA
{
enum class Reader { json };

} // namespace Dune

namespace Dune::Python
{
template <int dim, int dimworld, typename ScalarType>
class NURBSGrid;
//template< class Grid, typename In >
//inline static std::shared_ptr< Grid > readJSON ( In &input )
//{
//  JSONGridFactory< Grid > dgfFactory( input );
//  std::shared_ptr< Grid > grid( dgfFactory.grid() );
//  grid->loadBalance();
//  return grid;
//}

//template <std::integral auto dim, auto ReaderType, std::integral auto dimworld, typename ScalarType> requires (std::is_same_v<Reader, decltype(ReaderType)> or std::is_same_v<Dune::Python::IGA::Reader,decltype(ReaderType)>)
//inline static std::shared_ptr< Dune::IGA::NURBSGrid<dim, dimworld, ScalarType> > reader ( const std::tuple< Reader, std::string > &args )
//{
//  return nullptr;
//}
//
//template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
//inline static std::shared_ptr< Dune::IGA::NURBSGrid<dim, dimworld, ScalarType> > reader ( const std::tuple< Dune::Python::IGA::Reader, std::string > &args )
//{
//
//    switch (std::get<0>(args)) {
//      std::cout << "Reader::json" << std::endl;
//      case Dune::Python::IGA::Reader::json:
//        return readJSON<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>>(std::get<1>(args));
//
//      default:return nullptr;
//    }
//}
// we have to ahead of https://gitlab.dune-project.org/core/dune-grid/-/blob/releases/2.9/dune/python/grid/hierarchical.hh?ref_type=heads#L233 thus we provide  Dune::PriorityTag< 42 > and
//make sure this type is used if an iga grid is passed this is function is enabled and used by adl
template <template <auto,auto, typename> class Type, typename>
struct IsSpecializationTwoNonTypesAndType : std::false_type {};

template <template <auto,auto, typename> class Type, auto T, auto T2, typename S>
struct IsSpecializationTwoNonTypesAndType<Type, Type<T, T2,S>> : std::true_type {};

template <int dim, int dimworld, typename ScalarType>
struct Capabilities::HasGridFactory<Dune::IGA::NURBSGrid<dim,dimworld,ScalarType>> : public std::integral_constant< bool, false >
{};


   template< class Grid> requires(IsSpecializationTwoNonTypesAndType<Dune::IGA::NURBSGrid,Grid>::value)
inline static std::shared_ptr< Grid > reader ( const pybind11::dict &dict)
{

  std::cout << "Reader::Tuple" << std::endl;
  std::string file_path;
  Dune::Python::IGA::Reader reader=IGA::Reader::json;
  bool trim = true;
  std::array<int, 2> elevateDegree = {0, 0};
  std::array<int, 2> preKnotRefine = {0, 0};
  std::array<int, 2> postKnotRefine = {0, 0};
  if (dict.contains("reader"))
    reader =  dict["reader"].cast<Dune::Python::IGA::Reader>();

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
      if (dict.contains("elevate_degree"))
        postKnotRefine = dict["post_knot_refine"].cast<std::array<int, 2>>();

      static constexpr std::integral auto dim = Grid::dimension;
      static constexpr std::integral auto dimworld = Grid::dimensionworld;
      using ScalarType = typename Grid::ctype;
      return Dune::IGA::IbraReader < dim, dimworld, ScalarType
          > ::read(file_path, trim, elevateDegree, preKnotRefine, postKnotRefine);
      break;
    default:
      DUNE_THROW(Dune::NotImplemented, "Your requested reader is not implemeneted");
      return nullptr;
      break;
  }
}
}

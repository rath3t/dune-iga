// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include "dune/python/iga/gridenums.hh"
#include <dune/iga/hierarchicpatch/patchgrid.hh>
#include <dune/iga/io/ibrareader.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <dune/iga/trimmer/identitytrimmer/trimmer.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/grid/capabilities.hh>
#include <dune/python/grid/enums.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/pybind11/pybind11.h>

#if HAVE_DUNE_VTK
  #include <dune/iga/io/vtk/igadatacollector.hh>
  #include <dune/python/vtk/writer.hh>
  #include <dune/vtk/writers/vtkunstructuredgridwriter.hh>
#endif

namespace Dune {
template <int dim, int dimworld, template <int, int, typename> typename GridFamily_, typename ScalarType>
struct DGFGridFactory<IGA::PatchGrid<dim, dimworld, GridFamily_, ScalarType>>
{
};

  // DGFGridInfo
  // -----------

template <int dim, int dimworld, template <int, int, typename> typename GridFamily_, typename ScalarType>
  struct DGFGridInfo< IGA::PatchGrid<dim, dimworld, GridFamily_, ScalarType> >
  {
    static int refineStepsForHalf ()
    {
      return 1;
    }

    static double refineWeight ()
    {
      return 0.5;
    }
  };
}

namespace Dune::Python {

#ifndef DOXYGEN
// we have to ahead of
// https://gitlab.dune-project.org/core/dune-grid/-/blob/releases/2.9/dune/python/grid/hierarchical.hh?ref_type=heads#L233
// thus we use requires(IsSpecializationTwoNonTypesAndType<Dune::IGA::NURBSGrid,Grid>::value) to be sure this overload
// is used
// make sure this type is used if an iga grid is passed this is function is enabled and used by adl

template <int dim, int dimworld, template <int, int, typename> typename GridFamily_, typename ScalarType>
struct Capabilities::HasGridFactory<Dune::IGA::PatchGrid<dim, dimworld, GridFamily_, ScalarType>>
    : public std::integral_constant<bool, false>
{
};

template <template <auto, auto, template <int, int, typename> typename, typename> class Type, typename>
struct IsSpecializationTwoNonTypesTemplateAndType : std::false_type
{
};

template <template <auto, auto, template <int, int, typename> typename, typename> class Type, auto T, auto T2,
          template <int, int, typename> typename GF, typename S>
struct IsSpecializationTwoNonTypesTemplateAndType<Type, Type<T, T2, GF, S>> : std::true_type
{
};

} // namespace Dune::Python
#endif
namespace Dune::Python::IGA {
template <class Grid>
requires(IsSpecializationTwoNonTypesTemplateAndType<Dune::IGA::PatchGrid, Grid>::value)
inline static std::shared_ptr<Grid> reader(const pybind11::dict& dict) {
  std::string file_path;
  Dune::Python::IGA::Reader reader = IGA::Reader::json;
  bool trim                        = true;
  std::array<int, 2> degreeElevate = {0, 0};
  std::array<int, 2> preKnot       = {0, 0};
  if (dict.contains("reader"))
    reader = dict["reader"].cast<Dune::Python::IGA::Reader>();

  switch (reader) {
    case IGA::Reader::json: {
      if (dict.contains("file_path"))
        file_path = dict["file_path"].cast<std::string>();
      else
        DUNE_THROW(Dune::IOError, "No field in dict with name file_path. Unable to read grid");
      if (dict.contains("trim"))
        trim = dict["trim"].cast<bool>();
      if (dict.contains("degree_elevate"))
        degreeElevate = dict["degree_elevate"].cast<std::array<int, 2>>();
      if (dict.contains("pre_knot_refine"))
        preKnot = dict["pre_knot_refine"].cast<std::array<int, 2>>();

      using ScalarType                             = typename Grid::ctype;
      using GridFactory                            = Dune::GridFactory<Grid>;

      auto gridFactory = GridFactory();
      gridFactory.insertJson(file_path, trim, preKnot, degreeElevate);

      return gridFactory.createGrid();
    }
    default:
      DUNE_THROW(Dune::NotImplemented, "Your requested reader is not implemented");
  }
}

template <class Grid, class... options>
requires(IsSpecializationTwoNonTypesTemplateAndType<Dune::IGA::PatchGrid, Grid>::value)
void registerHierarchicalGrid(pybind11::module module, pybind11::class_<Grid, options...> cls) {
  using pybind11::operator""_a;

  static constexpr std::integral auto dimension      = Grid::dimension;
  static constexpr std::integral auto dimensionworld = Grid::dimensionworld;
  using ctype                                        = typename Grid::ctype;

  if constexpr (dimension == 2)
    module.def("reader", [](const pybind11::dict& args_) { return Dune::Python::IGA::reader<Grid>(args_); });

  Dune::Python::registerHierarchicalGrid(module, cls);

  using LeafGridView = typename Grid::LeafGridView;

  auto clsLeafView = insertClass<LeafGridView>(module, "LeafGrid", GenerateTypeName(cls, "LeafGridView"));
  if (clsLeafView.second)
    registerGridView(module, clsLeafView.first);

#if HAVE_DUNE_VTK
  if constexpr (dimension == 2) {
    pybind11::module::import("dune.vtk");

    using TrimmedWriterType =
        Dune::Vtk::UnstructuredGridWriter<LeafGridView, Dune::Vtk::DiscontinuousIgaDataCollector<LeafGridView>>;
    auto clsLeafViewWriter = insertClass<TrimmedWriterType>(module, "TrimmedVtkWriter",
                                                            GenerateTypeName(clsLeafView.first, "TrimmedVtkWriter"));
    if (clsLeafViewWriter.second)
      Dune::Vtk::registerVtkWriter<TrimmedWriterType>(module, clsLeafViewWriter.first);

    clsLeafView.first.def(
        "trimmedVtkWriter",
        [](const LeafGridView& self, int subSample = 0) {
          auto dataCollector =
              std::make_shared<Dune::Vtk::DiscontinuousIgaDataCollector<LeafGridView>>(self, subSample);
          return new Dune::Vtk::UnstructuredGridWriter(dataCollector, Vtk::FormatTypes::ASCII);
        },
        pybind11::arg("subSample") = 0);
  }
#endif

  using PatchData = Dune::IGA::NURBSPatchData<dimension, dimensionworld, ctype>;

  using ControlPointNetType = typename PatchData::ControlPointNetType;

  // cls.def(pybind11::init(
  //     [](const std::array<std::vector<double>, dimension>& knotSpans, const ControlPointNetType& controlPoints,
  //        const std::array<int, dimension>& order) { return new Grid(knotSpans, controlPoints, order); }));

  cls.def(pybind11::init([](const PatchData& nurbsPatchData) { return new Grid(nurbsPatchData); }));

  cls.def("patchData", [](Grid& self){return self.patchGeometryAtBack().patchData();});

  cls.def("globalRefine", [](Grid& self, int refCount) { self.globalRefine(refCount); }, pybind11::arg("refCount"));

  cls.def(
      "globalRefineInDirection",
      [](Grid& self, const std::array<int, dimension>& s) { self.globalRefineInDirection(s); }, pybind11::arg("s"));

  cls.def(
      "degreeElevate",
      [](Grid& self, const std::array<int, dimension>& elevationFactors, int lvl) {
        self.degreeElevate(elevationFactors, lvl);
      },
      pybind11::arg("elevationFactors"), pybind11::arg("lvl"));

  cls.def(
      "degreeElevateOnAllLevels",
      [](Grid& self, const std::array<int, dimension>& elevationFactors) {
        self.degreeElevateOnAllLevels(elevationFactors);
      },
      pybind11::arg("elevationFactors"));
}
} // namespace Dune::Python::IGA

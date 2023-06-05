// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <dune/python/pybind11/pybind11.h>
#include "dune/python/iga/reader.hh"
#include "dune/iga/nurbsgrid.hh"

PYBIND11_MODULE( _iga, m )
{
  pybind11::enum_< Dune::Python::IGA::Reader > reader( m, "reader" );
  reader.value( "json", Dune::Python::IGA::Reader::json );

  m.def("preBasis",[](const typename Dune::IGA::NURBSGrid<1,1,double>::GridView & gridView,int i=0 ) { return gridView.impl().preBasis();});
  m.def("preBasis",[](const typename Dune::IGA::NURBSGrid<1,2,double>::GridView & gridView,int i=0 ) { return gridView.impl().preBasis();});
  m.def("preBasis",[](const typename Dune::IGA::NURBSGrid<2,2,double>::GridView & gridView,int i=0 ) { return gridView.impl().preBasis();});
  m.def("preBasis",[](const typename Dune::IGA::NURBSGrid<2,3,double>::GridView & gridView,int i=0 ) { return gridView.impl().preBasis();});
  m.def("preBasis",[](const typename Dune::IGA::NURBSGrid<3,3,double>::GridView & gridView,int i=0 ) { return gridView.impl().preBasis();});

}

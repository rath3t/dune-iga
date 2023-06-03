// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <dune/python/pybind11/pybind11.h>
//#include "dune/iga/controlpoint.hh"

PYBIND11_MODULE( _iga, m )
{
  // enumeration types from dune-grid

//Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<dimension + 1>()), [&](const auto i) {
//std::get<i>(*entityVector_.get()).reserve(patch.size(i));
//for (unsigned int j = 0; j < patch.size(i); ++j) {
//std::get<i>(*entityVector_.get())
//.emplace_back(NURBSGridEntity<i, dimension, GridImpl>(*this, j, currentPatchId));
//}
//});
//
//using FEreq   = ControlPoint<;
//auto includes = Dune::Python::IncludeFiles{"ikarus/finiteElements/feRequirements.hh"};
//auto lv       = Dune::Python::insertClass<FEreq>(
//    m, "FErequirements", Dune::Python::GenerateTypeName("FErequirements<Eigen::Ref<Eigen::VectorXd>>"),
//    includes)
//    .first;
//lv.def(py::init());
//lv.def("addAffordance", [](FEreq& self, const ScalarAffordances& affordances) { self.addAffordance(affordances); });
//lv.def("addAffordance", [](FEreq& self, const VectorAffordances& affordances) { self.addAffordance(affordances); });
//lv.def("addAffordance", [](FEreq& self, const MatrixAffordances& affordances) { self.addAffordance(affordances); });
//lv.def(
//"insertGlobalSolution",
//[](FEreq& self, FESolutions solType, Ref<VectorXd> solVec) {
//self.insertGlobalSolution(std::move(solType), solVec);
//},
//"solutionType"_a, "solutionVector"_a.noconvert());
//lv.def(
//"getGlobalSolution", [](FEreq& self, FESolutions solType) { return self.getGlobalSolution(std::move(solType)); },
//py::return_value_policy::reference_internal);
//lv.def(
//"insertParameter",
//[](FEreq& self, FEParameter parType, ValueWrapper<double>& parVal) {
//self.insertParameter(std::move(parType), parVal.val);
//},
//py::keep_alive<1, 3>(), "FEParameter"_a, "parameterValue"_a.noconvert());
//
//lv.def("getParameter", [](const FEreq& self, FEParameter parType) { return self.getParameter(std::move(parType)); });

}

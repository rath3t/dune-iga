// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "dune/python/iga/gridenums.hh"
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE(_iga, m) {
  pybind11::enum_<Dune::Python::IGA::Reader> reader(m, "reader");
  reader.value("json", Dune::Python::IGA::Reader::json);
}

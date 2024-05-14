# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

find_package(Eigen3 3.3.9 REQUIRED)
include(AddEigenFlags)

find_package(nlohmann_json 3.11.2 REQUIRED)
include(AddnLohmannJsonFlags)

find_package(PkgConfig REQUIRED)
pkg_check_modules(Clipper2Z REQUIRED IMPORTED_TARGET Clipper2Z)

find_library(Clipper2Z_LIB Clipper2Z HINTS ${PKG_Clipper2Z_LIBDIR} REQUIRED)

if(Clipper2Z_LIB)
  message(STATUS "Clipper2Z_LIB_FOUND" ${Clipper2Z_LIB_FOUND})
  set(Clipper2Z_LIB_FOUND TRUE)
endif(Clipper2Z_LIB)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Clipper2Z DEFAULT_MSG Clipper2Z_INCLUDEDIR Clipper2Z_LIB)

include(AddClipperLibFlags)

find_package(earcut_hpp REQUIRED)
include(AddEarCutFlags)

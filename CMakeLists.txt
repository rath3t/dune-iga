# SPDX-FileCopyrightText: 2022 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-2.1-or-later

cmake_minimum_required(VERSION 3.16)
project(dune-iga CXX)

set(CMAKE_CXX_STANDARD 20)

if(NOT
   (dune-common_DIR
    OR dune-common_ROOT
    OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
  string(REPLACE ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
                 ${PROJECT_BINARY_DIR})
endif()

# find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
     ${dune-common_MODULE_PATH})

find_package(PkgConfig REQUIRED)
pkg_check_modules(Clipper2Lib REQUIRED IMPORTED_TARGET Clipper2)
find_package(nlohmann_json 3.11.2 REQUIRED)
find_package(earcut_hpp REQUIRED)
find_package(autodiff REQUIRED)
find_package(Matplot++ REQUIRED)
find_package(spdlog REQUIRED)


# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory(src)
add_subdirectory(dune)
add_subdirectory(cmake/modules)
add_subdirectory(cmake/FormatTarget)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

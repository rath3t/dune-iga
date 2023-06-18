# SPDX-FileCopyrightText: 2023 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.16)
project(
  dune-iga
  VERSION 0.1.7
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

option(DUNE_ENABLE_PYTHONBINDINGS "Enable Python bindings" ON)
option(DUNE_PYTHON_ALLOW_GET_PIP "Allow dune-common to install pip into venv"
       ON)
add_definitions(-DBUILD_SHARED_LIBS=1)
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

# silence alugrid
set(ENV{ALUGRID_VERBOSITY_LEVEL} 0)

# include the dune macros

include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

if(ikarus_FOUND
   AND Eigen3_FOUND
   AND autodiff_FOUND)
  add_subdirectory(src)
endif()
if(DUNE_ENABLE_PYTHONBINDINGS)
  add_subdirectory(python)
endif()
add_subdirectory(dune)
add_subdirectory(cmake)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

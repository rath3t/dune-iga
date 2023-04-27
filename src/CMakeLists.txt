# SPDX-FileCopyrightText: 2022 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-2.1-or-later
add_definitions(-DDUNE_LOCALFEFUNCTIONS_USE_EIGEN=1)

add_executable("dune-iga_example" dune-iga_example.cc)

target_link_dune_default_libraries("dune-iga_example")
target_link_libraries(dune-iga_example PRIVATE  PkgConfig::Clipper2Lib nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp Matplot++::matplot)

# Copy the example grid used for testing into the build dir
file(COPY auxiliaryFiles/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/auxiliaryFiles)
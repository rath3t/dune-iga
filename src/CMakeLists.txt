# SPDX-FileCopyrightText: 2022 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-2.1-or-later
add_definitions(-DDUNE_LOCALFEFUNCTIONS_USE_EIGEN=1)

add_executable("Example-1" dune-iga_ex1.cc)
target_link_dune_default_libraries("Example-1")
target_link_libraries(Example-1 PRIVATE  PkgConfig::Clipper2Lib nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp Matplot++::matplot)

add_executable("Example-1a" dune-iga_ex1a.cc)
target_link_dune_default_libraries("Example-1a")
target_link_libraries(Example-1a PRIVATE  PkgConfig::Clipper2Lib nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp Matplot++::matplot)

add_executable("Example-2" dune-iga_ex2.cc)
target_link_dune_default_libraries("Example-2")
target_link_libraries(Example-2 PRIVATE  PkgConfig::Clipper2Lib nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp Matplot++::matplot)

add_executable("Example-3" dune-iga_ex3.cc)
target_link_dune_default_libraries("Example-3")
target_link_libraries(Example-3 PRIVATE  PkgConfig::Clipper2Lib nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp Matplot++::matplot)

add_executable("Example-4" dune-iga_ex4.cc)
target_link_dune_default_libraries("Example-4")
target_link_libraries(Example-4 PRIVATE  PkgConfig::Clipper2Lib nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp Matplot++::matplot)

add_executable("Example-5" dune-iga_ex5.cc)
target_link_dune_default_libraries("Example-5")
target_link_libraries(Example-5 PRIVATE  PkgConfig::Clipper2Lib nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp Matplot++::matplot)


# Copy the example grid used for testing into the build dir
file(COPY auxiliaryFiles/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/auxiliaryFiles)
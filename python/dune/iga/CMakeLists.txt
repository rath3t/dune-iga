# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_subdirectory(basis)

add_python_targets(iga __init__ _igagrids _nurbspatchdata _boundarypatch _nurbsalgorithms generator)

dune_add_pybind11_module(NAME _iga)
set_property(
  TARGET _iga
  PROPERTY LINK_LIBRARIES dunecommon dunegeometry Clipper2Z nlohmann_json::nlohmann_json
           earcut_hpp::earcut_hpp
  APPEND
)
    target_link_libraries(
            _iga PRIVATE Clipper2Z nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp
    )

if(SKBUILD)
  install(TARGETS _iga LIBRARY DESTINATION python/dune/iga)
endif()

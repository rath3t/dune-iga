# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

install(FILES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga)

add_subdirectory(geometrykernel)
add_subdirectory(hierarchicpatch)
add_subdirectory(io)
add_subdirectory(splines)
add_subdirectory(parameterspace)
add_subdirectory(integrationrules)
add_subdirectory(utils)

add_subdirectory(test)

set(HEADERS nurbsbasis.hh patchgrid.hh)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/)

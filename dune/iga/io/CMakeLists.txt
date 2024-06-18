# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_subdirectory(vtk)

set(HEADERS griddrawer.hh ibrageometry.hh ibrajsonreader.hh ibrareader.hh createunstructuredgrid.hh)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/io)

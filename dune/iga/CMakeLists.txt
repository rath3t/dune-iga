# SPDX-FileCopyrightText: 2022 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-FileCopyrightText: 2022 The dune-iga
# developers mueller@ibb.uni-stuttgart.de
#
# SPDX-License-Identifier: LGPL-2.1-or-later

# install headers
file(
  GLOB_RECURSE headers
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_DEPENDS *.hh)
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga)

add_subdirectory(test)

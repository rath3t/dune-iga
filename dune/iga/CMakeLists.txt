# SPDX-FileCopyrightText: 2023 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
#install(
#  FILES
#  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga)

add_subdirectory(test)
#add_subdirectory(geometry)
add_subdirectory(io)
add_subdirectory(trim)
add_subdirectory(utils)
add_subdirectory(hierarchicpatch)

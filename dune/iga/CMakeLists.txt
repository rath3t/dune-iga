# SPDX-FileCopyrightText: 2023 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES bsplinealgorithms.hh
        controlpoint.hh
        gridcapabilities.hh
        nurbsalgorithms.hh
        nurbsbasis.hh
        nurbsgeometry.hh
        nurbsgrid.hh
        nurbsgridentity.hh
        nurbsgridindexsets.hh
        nurbsgridleafiterator.hh
        nurbsgridtraits.hh
        nurbsidset.hh
        nurbsintersection.hh
        nurbsleafgridview.hh
        nurbslocalgeometry.hh
        nurbspatch.hh
        nurbspatchdata.hh
        nurbspatchgeometry.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga)

add_subdirectory(test)
add_subdirectory(geometry)
add_subdirectory(io)
add_subdirectory(trim)
add_subdirectory(utils)

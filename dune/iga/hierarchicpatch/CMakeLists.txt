# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

set(HEADERS
    concepts.hh
    gridcapabilities.hh
    patchgrid.hh
    patchgridentity.hh
    patchgridfactory.hh
    patchgridfwd.hh
    patchgridgeometry.hh
    patchgridintersections.hh
    patchgridview.hh
    traits.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/hierarchicpatch)

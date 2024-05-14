# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

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

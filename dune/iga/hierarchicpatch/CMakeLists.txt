# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
        patchgridentity.hh
        patchgridentityseed.hh
        patchgridhierarchiciterator.hh
        patchgridgeometry.hh
        ../trimmer/identitytrimmer/patchgridlocalgeometry.hh
        ../trimmer/defaulttrimmer/patchgridintersectioniterator.hh
        patchgridintersections.hh
        patchgridindexsets.hh
        patchgridleafiterator.hh
        ../trimmer/identitytrimmer/patchgridleveliterator.hh
        patchgridview.hh
)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/hierarchicpatch)

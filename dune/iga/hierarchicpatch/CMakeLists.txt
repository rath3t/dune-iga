# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
        hierachicpatchgridentity.hh
        hierachicpatchgridentityseed.hh
  hierachicpatchgridhierarchiciterator.hh
  hierachicpatchgridgeometry.hh
        hierachicpatchgridlocalgeometry.hh
        geometrykernel/nurbspatchgeometrylocalview.hh
  hierachicpatchgridintersectioniterator.hh
  hierachicpatchgridintersections.hh
  hierachicpatchgridindexsets.hh
  hierachicpatchgridleafiterator.hh
  hierachicpatchgridleveliterator.hh
)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/hierarchicpatch)

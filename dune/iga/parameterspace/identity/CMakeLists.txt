# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

set(HEADERS
    idset.hh
    patchgridentityseed.hh
    patchgridhierarchiciterator.hh
    patchgridindexsets.hh
    patchgridintersectioniterator.hh
    patchgridleafiterator.hh
    patchgridleveliterator.hh
    patchgridlocalgeometry.hh
    parameterspace.hh
)

install(FILES ${HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/parameterspace/identitytrimmer
)

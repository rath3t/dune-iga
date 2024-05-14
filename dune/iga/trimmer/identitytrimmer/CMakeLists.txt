# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers

set(HEADERS
        idSet.hh
        patchgridentityseed.hh
        patchgridhierarchiciterator.hh
        patchgridindexsets.hh
        patchgridintersectioniterator.hh
        patchgridleafiterator.hh
        patchgridleveliterator.hh
        patchgridlocalgeometry.hh
        trimmer.hh
)

install(FILES ${HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/trimmer/identitytrimmer
)

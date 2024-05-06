# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_subdirectory(trimmingutils)

# install headers

set(HEADERS
        createlevel.hh
        elementtrimdata.hh
        entitycontainer.hh
        idset.hh
        patchgridhierarchiciterator.hh
        patchgridindexsets.hh
        patchgridintersectioniterator.hh
        patchgridintersections.hh
        patchgridleafiterator.hh
        patchgridleveliterator.hh
        patchtrimdata.hh
        referenceelement.hh
        trimmedentity.hh
        trimmedlocalgeometry.hh
        trimmer.hh
)

install(FILES ${HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/trimmer/defaulttrimmer
)

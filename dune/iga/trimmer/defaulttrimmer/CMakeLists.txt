# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_subdirectory(trimmingutils)
add_subdirectory(integrationrules)

# install headers

set(HEADERS
        createentities.hh
        createlevel.hh
        elementtrimdata.hh
        entitycontainer.hh
        entityinfo.hh
        idset.hh
        patchgridentityseed.hh
        patchgridhierarchiciterator.hh
        patchgridindexsets.hh
        patchgridintersectioniterator.hh
        patchgridintersections.hh
        patchgridleafiterator.hh
        patchgridleveliterator.hh
        patchtrimdata.hh
        referenceelement.hh
        trimelement.hh
        trimmedentity.hh
        trimmedlocalgeometry.hh
        trimmer.hh
        trimmerpreferences.hh
)

install(FILES ${HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/trimmer/defaulttrimmer
)

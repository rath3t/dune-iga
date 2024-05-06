# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_subdirectory(defaulttrimmer)
add_subdirectory(identitytrimmer)


# install headers
set(HEADERS
        concepts.hh
        localgeometryvariant.hh
)

install(FILES ${HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/trimmer
)


# SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(modules
    DuneIgaMacros.cmake
    AddAutoDiffFlags.cmake
    AddClipperLibFlags.cmake
    AddEarCutFlags.cmake
    AddEigenFlags.cmake
    AddnLohmannJsonFlags.cmake
)

install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})

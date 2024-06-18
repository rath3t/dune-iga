# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS ../../dune/**/*.hh)
list(APPEND HEADERS modules.txt mainpage.txt)
add_doxygen_target(DEPENDS ${HEADERS})

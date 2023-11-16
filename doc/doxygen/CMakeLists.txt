# SPDX-FileCopyrightText: 2023 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS ../../dune/**/*.hh)
message(STATUS ${HEADERS})
add_doxygen_target(DEPENDS ${HEADERS})

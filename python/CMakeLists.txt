# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
set(CXX_MAX_STANDARD 20)

add_subdirectory(dune)

dune_python_configure_bindings(
  PATH "." CMAKE_METADATA_FLAGS DUNE_OPTS_FILE HAVE_DUNE_ALUGRID HAVE_CLIPPERLIB2Z CXX_MAX_STANDARD
)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()

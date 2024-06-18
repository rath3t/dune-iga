# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_subdirectory(iga)

if(DUNE_ENABLE_PYTHONBINDINGS)
  add_subdirectory(python)
endif()

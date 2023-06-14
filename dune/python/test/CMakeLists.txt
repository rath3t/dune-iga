# SPDX-FileCopyrightText: 2023 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later

configure_file(setpath.py.in ${CMAKE_CURRENT_SOURCE_DIR}/setpath.py)

dune_python_add_test(
  NAME
  pytest1
  SCRIPT
  readGrid.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python)
dune_python_add_test(
  NAME
  pytest2
  SCRIPT
  poisson.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python)

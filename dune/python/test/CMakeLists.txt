# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

dune_python_add_test(
  NAME
  testidentitygrid
  SCRIPT
  testidentitygrid.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python
)


dune_python_add_test(
  NAME
  testtrimmedgrid
  SCRIPT
  testtrimmedgrid.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python
)

dune_python_add_test(
  NAME
  poissontest
  SCRIPT
  poissontest.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python
)

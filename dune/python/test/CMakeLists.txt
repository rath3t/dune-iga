# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

configure_file(setpath.py.in ${CMAKE_CURRENT_SOURCE_DIR}/setpath.py)

dune_python_add_test(NAME pytest1
                     SCRIPT readGrid.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python)
#dune_python_add_test(NAME pytest2
#                     SCRIPT test2.py
#                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                     LABELS quick)

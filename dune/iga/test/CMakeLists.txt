# SPDX-FileCopyrightText: 2022 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-2.1-or-later

set(TESTS dune-iga)
option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
foreach(_test ${TESTS})
  dune_add_test(SOURCES ${_test}.cc TIMEOUT 10000)
  target_compile_options(${_test} PRIVATE -g)
  if(ENABLE_TEST_COVERAGE)
    message("Enable Test Coverage")
    target_compile_options(${_test} PUBLIC -O0 -g -fprofile-arcs
                                           -ftest-coverage)
    target_link_options(${_test} PUBLIC -fprofile-arcs -ftest-coverage)

  endif()
endforeach()

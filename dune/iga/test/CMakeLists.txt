# SPDX-FileCopyrightText: 2023 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTS gridTests trimmedGridTests)
option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
foreach(_test ${TESTS})
  dune_add_test(SOURCES ${_test}.cc TIMEOUT 10000 LABELS cpp)
  target_compile_options(${_test} PRIVATE -g)

  # find_package(Matplot++ REQUIRED)

  target_link_libraries(
    ${_test} PRIVATE PkgConfig::Clipper2Lib nlohmann_json::nlohmann_json
                     earcut_hpp::earcut_hpp)

  target_compile_options(
    ${_test}
    PUBLIC -Wno-unused-variable -Wno-unused-parameter
           -Wno-unused-but-set-variable -Wno-deprecated-copy
           -ftemplate-backtrace-limit=0
    # -fsanitize=address
  )
  # target_link_options(${_test} PUBLIC -fsanitize=address -static-libasan)

  set_source_files_properties(${_test}.cc COMPILE_FLAGS -fno-var-tracking)

  if(ENABLE_TEST_COVERAGE)
    message("Enable Test Coverage")
    target_compile_options(${_test} PUBLIC -O0 -g -fprofile-arcs
                                           -ftest-coverage)
    target_link_options(${_test} PUBLIC -fprofile-arcs -ftest-coverage)

  endif()
endforeach()

# Copy the example grid used for testing into the build dir
file(COPY auxiliaryFiles/
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/auxiliaryFiles)

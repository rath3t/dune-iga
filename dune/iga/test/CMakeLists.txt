# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later


set(TESTS
        gridtests
        gridteststrimmed
        testbsplines
        testnurbsbasis
        # factorytest
        trimtest
        testfindintersection
        testibrareader
        testpatchgeometry
        testutillities
        testintersections
        testintegrationrules
        testvtkdatacollector
        testnurbsbasistrimmed
)
option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)

set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Threads REQUIRED)

foreach (_test ${TESTS})
    dune_add_test(SOURCES ${_test}.cpp TIMEOUT 10000 LABELS cpp)
    # add_dune_PkgConfig_Clipper2Lib_flags(${_test})
    target_compile_options(${_test} PRIVATE -g)
    target_link_libraries(${_test} PUBLIC Threads::Threads)

    target_link_libraries(
            ${_test} PRIVATE Clipper2Z nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp
            Matplot++::matplot
    )

    target_compile_options(
            ${_test} PUBLIC -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable
            -Wno-deprecated-copy -ftemplate-backtrace-limit=0 # -fsanitize=address
    )
    # target_link_options(${_test} PUBLIC -fsanitize=address -static-libasan)

    set_source_files_properties(${_test}.cc COMPILE_FLAGS -fno-var-tracking)

    if (ENABLE_TEST_COVERAGE)
        message("Enable Test Coverage")
        target_compile_options(${_test} PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
        target_link_options(${_test} PUBLIC -fprofile-arcs -ftest-coverage)

    endif ()
endforeach ()

dune_add_test(SOURCES concepttests.cpp COMPILE_ONLY LABELS cpp)

target_link_libraries(
        concepttests PRIVATE Clipper2Z nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp
)

target_compile_options(
        concepttests PUBLIC -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable
        -Wno-deprecated-copy -ftemplate-backtrace-limit=0
)

# Copy the example grid used for testing into the build dir
file(COPY auxiliaryfiles/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/auxiliaryfiles)

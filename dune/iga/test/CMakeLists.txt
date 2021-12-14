set(TESTS
  dune-iga
        staticlinearalgebratest
)

foreach(_test ${TESTS})
  dune_add_test(SOURCES ${_test}.cc)
  target_compile_options(${_test} PRIVATE -g -fconcepts-diagnostics-depth=5)
endforeach()

find_package(Eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})
target_link_libraries(staticlinearalgebratest PRIVATE Eigen3)

if( Eigen3_FOUND )
  message("Eigen was found")
  else()
  message("Eigen was not found")
  endif()
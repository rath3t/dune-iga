set(TESTS
  dune-iga
)

foreach(_test ${TESTS})
  dune_add_test(SOURCES ${_test}.cc)
  target_compile_options(${_test} PRIVATE -g -fconcepts-diagnostics-depth=5)
endforeach()


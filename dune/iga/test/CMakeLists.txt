set(TESTS
  dune-iga
  dune-iga_trimming_test
        examples
)

foreach(_test ${TESTS})
  dune_add_test(SOURCES ${_test}.cc)
  target_compile_options(${_test} PRIVATE -g -fconcepts-diagnostics-depth=5)
  target_link_libraries(${_test} PUBLIC nlohmann_json::nlohmann_json )
  target_link_libraries(${_test} PUBLIC pybind11::module pybind11::lto pybind11::windows_extras )
endforeach()


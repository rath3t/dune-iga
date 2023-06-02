# SPDX-FileCopyrightText: 2023 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later
add_definitions(-DDUNE_LOCALFEFUNCTIONS_USE_EIGEN=1)
set(CMAKE_CXX_STANDARD 20)

file(GLOB programSourceFiles CONFIGURE_DEPENDS *.cc)
add_custom_target(AllExamples) # Create target that includes all examples
foreach(programSourceFile ${programSourceFiles})
  get_filename_component(programName ${programSourceFile} NAME_WLE)
  add_executable(${programName} ${programSourceFile})
  set_target_properties(${programName} PROPERTIES CXX_STANDARD 20)
  add_dune_pythonlibs_flags(${programName})
  target_compile_options(
    ${programName}
    PUBLIC -Wno-unused-variable -Wno-unused-parameter
           -Wno-unused-but-set-variable -Wno-deprecated-copy
           -ftemplate-backtrace-limit=0)

  target_link_libraries(${programName} PUBLIC ikarus)
  target_link_libraries(
    ${programName} PUBLIC Eigen3::Eigen PkgConfig::Clipper2Lib
                          nlohmann_json::nlohmann_json earcut_hpp::earcut_hpp)
  add_dependencies(AllExamples ${programName})
  find_package(PythonLibs REQUIRED)
  include_directories(${PYTHON_INCLUDE_DIRS})
  target_link_libraries(${programName} PUBLIC ${PYTHON_LIBRARIES})
endforeach()

# Copy the example grid used for testing into the build dir
file(COPY auxiliaryFiles/
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/auxiliaryFiles)

cmake_minimum_required(VERSION 3.16)
project(dune-iga CXX)

set (CMAKE_CXX_STANDARD 20)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})


include(FetchContent)
FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG        v3.10.5
)
#
FetchContent_MakeAvailable(nlohmann_json)
include_directories("${PROJECT_SOURCE_DIR}/external_libraries")

FetchContent_Declare(
        pybind
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.9
)
#
FetchContent_MakeAvailable(pybind)

#include the dune macros
include(DuneMacros)
# start a dune project with information from dune.module

dune_project()
dune_enable_all_packages()
#target_link_libraries(${ProjectName} PRIVATE nlohmann_json::nlohmann_json)

add_subdirectory("src")
add_subdirectory("dune")
add_subdirectory("doc")
#add_subdirectory("cmake/modules")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

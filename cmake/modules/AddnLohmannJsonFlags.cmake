# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_SPDLOG for config.h
set(HAVE_NLOHMANNJSON ${nlohmann_json_FOUND})

# register all spdlog related flags
if(nlohmann_json_FOUND)
  dune_register_package_flags(LIBRARIES nlohmann_json::nlohmann_json COMPILE_DEFINITIONS "ENABLE_NLOHMANNJSON=1")
endif()

# add function to link against the spdlog library
function(add_dune_nlohmann_json_flags _targets)
  if(nlohmann_json_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC nlohmann_json::nlohmann_json)
      target_compile_definitions(${_target} PUBLIC ENABLE_NLOHMANNJSON=1)
    endforeach(_target)
  endif()
endfunction(add_dune_nlohmann_json_flags)

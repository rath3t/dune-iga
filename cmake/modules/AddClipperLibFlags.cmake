# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_SPDLOG for config.h
set(HAVE_CLIPPERLIB2 ${Clipper2_FOUND})

# register all spdlog related flags
if(Clipper2_FOUND)
  dune_register_package_flags(LIBRARIES PkgConfig::Clipper2Lib COMPILE_DEFINITIONS "ENABLE_CLIPPERLIB2=1")
endif()

# add function to link against the spdlog library
function(add_dune_Clipper2_flags _targets)
  if(Clipper2_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC PkgConfig::Clipper2Lib)
      target_compile_definitions(${_target} PUBLIC ENABLE_CLIPPERLIB2=1)
    endforeach(_target)
  endif()
endfunction(add_dune_Clipper2_flags)

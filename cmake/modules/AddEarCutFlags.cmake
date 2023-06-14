# SPDX-FileCopyrightText: 2023 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_EARCUT for config.h
set(HAVE_EARCUT ${earcut_hpp_FOUND})

# register all earcut_hpp related flags
if(earcut_hpp_FOUND)
  dune_register_package_flags(LIBRARIES earcut_hpp::earcut_hpp
                              COMPILE_DEFINITIONS "ENABLE_EARCUT=1")
endif()

# add function to link against the earcut_hpp library
function(add_dune_earcut_hpp_flags _targets)
  if(earcut_hpp_FOUND)
    message(STATUS "earcut_hpp_FOUND ${earcut_hpp_FOUND}")
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC earcut_hpp::earcut_hpp)
      target_compile_definitions(${_target} PUBLIC ENABLE_EARCUT=1)
    endforeach(_target)
  endif()
endfunction(add_dune_earcut_hpp_flags)

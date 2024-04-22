# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_CLIPPERLIB2 for config.h set(HAVE_CLIPPERLIB2Z ${Clipper2ZLib_FOUND})

# register all Clipper2 related flags
if(Clipper2Lib_FOUND)
  dune_register_package_flags(LIBRARIES Clipper2Z COMPILE_DEFINITIONS "HAVE_CLIPPERLIB2Z=1")
endif()

# add function to link against the Clipper2 library
function(add_dune_clipper2z_flags _targets)
  if(Clipper2Lib_FOUND)
    message(STATUS "Clipper2_FOUND ${Clipper2_FOUND}")
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC Clipper2Z)
      target_compile_definitions(${_target} PRIVATE USINGZ HAVE_CLIPPERLIB2Z=1)
    endforeach(_target)
  endif()
endfunction(add_dune_clipper2z_flags)

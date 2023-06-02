# SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
#
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(dunelocalfeformatFormat)

if(ADD_FORMATTARGET)
  include(CPM.cmake)
  cpmaddpackage(
    NAME
    Format.cmake
    VERSION
    1.7.3
    GITHUB_REPOSITORY
    TheLartians/Format.cmake
    OPTIONS
    "FORMAT_SKIP_CMAKE NO"
    "CMAKE_FORMAT_EXCLUDE config.h.cmake")
endif()

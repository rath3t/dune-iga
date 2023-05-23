# SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
#
# SPDX-License-Identifier: LGPL-2.1-or-later

# File for module specific CMake tests.
find_package(spdlog)
include(AddSpdlogFlags)

find_package(Matplot++)
include(AddMatplotppFlags)

find_package(Eigen3 3.3.9 REQUIRED)
include(AddEigenFlags)

find_package(nlohmann_json 3.11.2 REQUIRED)
include(AddnLohmannJsonFlags)

find_package(PkgConfig REQUIRED)
pkg_check_modules(Clipper2Lib REQUIRED IMPORTED_TARGET Clipper2)
include(AddClipperLibFlags)


find_package(earcut_hpp REQUIRED)
include(AddEarCutFlags)

# SPDX-FileCopyrightText: 2022 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-2.1-or-later

add_executable("dune-iga_example" dune-iga_example.cc)

target_link_dune_default_libraries("dune-iga_example")

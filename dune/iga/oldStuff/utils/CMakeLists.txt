# SPDX-FileCopyrightText: 2023 The dune-iga developers
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES concepts.hh fillquadraturerule.hh igahelpers.hh
        ../hierarchicpatch/utils/linearalgebra.hh mdnet.hh
        ../../utils/typetraits.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/utils)

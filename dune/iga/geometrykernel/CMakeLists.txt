# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

set(HEADERS
    algorithms.hh
    basealgorithms.hh
    closestpointprojection.hh
    controlpoint.hh
    findintersection.hh
    geohelper.hh
    makecirculararc.hh
    makesurfaceofrevolution.hh
    nurbspatchgeometry.hh
    nurbspatchgeometrylocalview.hh
    nurbspatchtransform.hh
    slicecurve.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/geometrykernel)

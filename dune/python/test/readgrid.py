# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
import sys
import math


import os, logging

# enable DUNE_SAVE_BUILD to test various output options
os.environ["DUNE_LOG_LEVEL"] = "warning"
os.environ["DUNE_SAVE_BUILD"] = "console"


from dune.iga import (
    IGAGrid,
    ControlPointNet,
    ControlPoint,
    NurbsPatchData,
    makeCircularArc,
    makeSurfaceOfRevolution,
)
from dune.grid import structuredGrid

from dune.iga import reader as readeriga
from dune.iga.basis import defaultGlobalBasis, Power, Lagrange, Nurbs
from dune.common import FieldVector
from dune.grid import gridFunction

if __name__ == "__main__":
    # one dimensional test
    cp = ControlPoint((0, 0, 0), 1)
    cp2 = ControlPoint((0, 0, 3), 1)
    netC = (cp, cp2)
    net = ControlPointNet(netC)
    nurbsPatchData = NurbsPatchData(((0, 0, 1, 1)), net, (1))
    gridView = IGAGrid(nurbsPatchData)
    reader = (readeriga.json, "../../iga/test/auxiliaryfiles/element_trim.ibra")
    refinements = 5
    gridView = IGAGrid(reader, dimgrid=2, dimworld=2)
    gridView.hierarchicalGrid.globalRefine(refinements)
    fv = FieldVector(
        10 * [0.1]
    )  # this is here due to https://gitlab.dune-project.org/core/dune-common/-/issues/340

    globalBasis = defaultGlobalBasis(gridView, Power(Nurbs(), 2))

    vtkWriter = gridView.trimmedVtkWriter()
    print(type(vtkWriter))

    gf1 = gridView.function(
        lambda e, x: math.sin(
            math.pi * (e.geometry.toGlobal(x)[0] + e.geometry.toGlobal(x)[1])
        )
    )

    gridView.writeVTK("test", pointdata=[gf1])

    @gridFunction(gridView)
    def g(x):
        return [math.sin(2 * math.pi * x[0] * x[1]), x[0] * x[1]] * 5

    vtkWriter.addPointData(g, name="g")
    vtkWriter.write(name="TestGrid")

    cp = ControlPoint((1, 2), 3)
    cp2 = cp + cp

    assert cp.coords[0] == 1
    assert cp.coords[1] == 2
    assert cp.weight == 3

    netC = ((cp, cp2), (cp, cp))
    net = ControlPointNet(netC)
    print(net.get((0, 1)))
    print(net.get((1, 0)))
    net.set((1, 0), ControlPoint((1, 4), 3))
    assert net.get((1, 0)).coords[0] == 1
    assert net.get((1, 0)).coords[1] == 4

    nurbsPatchData = NurbsPatchData(((0.0, 0, 1, 1), (0, 0, 1, 1)), net, (1, 1))

    gridView = IGAGrid(nurbsPatchData)
    assert gridView.size(0) == 1
    assert gridView.size(1) == 4
    assert gridView.size(2) == 4
    globalBasis = defaultGlobalBasis(gridView, Power(Nurbs(), 2))
    assert len(net) * 2 == len(globalBasis)
    newnurbsPatchData = nurbsPatchData.degreeElevate(0, 3)
    # degree elevate changes nothing
    gridView = IGAGrid(newnurbsPatchData)
    assert gridView.size(0) == 1
    assert gridView.size(1) == 4
    assert gridView.size(2) == 4
    newnurbsPatchData = nurbsPatchData.knotRefinement((0.5,), 0)
    newnurbsPatchData = newnurbsPatchData.knotRefinement((0.5,), 1)
    gridView = IGAGrid(newnurbsPatchData)
    assert gridView.size(0) == 4
    assert gridView.size(2) == 9

    r = 1.0
    R = 2.0
    circleData = makeCircularArc(r)
    nurbsPatchDataArc = makeSurfaceOfRevolution(circleData, (R, 0, 0), (0, 1, 0), 360.0)
    gridView = IGAGrid(nurbsPatchDataArc)
    vtkWriter = gridView.vtkWriter(subsampling=4)
    vtkWriter.write(name="Torus")

    try:
        nurbsPatchDataArc = makeSurfaceOfRevolution(
            newnurbsPatchData, (R, 0, 0), (0, 1, 0), 360.0
        )
        raise Exception(f"makeSurfaceOfRevolution should have raised an Value error")
    except ValueError:
        pass

    reader = (readeriga.json, "../../iga/test/auxiliaryfiles/element.ibra")
    gridView = IGAGrid(reader, dimgrid=2, dimworld=2)

    assert gridView.size(0) == 1
    assert gridView.size(1) == 4
    assert gridView.size(2) == 4
    gridView.hierarchicalGrid.globalRefine(1)
    gridView = gridView.hierarchicalGrid.leafView
    assert gridView.size(0) == 4
    assert gridView.size(2) == 9
    gridView.hierarchicalGrid.globalRefine(1)
    gridView = gridView.hierarchicalGrid.leafView
    assert gridView.size(0) == 16
    assert gridView.size(2) == 25

    assert gridView.dimGrid == 2
    assert gridView.dimWorld == 2

    patchD = gridView.hierarchicalGrid.patchData()
    assert patchD.degree[0] == 1
    assert patchD.degree[1] == 1

    # read and refine
    inputParameter = dict(
        file_path="../../iga/test/auxiliaryfiles/element.ibra",
        reader=readeriga.json,
        degree_elevate=(1, 1),
    )
    gridView2 = IGAGrid(inputParameter, dimgrid=2, dimworld=2)
    # degree elevation shouldn't change anything
    assert gridView2.size(0) == 1
    assert gridView2.size(1) == 4
    assert gridView2.size(2) == 4

    patchD = gridView2.hierarchicalGrid.patchData()

    assert patchD.degree[0] == 2
    assert patchD.degree[1] == 2

    inputParameter = dict(
        file_path="../../iga/test/auxiliaryfiles/element.ibra",
        reader=readeriga.json,
        pre_knot_refine=(1, 1),
    )
    gridView3 = IGAGrid(inputParameter, dimgrid=2, dimworld=2)
    # degree elevation shouldn't change anything
    assert gridView3.size(0) == 4
    assert gridView3.size(2) == 9

    inputParameter = dict(
        file_path="../../iga/test/auxiliaryfiles/element.ibra",
        reader=readeriga.json,
    )
    gridView4 = IGAGrid(inputParameter, dimgrid=2, dimworld=2)
    gridView4.hierarchicalGrid.globalRefine(1)
    assert gridView4.size(0) == 4
    assert gridView4.size(2) == 9

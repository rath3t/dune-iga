# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception


import sys
import setpath
import math

from mpi4py import MPI
setpath.set_path()
print(sys.path)

# from dune.iga import gridReader

#DUNE_LOG_LEVEL
#DUNE_FORCE_BUILD
#DUNE_SAVE_BUILD
from dune.iga import IGAGrid,ControlPointNet,ControlPoint,NurbsPatchData
from dune.grid import structuredGrid
#from dune.functions import Power,Lagrange,defaultGlobalBasis
#import dune.functions
from dune.iga import reader as readeriga
from dune.iga.basis import defaultGlobalBasis,Power,Lagrange,Nurbs
from dune.common import FieldVector
if __name__ == "__main__":
    reader = (readeriga.json, "../../iga/test/auxiliaryFiles/element.ibra")
    refineMents=5
    gridView = IGAGrid(reader, dimgrid=2,dimworld=2)
    gridView.hierarchicalGrid.globalRefine(refineMents)
    #print(help(gridView))
    fv = FieldVector(10*[0.1])

    # nurbsBasis= gridView.preBasis()
    globalBasis = defaultGlobalBasis(
        gridView, Power(Nurbs(), 2)
     )

    # assert 2*(4**refineMents)==len(globalBasis)*2

    #vtkWriter = gridView.trimmedVtkWriter()

    gf1 = gridView.function(lambda e,x: math.sin(math.pi*(e.geometry.toGlobal(x)[0]+e.geometry.toGlobal(x)[1])))

    gridView.writeVTK("test",pointdata=[gf1],name="sinFunc")
    print(help(arg))
    vtkWriter.addCellData()
    vtkWriter.write("TestGrid")

    # basisLagrange12 = testData(
    #     grid, Power(Lagrange(order=1),2)
    # )

    #
    # basisLagrange12 = dune.functions.defaultGlobalBasis(
    #     grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    # )

    cp=ControlPoint((1,2),3)
    cp2=cp+cp

    assert cp.coords[0]==1
    assert cp.coords[1]==2
    assert cp.weight==3

    netC=((cp,cp2),(cp,cp))
    net=ControlPointNet(netC)
    print(net.get((0,1)))
    print(net.get((1,0)))
    net.set((1,0),ControlPoint((1,4),3))
    assert net.get((1,0)).coords[0]==1
    assert net.get((1,0)).coords[1]==4

    nurbsPatchData=NurbsPatchData(((0,0,1,1),(0,0,1,1)),net,(1,1))

    gridView=IGAGrid(nurbsPatchData)
    assert gridView.size(0)==1
    assert gridView.size(1)==4
    assert gridView.size(2)==4
    # basisLagrange12 = testData(
    #     gridView, (1,2)
    # )

    nurbsBasis= nurbsPatchData.asBasis()
    globalBasis = defaultGlobalBasis(
            gridView, Power(nurbsBasis, 2)
        )

    assert len(net)==len(globalBasis)*2
    print("YEEEES")

    if True:
        reader = (readeriga.json, "../../iga/test/auxiliaryFiles/element.ibra")
        gridView = IGAGrid(reader, dimgrid=2,dimworld=2)

        assert gridView.size(0)==1
        assert gridView.size(1)==4
        assert gridView.size(2)==4
        gridView.hierarchicalGrid.globalRefine(1)
        gridView = gridView.hierarchicalGrid.leafView
        assert gridView.size(0)==4
        assert gridView.size(2)==9
        gridView.hierarchicalGrid.globalRefine(1)
        gridView = gridView.hierarchicalGrid.leafView
        assert gridView.size(0)==16
        assert gridView.size(2)==25

        assert gridView.dimGrid==2
        assert gridView.dimWorld==2

        # read and refine
        inputParameter= dict(
            file_path="../../iga/test/auxiliaryFiles/element.ibra",
            reader= readeriga.json,
            elevate_degree= (1,1)    )
        gridView2 = IGAGrid(inputParameter, dimgrid=2,dimworld=2)
        # degree elevation shouldn't change anything
        assert gridView2.size(0)==1
        assert gridView2.size(1)==4
        assert gridView2.size(2)==4

        inputParameter= dict(
            file_path="../../iga/test/auxiliaryFiles/element.ibra",
            reader= readeriga.json,
            pre_knot_refine= (1,1)    )
        gridView3 = IGAGrid(inputParameter, dimgrid=2,dimworld=2)
        # degree elevation shouldn't change anything
        assert gridView3.size(0)==4
        assert gridView3.size(2)==9

        inputParameter= dict(
            file_path="../../iga/test/auxiliaryFiles/element.ibra",
            reader= readeriga.json,
            post_knot_refine= (1,1)    )
        gridView4 = IGAGrid(inputParameter, dimgrid=2,dimworld=2)
        # degree elevation shouldn't change anything
        assert gridView4.size(0)==4
        assert gridView4.size(2)==9

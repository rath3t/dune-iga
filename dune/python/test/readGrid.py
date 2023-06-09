# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception


import setpath


setpath.set_path()


# from dune.iga import gridReader

#DUNE_LOG_LEVEL
#DUNE_FORCE_BUILD
#DUNE_SAVE_BUILD
import mpi4py
from dune.iga import IgaGrid,ControlPoint,ControlPointNet,NurbsPatchData

from dune.functions import Power,Lagrange,defaultGlobalBasis
#import dune.functions
from dune.iga import reader as readeriga
# from dune.iga.basis import Nurbs,igaDefaultGlobalBasis,Power,Lagrange
if __name__ == "__main__":
    lowerLeft = []
    upperRight = []
    elements = []
    for i in range(2):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    # grid = dune.grid.structuredGrid(lowerLeft, upperRight, elements)
    #
    # basisLagrange12 = dune.functions.defaultGlobalBasis(
    #     grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    # )
    if True:
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

        gridView=IgaGrid(nurbsPatchData)
        basisLagrange12 = defaultGlobalBasis(
            gridView, Power(Lagrange(order=1), 2)
        )

    # nurbsBasis= nurbsPatchData.asBasis()
        # globalBasis = dune.functions.defaultGlobalBasis(
        #     grid, dune.functions.Power(nurbsBasis, 2)
        # )
        #
        # assert len(net)==len(globalBasis)*2


    if True:
        reader = (readeriga.json, "../../iga/test/auxiliaryFiles/element.ibra")
        gridView = IgaGrid(reader, dimgrid=2,dimworld=2)

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
        gridView2 = IgaGrid(inputParameter, dimgrid=2,dimworld=2)
        # degree elevation shouldn't change anything
        assert gridView2.size(0)==1
        assert gridView2.size(1)==4
        assert gridView2.size(2)==4

        inputParameter= dict(
            file_path="../../iga/test/auxiliaryFiles/element.ibra",
            reader= readeriga.json,
            pre_knot_refine= (1,1)    )
        gridView3 = IgaGrid(inputParameter, dimgrid=2,dimworld=2)
        # degree elevation shouldn't change anything
        assert gridView3.size(0)==4
        assert gridView3.size(2)==9

        inputParameter= dict(
            file_path="../../iga/test/auxiliaryFiles/element.ibra",
            reader= readeriga.json,
            post_knot_refine= (1,1)    )
        gridView4 = IgaGrid(inputParameter, dimgrid=2,dimworld=2)
        # degree elevation shouldn't change anything
        assert gridView4.size(0)==4
        assert gridView4.size(2)==9

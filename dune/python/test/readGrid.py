# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception


import setpath


setpath.set_path()


# from dune.iga import gridReader



from dune.iga import igaGrid
from dune.grid import reader
from dune.iga import reader as readeriga
if __name__ == "__main__":

    reader = (readeriga.json, "../../iga/test/auxiliaryFiles/element.ibra")
    gridView = igaGrid(reader, dimgrid=2,dimworld=2)

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
    gridView2 = igaGrid(inputParameter, dimgrid=2,dimworld=2)
    # degree elevation shouldn't change anything
    assert gridView2.size(0)==1
    assert gridView2.size(1)==4
    assert gridView2.size(2)==4

    inputParameter= dict(
        file_path="../../iga/test/auxiliaryFiles/element.ibra",
        reader= readeriga.json,
        pre_knot_refine= (1,1)    )
    gridView3 = igaGrid(inputParameter, dimgrid=2,dimworld=2)
    # degree elevation shouldn't change anything
    assert gridView3.size(0)==4
    assert gridView3.size(2)==9

    inputParameter= dict(
        file_path="../../iga/test/auxiliaryFiles/element.ibra",
        reader= readeriga.json,
        post_knot_refine= (1,1)    )
    gridView4 = igaGrid(inputParameter, dimgrid=2,dimworld=2)
    # degree elevation shouldn't change anything
    assert gridView4.size(0)==4
    assert gridView4.size(2)==9

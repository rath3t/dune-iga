# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception


import setpath


setpath.set_path()


from dune.iga import gridReader



from dune.iga import igaGrid
from dune.grid import reader
if __name__ == "__main__":
    print("__name__-1")
    # jRr= gridReader(gridDim=2,worldDim=2)
    print("BLA-1")
    # grid=jRr.read("../../iga/test/auxiliaryFiles/element.ibra")

    # gv=grid.leafView
    # print("BLA-2")
    # alternative syntax for the same result exploiting dune grid file reader
    reader = (reader.dgf, "../../iga/test/auxiliaryFiles/element.ibra")
    print("BLA-3")
    gridView = igaGrid(reader, dimgrid=2,dimworld=2)
    print(help(gridView))
    # print("BLA-2")
    print(help(gridView))
    # print(help(jRr.read))
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


    for e in gridView.elements:
        print(e.geometry.corners[0])

    assert gridView.dimGrid==2
    assert gridView.dimWorld==2

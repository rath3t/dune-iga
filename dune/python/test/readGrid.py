# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception


import setpath


setpath.set_path()


from dune.iga import gridReader



from dune.iga import igaGrid
from dune.grid import reader
if __name__ == "__main__":
    # jRr= gridReader(gridDim=2,worldDim=2)
    # grid=jRr.read("../../iga/test/auxiliaryFiles/element.ibra")

    # alternative syntax for the same result exploiting dune grid file reader
    reader = (reader.dgf, "../../iga/test/auxiliaryFiles/element.ibra")
    grid = igaGrid(reader, dimgrid=2,dimworld=2)


    assert grid.size(0)==1
    assert grid.size(1)==4
    assert grid.size(2)==4
    grid.globalRefine(1)
    assert grid.size(0)==4
    assert grid.size(2)==9
    print(help(grid))
    gv = grid.leafView

    for e in elements:
        print(e)

    print(help(gv))

    assert gv.dimGrid==2
    assert gv.dimWorld==2

    #
    # reader = (dune.grid.reader.dgf, codeVecFunc)
    # grid = igaGrid(reader, dimgrid=2,worldDim=2)
    #
    # assert grid.size(0)==1
    # assert grid.size(1)==4
    # assert grid.size(2)==4

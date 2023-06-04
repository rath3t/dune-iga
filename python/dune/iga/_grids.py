from __future__ import absolute_import, division, print_function, unicode_literals

import sys, os
import logging
logger = logging.getLogger(__name__)

from ._iga import reader

def igaGrid(constructor, dimgrid=None, dimworld=None):
    print("BLA0")
    """
    Create an ALUGrid instance.

    Note: This functions has to be called on all cores and the parameters passed should be the same.
          Otherwise unexpected behavior will occur.

    Parameters:
    -----------

        constructor  means of constructing the grid, i.e. a grid reader or a
                     dictionary holding macro grid information
        dimgrid      dimension of grid, i.e. 2 or 3
        dimworld     dimension of world, i.e. 2 or 3 and >= dimension

    Returns:
    --------

    An IGAGrid instance with given refinement (conforming or nonconforming) and element type (simplex or cube).
    """
    from dune.grid.grid_generator import module, getDimgrid
    print("BLA0")
    if not dimgrid:
        dimgrid = getDimgrid(constructor)

    if dimworld is None:
        dimworld = dimgrid

    typeName = "Dune::IGA::NURBSGrid< " + str(dimgrid) + ", " + str(dimworld) + ",double>"

    includes = ["dune/iga/nurbsgrid.hh"]
    print("BLA0")
    gridModule = module(includes, typeName)
    print("BLA")
    print(help(constructor))
    if type(constructor) is dict:
        readGrid = gridModule.reader(constructor)
    elif type(constructor) is tuple:
        if(len(constructor))==2:
            readGrid = gridModule.reader(dict(reader=constructor[0],file_path=constructor[1]))
        else:
            raise Exception("Only tuple of size two are allowed   (readeriga.json, filename)")
    else:
        readGrid = gridModule.reader(constructor)
    # print(help(readGrid))
    readGrid = gridModule.LeafGrid(readGrid)
    # print(help(readGrid))
    print("help(readGrid)")
    return readGrid


grid_registry = {
    "IGA"        : igaGrid,

}

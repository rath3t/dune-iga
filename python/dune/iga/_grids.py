from __future__ import absolute_import, division, print_function, unicode_literals

import sys, os
import logging
logger = logging.getLogger(__name__)



def igaGrid(constructor, dimgrid=None, dimworld=None):
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

    if not dimgrid:
        dimgrid = getDimgrid(constructor)

    if dimworld is None:
        dimworld = dimgrid

    typeName = "Dune::IGA::NURBSGrid< " + str(dimgrid) + ", " + str(dimworld) + ">"

    includes = ["dune/iga/nurbsgrid.hh", "dune/iga/io/ibra/ibrareader.hh","dune/iga/nurbsleafgridview.hh","dune/iga/nurbsgridleafiterator.hh"]
    gridModule = module(includes, typeName)

    gridView = gridModule.LeafGrid(gridModule.reader(constructor))

    return gridView


grid_registry = {
    "IGA"        : igaGrid,

}

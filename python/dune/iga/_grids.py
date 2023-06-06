from __future__ import absolute_import, division, print_function, unicode_literals

import sys, os
import logging
logger = logging.getLogger(__name__)

from ._iga import reader
from dune.functions import Tree
from dune.functions import defaultGlobalBasis as defaultGlobalBasisBase
from dune.functions import preBasisTypeName as preBasisTypeNameBase
from dune.functions import indexMergingStrategy
class Nurbs(Tree):
    def __init__(self, order, dimRange=1):
        Tree.__init__(self, "Nurbs")
        self.dimRange = dimRange

    def __repr__(self):
            return "Nurbs"


def preBasisTypeName(tree, gridViewTypeName):
    assert isinstance(tree, Tree)
    if isinstance(tree, Nurbs):
        scalarPreBasis = "Dune::Functions::NurbsPreBasis< " + gridViewTypeName + " , " + str(tree.order) + " >"
        if tree.dimRange != 1:
            IMS = indexMergingStrategy(False, "interleaved")
            return "Dune::Functions::PowerPreBasis< " + IMS + " , " + scalarPreBasis + " , " + str(tree.dimRange) + " >"
        else:
            return scalarPreBasis
    else
        return preBasisTypeNameBase(tree, gridViewTypeName)


def defaultGlobalBasis(gridView, tree):
    from dune.functions import load

    headers = ["powerbasis", "compositebasis", "lagrangebasis", "subspacebasis", "defaultglobalbasis"]

    includes = []
    includes += list(gridView.cppIncludes)
    includes += ["dune/functions/functionspacebases/" + h + ".hh" for h in headers]
    includes += ["dune/iga/nurbsbasis.hh"]

    typeName = "Dune::Functions::DefaultGlobalBasis< " + preBasisTypeName(tree, gridView.cppTypeName) + " >"

    return load(includes, typeName).GlobalBasis(gridView)




def igaGrid(constructor, dimgrid=None, dimworld=None):
    print("BLA0")
    """
    Create an IGAGrid instance.

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

    gridModule = module(includes, typeName)

    if type(constructor) is dict:
        readGrid = gridModule.reader(constructor)
    elif type(constructor) is tuple:
        if(len(constructor))==2:
            readGrid = gridModule.reader(dict(reader=constructor[0],file_path=constructor[1]))
        else:
            raise Exception("Only tuple of size two are allowed   (readeriga.json, filename)")
    else:
        readGrid = gridModule.reader(constructor)
    gridView = gridModule.LeafGrid(readGrid)
    return gridView


grid_registry = {
    "IGA"        : igaGrid,

}

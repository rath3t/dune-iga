# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later


def IGAGrid(constructor, dimgrid=None, dimworld=None):
    """
    Create an IGAGrid instance.

    Parameters:
    -----------

        constructor  means of constructing the grid, i.e. a grid reader or a
                     dictionary holding macro grid information
        dimgrid      dimension of grid, i.e. 2 or 3
        dimworld     dimension of world, i.e. 2 or 3 and >= dimension

    Returns:
    --------

    An IGAGrid instance
    """

    if hasattr(constructor, "patchDim"):
        dimgrid = constructor.patchDim

    if hasattr(constructor, "dimworld"):
        dimworld = constructor.dimworld

    if dimgrid == None and dimworld == None:
        raise Exception(
            "If you don't pass the patch data you have to pass dimgrid and dimworld"
        )

    typeName = (
        "Dune::IGA::NURBSGrid< " + str(dimgrid) + ", " + str(dimworld) + ",double>"
    )

    includes = ["dune/python/iga/grid.hh"]
    from dune.generator.generator import SimpleGenerator
    from dune.common.hashit import hashIt

    generator = SimpleGenerator("HierarchicalGrid", "Dune::Python::IGA")
    moduleName = "NURBSGrid_" + hashIt(typeName)
    kwargs = dict()
    kwargs["dynamicAttr"] = True
    kwargs["holder"] = "std::shared_ptr"

    gridModule = generator.load(
        includes=includes, typeName=typeName, moduleName=moduleName, **kwargs
    )

    if type(constructor) is dict:
        readGrid = gridModule.reader(constructor)
    elif type(constructor) is tuple:
        if (len(constructor)) == 2:
            readGrid = gridModule.reader(
                dict(reader=constructor[0], file_path=constructor[1])
            )
        else:
            raise Exception(
                "Only tuple of size two are allowed   (readeriga.json, filename)"
            )
    else:
        readGrid = gridModule.HierarchicalGrid(constructor)
    gridView = gridModule.LeafGrid(readGrid)
    return gridView


grid_registry = {
    "IGA": IGAGrid,
}

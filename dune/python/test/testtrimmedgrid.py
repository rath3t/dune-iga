import os

import dune.iga
from dune.iga import IGAGrid, IGAGridType, ControlPointNet, ControlPoint, NurbsPatchData
from dune.iga import reader as readeriga

if __name__ == "__main__":
    filedir = os.path.dirname(__file__)
    filename = os.path.join(filedir, f"auxiliaryfiles/element_trim.ibra")
    inputParameter = dict(
        file_path=filename,
        reader=readeriga.json,
        pre_knot_refine=(1, 1),
    )
    gridView_untrimmed = IGAGrid(
        inputParameter, dimgrid=2, dimworld=2, gridType=IGAGridType.Identity
    )

    gridView_trimmed = IGAGrid(
        inputParameter, dimgrid=2, dimworld=2, gridType=IGAGridType.Default
    )

    dune.iga.registerParameterSpacePreferences(targetAccuracy=0.001)

    ## test grids

    assert gridView_untrimmed.size(0) == 4
    assert gridView_untrimmed.size(2) == 9

    assert gridView_trimmed.size(0) == 4
    assert gridView_trimmed.size(2) == 10

    gridView_trimmed.hierarchicalGrid.globalRefine(1)

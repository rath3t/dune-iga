import os
import logging

os.environ["DUNE_LOG_LEVEL"] = "warning"
os.environ["DUNE_SAVE_BUILD"] = "console"
from dune.iga import IGAGrid, IGAGridType, ControlPointNet, ControlPoint, NurbsPatchData
from dune.iga import reader as readeriga

if __name__ == "__main__":

    inputParameter = dict(
        file_path="../../iga/test/auxiliaryfiles/element_trim.ibra",
        reader=readeriga.json,
        pre_knot_refine=(1, 1),
    )
    gridView_untrimmed = IGAGrid(
        inputParameter, dimgrid=2, dimworld=2, gridType=IGAGridType.Identity
    )

    gridView_trimmed = IGAGrid(
        inputParameter, dimgrid=2, dimworld=2, gridType=IGAGridType.Default
    )

    ## test grids

    assert gridView_untrimmed.size(0) == 4
    assert gridView_untrimmed.size(2) == 9

    assert gridView_trimmed.size(0) == 4
    assert gridView_trimmed.size(2) == 10



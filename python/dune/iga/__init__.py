# from ._grids import *
from ._grids import *
from ._iga import *
from ._nurbspatchdata import *
from .basis import *
registry = dict()
registry["grid"] = grid_registry

# registry["globalBasis"] = {
#      "IGA" : defaultGlobalBasis
#  }

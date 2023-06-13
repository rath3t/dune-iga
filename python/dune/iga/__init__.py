from ._iga import *
from ._nurbspatchdata import *
from ._igagrids import *
from .basis import *
registry = dict()
registry["grid"] = grid_registry
#
registry["globalBasis"] = {
     "IGA" : defaultGlobalBasis
 }

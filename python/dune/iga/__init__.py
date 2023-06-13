from ._iga import *
from ._nurbspatchdata import *
from ._igagrids import *
registry = dict()
registry["grid"] = grid_registry
#
# registry["globalBasis"] = {
#      "IGA" : igaDefaultGlobalBasis
#  }

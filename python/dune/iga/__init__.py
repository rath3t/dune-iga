# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
from ._iga import *
from ._nurbspatchdata import *
from ._igagrids import *
from .basis import *
from ._boundarypatch import *
from ._nurbsalgorithms import *

registry = dict()
registry["grid"] = grid_registry
#
registry["globalBasis"] = {"IGA": defaultGlobalBasis}

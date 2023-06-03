# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
import dune



def gridReader(gridDim,worldDim) :
    generator = SimpleGenerator("IbraReader", "Dune::Python")

    element_type = f"Dune::IGA::IbraReader<{gridDim},{worldDim}>"

    includes = []
    includes += ["dune/python/iga/ibrareader.hh"]
    moduleName = "IbraReader_" + hashIt(element_type)
    module = generator.load(
        includes=includes,
        typeName=element_type,
        moduleName=moduleName
    )
    return module.IbraReader()

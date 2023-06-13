# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
from dune.common.hashit import hashIt
from dune.common import FieldVector
from dune.generator.generator import SimpleGenerator


def ControlPoint(coords, weight=1):
    fv = FieldVector(coords)
    generator = SimpleGenerator("ControlPoint", "Dune::Python")

    element_type = f"Dune::IGA::ControlPoint<{fv.cppTypeName}>"

    includes = []
    includes += ["dune/python/iga/nurbspatchdata.hh"]
    moduleName = "ControlPoint_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.ControlPoint(fv, weight)


def ControlPointNet(controlPoints):
    generator = SimpleGenerator("MultiDimensionNet", "Dune::Python")

    element_type = f"Dune::IGA::MultiDimensionNet<{len(controlPoints)},{controlPoints[0][0].cppTypeName}>"

    includes = []
    includes += ["dune/python/iga/nurbspatchdata.hh"]
    moduleName = "NurbsPatchData_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.MultiDimensionNet(controlPoints)


def NurbsPatchData(knotSpans, controlPointNet, degree):
    generator = SimpleGenerator("NurbsPatchData", "Dune::Python")

    worldDim = len(controlPointNet.get((0, 0)).coords)
    dim = controlPointNet.netDim
    # scalarType=print(typestr(controlPointNet.get((0,0)).weight))
    element_type = f"Dune::IGA::NURBSPatchData<{worldDim},{dim},double>"

    includes = []
    includes += ["dune/python/iga/nurbspatchdata.hh"]
    moduleName = "NurbsPatchData_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.NurbsPatchData(knotSpans, controlPointNet, degree)

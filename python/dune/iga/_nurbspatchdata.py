# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
from dune.common.hashit import hashIt
from dune.common import FieldVector
from .generator import MySimpleGenerator


def ControlPoint(coords, weight=1):
    fv = FieldVector(coords)
    generator = MySimpleGenerator("ControlPoint", "Dune::Python")

    element_type = f"Dune::IGA::ControlPoint<{fv.cppTypeName}>"

    includes = []
    includes += ["dune/python/iga/nurbspatchdata.hh"]
    moduleName = "ControlPoint_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.ControlPoint(fv, weight)


def ControlPointNet(controlPoints):
    generator = MySimpleGenerator("MultiDimensionalNet", "Dune::Python")

    try:
        controlPointType= controlPoints[0][0][0].cppTypeName
        netDim=3
    except:
        try:
            controlPointType= controlPoints[0][0].cppTypeName
            netDim=2
        except:
            try:
                controlPointType= controlPoints[0].cppTypeName
                netDim=1
            except:
                raise Exception("Controlpoint type not deducable from list")
    element_type = f"Dune::IGA::MultiDimensionalNet<{netDim},{controlPointType}>"

    includes = []
    includes += ["dune/python/iga/nurbspatchdata.hh"]
    moduleName = "ControlPointNet_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.MultiDimensionalNet(controlPoints)


def NurbsPatchData(knotSpans, controlPointNet, degree):
    generator = MySimpleGenerator("NurbsPatchData", "Dune::Python")

    worldDim = len(controlPointNet.directGet( 0).coords)
    dim = controlPointNet.netDim
    element_type = f"Dune::IGA::NURBSPatchData<{dim},{worldDim},double>"

    includes = []
    includes += ["dune/python/iga/nurbspatchdata.hh"]
    moduleName = "NurbsPatchData_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    if isinstance(knotSpans,tuple):
        if isinstance(knotSpans[0],float) or  isinstance(knotSpans[0], int):
            knotSpans= list([list(knotSpans),])
    if isinstance(degree,int):
        degree=(degree,)
    return module.NurbsPatchData(knotSpans, controlPointNet, degree)


def NurbsPatchDataDefault(dim, worldDim):
    generator = MySimpleGenerator("NurbsPatchData", "Dune::Python")

    element_type = f"Dune::IGA::NURBSPatchData<{dim},{worldDim},double>"
    includes = []
    includes += ["dune/python/iga/nurbspatchdata.hh"]
    moduleName = "NurbsPatchData_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.NurbsPatchData()

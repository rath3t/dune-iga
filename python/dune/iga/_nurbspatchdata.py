from dune.common.hashit import hashIt
from dune.common import FieldVector
from dune.generator.generator import SimpleGenerator
def ControlPoint(coords,weight=1):
    fv= FieldVector(coords)
    generator = SimpleGenerator("ControlPoint", "Dune::Python")

    element_type = f"Dune::IGA::ControlPoint<{fv.cppTypeName}>"

    includes = []
    includes += ["dune/python/iga/nurbspatchdata.hh"]
    moduleName = "ControlPoint_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.ControlPoint(fv,weight)



def ControlPointNet(controlPoints):
    generator = SimpleGenerator("MultiDimensionNet", "Dune::Python")

    element_type = f"Dune::IGA::MultiDimensionNet<{len(controlPoints)},{controlPoints[0][0].cppTypeName}>"

    includes = []
    includes += ["dune/python/iga/nurbspatchdata.hh"]
    moduleName = "NurbsPatchData_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.ControlPointNet(controlPoints)






def NurbsPatchData(knotSpans,controlPointNet,degree):
    generator = MySimpleGenerator("NurbsPatchData", "Dune::Python")


    element_type = f"Dune::IGA::NURBSPatchData<{pbfName}>"

    includes = []
    includes += ["dune/iga/nurbspatchdata.hh"]
    moduleName = "NurbsPatchData_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    basis = defaultGlobalBasis(gv, tree)
    return module.NurbsPatchData(knotSpans,controlPointNet,degree)










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

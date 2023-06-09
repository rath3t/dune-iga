# from __future__ import absolute_import, division, print_function, unicode_literals
#
# from dune.functions.tree import Composite, DG, Lagrange, Power, Tree
#
# class Nurbs(Tree):
#     def __init__(self, dimRange=1):
#         Tree.__init__(self, "Nurbs")
#         self.dimRange = dimRange
#
#     def __repr__(self):
#         if self.dimRange == 1:
#             return "Nurbs"
#         else:
#             return "Nurbs^" + str(self.dimRange)
#
#
# duneFunctionsLayouts = {"lexicographic": "Lexicographic", "interleaved": "Interleaved"}
#
# def indexMergingStrategy(blocked, layout):
#     return "Dune::Functions::BasisBuilder::" + ("Blocked" if blocked else "Flat") + duneFunctionsLayouts[layout]
#
#
# def preBasisTypeName(tree, gridViewTypeName):
#     assert isinstance(tree, Tree)
#     if isinstance(tree, Lagrange):
#         scalarPreBasis = "Dune::Functions::LagrangePreBasis< " + gridViewTypeName + " , " + str(tree.order) + " >"
#         if tree.dimRange != 1:
#             IMS = indexMergingStrategy(False, "interleaved")
#             return "Dune::Functions::PowerPreBasis< " + IMS + " , " + scalarPreBasis + " , " + str(tree.dimRange) + " >"
#         else:
#             return scalarPreBasis
#     elif isinstance(tree, Nurbs):
#         scalarPreBasis = "Dune::Functions::NurbsPreBasis< " + gridViewTypeName + " >"
#         if tree.dimRange != 1:
#             IMS = indexMergingStrategy(False, "interleaved")
#             return "Dune::Functions::PowerPreBasis< " + IMS + " , " + scalarPreBasis + " , " + str(tree.dimRange) + " >"
#         else:
#             return scalarPreBasis
#     elif isinstance(tree, DG):
#         raise Exception(repr(tree) + " not supported by dune-functions.")
#     elif isinstance(tree, Composite):
#         IMS = indexMergingStrategy(tree.blocked, tree.layout)
#         ChildPreBases = " , ".join(preBasisTypeName(c, gridViewTypeName) for c in tree.children)
#         return "Dune::Functions::CompositePreBasis< " + IMS + " , " + ChildPreBases + " >"
#     elif isinstance(tree, Power):
#         IMS = indexMergingStrategy(tree.blocked, tree.layout)
#         ChildPreBasis = preBasisTypeName(tree.children[0], gridViewTypeName)
#         return "Dune::Functions::PowerPreBasis< " + IMS + " , " + ChildPreBasis + " , " + str(tree.exponent) + " >"
#     else:
#         raise Exception("Unknown type of tree: " + repr(tree))
#
#
#
#
# from dune.generator.generator import SimpleGenerator
# from dune.common.hashit import hashIt
# from dune.functions import load
# def defaultGlobalBasis(gridView, tree):
#
#     generator = SimpleGenerator("GlobalBasis", "Dune::Python")
#
#     headers = ["powerbasis", "compositebasis", "lagrangebasis", "subspacebasis", "defaultglobalbasis"]
#
#     includes =  []
#     #includes += list(gridView.cppIncludes)
#     includes += ["dune/functions/functionspacebases/" + h + ".hh" for h in headers]
#     includes += ["dune/iga/nurbsbasis.hh"]
#     includes += ["dune/iga/nurbsgrid.hh"]
#     # includes += ["dune/iga/functionspacebases/nurbsbasis.hh"]
#     # includes +=  ["dune/python/functions/globalbasis.hh"]
#     # print(preBasisTypeName(tree, gridView.cppTypeName))
#     print(includes)
#     element_type  = "Dune::Functions::DefaultGlobalBasis< " + preBasisTypeName(tree, gridView.cppTypeName) + " >"
#     print(element_type)
#     moduleName = "globalBasis_" + hashIt(element_type)
#
#     # module = generator.load(
#     #     includes=includes, typeName=element_type, moduleName=moduleName
#     # )
#
#     # return module.GlobalBasis(gridView)
#     return load( includes=includes, typeName=element_type).GlobalBasis(gridView)

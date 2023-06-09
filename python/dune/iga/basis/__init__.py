from __future__ import absolute_import, division, print_function, unicode_literals


class Tree(object):
    def __init__(self, name, children=None):
        self.name = name
        self.children = []
        if children is not None:
            assert(all(isinstance(c, Tree) for c in children))
            self.children = list(children)

    def __mul__(self, other):
        return Composite(self, other)

    def __pow__(self, p):
        return Power(self, p)


class Lagrange(Tree):
    def __init__(self, order, dimRange=1):
        Tree.__init__(self, "Lagrange")
        self.order = order
        self.dimRange = dimRange

    def __repr__(self):
        if self.dimRange == 1:
            return "Lagrange<" + str(self.order) + ">"
        else:
            return "Lagrange<" + str(self.order) + ">^" + str(self.dimRange)


class DG(Tree):
    def __init__(self, order, dimRange=1):
        Tree.__init__(self, "DG")
        self.order = order
        self.dimRange = dimRange

    def __repr__(self):
        if self.dimRange == 1:
            return "DG<" + str(self.order) + ">"
        else:
            return "DG<" + str(self.order) + ">^" + str(self.dimRange)


class Composite(Tree):
    def __init__(self, *args, **kwargs):
        assert len(args) > 0
        Tree.__init__(self, "Composite", args)
        self.blocked = kwargs.get("blocked", False)
        self.layout = kwargs.get("layout", "lexicographic")

    def __repr__(self):
        return "(" + " * ".join(repr(c) for c in self.children) + ")"


class Power(Tree):
    def __init__(self, children, exponent, **kwargs):
        assert children is not None
        Tree.__init__(self, "Power", [children])
        assert len(self.children) == 1
        self.exponent = exponent
        self.blocked = kwargs.get("blocked", False)
        self.layout = kwargs.get("layout", "lexicographic")

    def __repr__(self):
        if self.exponent == 1:
            return repr(self.children[0])
        else:
            return "[" + repr(self.children[0]) + "]^" + str(self.exponent)


def toFem(tree, dimGrid):
    assert isinstance(tree, Tree)
    if isinstance(tree, Lagrange):
        return "LagrangeDiscreteFunctionSpace< FunctionSpace< " + str(dimGrid) + ", " + str(tree.dimRange) + " >, " + str(tree.order) + " >"
    elif isinstance(tree, DG):
        return "DiscontinuousGalerkinSpace< FunctionSpace< " + str(dimGrid) + ", " + str(tree.dimRange) + " >, " + str(tree.order) + " >"
    elif isinstance(tree, Composite):
        return "TupleDiscreteFunctionSpace< " + ", ".join(toFem(c, dimGrid) for c in tree.children) + " >"
    elif isinstance(tree, Power):
        return "PowerDiscreteFunctionSpace< " + toFem(tree.children[0], dimGrid) + ", " + str(tree.exponent) + " >"
    else:
        raise Exception("Unknown type of tree: " + repr(tree))

class Nurbs(Tree):
    def __init__(self, dimRange=1):
        Tree.__init__(self, "Nurbs")
        self.dimRange = dimRange

    def __repr__(self):
        if self.dimRange == 1:
            return "Nurbs"
        else:
            return "Nurbs^" + str(self.dimRange)


duneFunctionsLayouts = {"lexicographic": "Lexicographic", "interleaved": "Interleaved"}

def indexMergingStrategy(blocked, layout):
    return "Dune::Functions::BasisBuilder::" + ("Blocked" if blocked else "Flat") + duneFunctionsLayouts[layout]


def preBasisTypeName(tree, gridViewTypeName):
    assert isinstance(tree, Tree)
    if isinstance(tree, Lagrange):
        scalarPreBasis = "Dune::Functions::LagrangePreBasis< " + gridViewTypeName + " , " + str(tree.order) + " >"
        if tree.dimRange != 1:
            IMS = indexMergingStrategy(False, "interleaved")
            return "Dune::Functions::PowerPreBasis< " + IMS + " , " + scalarPreBasis + " , " + str(tree.dimRange) + " >"
        else:
            return scalarPreBasis
    elif isinstance(tree, Nurbs):
        scalarPreBasis = "Dune::Functions::NurbsPreBasis< " + gridViewTypeName + " >"
        if tree.dimRange != 1:
            IMS = indexMergingStrategy(False, "interleaved")
            return "Dune::Functions::PowerPreBasis< " + IMS + " , " + scalarPreBasis + " , " + str(tree.dimRange) + " >"
        else:
            return scalarPreBasis
    elif isinstance(tree, DG):
        raise Exception(repr(tree) + " not supported by dune-functions.")
    elif isinstance(tree, Composite):
        IMS = indexMergingStrategy(tree.blocked, tree.layout)
        ChildPreBases = " , ".join(preBasisTypeName(c, gridViewTypeName) for c in tree.children)
        return "Dune::Functions::CompositePreBasis< " + IMS + " , " + ChildPreBases + " >"
    elif isinstance(tree, Power):
        IMS = indexMergingStrategy(tree.blocked, tree.layout)
        ChildPreBasis = preBasisTypeName(tree.children[0], gridViewTypeName)
        return "Dune::Functions::PowerPreBasis< " + IMS + " , " + ChildPreBasis + " , " + str(tree.exponent) + " >"
    else:
        raise Exception("Unknown type of tree: " + repr(tree))




from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
# from dune.functions import load
def igaDefaultGlobalBasis(gridView, tree):

    generator = SimpleGenerator("IGAGlobalBasis", "Dune::Python::IGA")

    headers = ["powerbasis", "compositebasis", "lagrangebasis", "subspacebasis", "defaultglobalbasis"]

    includes =  []
    includes += ["dune/python/iga/globalBasis.hh"]
    includes += ["dune/functions/functionspacebases/" + h + ".hh" for h in headers]
    # includes += ["dune/iga/nurbsbasis.hh"]
    # includes += ["dune/iga/nurbsgrid.hh"]
    # includes += ["dune/iga/functionspacebases/nurbsbasis.hh"]
    # includes +=  ["dune/python/functions/globalbasis.hh"]
    # print(preBasisTypeName(tree, gridView.cppTypeName))
    print(includes)
    element_type  = "Dune::Functions::DefaultGlobalBasis< " + preBasisTypeName(tree, gridView.cppTypeName) + " >"
    print(element_type)
    moduleName = "igaGlobalBasis_" + hashIt(element_type)
    kwargs=dict()
    kwargs["dynamicAttr"] = True
    kwargs["holder"] = "std::shared_ptr"
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName,**kwargs
    )

    return module.IGAGlobalBasis(gridView)
    # return load( includes=includes, typeName=element_type).GlobalBasis(gridView)

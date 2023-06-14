# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.functions import Power, Lagrange, DG, Composite, Tree


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
    return (
        "Dune::Functions::BasisBuilder::"
        + ("Blocked" if blocked else "Flat")
        + duneFunctionsLayouts[layout]
    )


def preBasisTypeName(tree, gridViewTypeName):
    assert isinstance(tree, Tree)
    if isinstance(tree, Lagrange):
        scalarPreBasis = (
            "Dune::Functions::LagrangePreBasis< "
            + gridViewTypeName
            + " , "
            + str(tree.order)
            + " >"
        )
        if tree.dimRange != 1:
            IMS = indexMergingStrategy(False, "interleaved")
            return (
                "Dune::Functions::PowerPreBasis< "
                + IMS
                + " , "
                + scalarPreBasis
                + " , "
                + str(tree.dimRange)
                + " >"
            )
        else:
            return scalarPreBasis
    elif isinstance(tree, Nurbs):
        scalarPreBasis = "Dune::Functions::NurbsPreBasis< " + gridViewTypeName + " >"
        if tree.dimRange != 1:
            IMS = indexMergingStrategy(False, "interleaved")
            return (
                "Dune::Functions::PowerPreBasis< "
                + IMS
                + " , "
                + scalarPreBasis
                + " , "
                + str(tree.dimRange)
                + " >"
            )
        else:
            return scalarPreBasis
    elif isinstance(tree, DG):
        raise Exception(repr(tree) + " not supported by dune-functions.")
    elif isinstance(tree, Composite):
        IMS = indexMergingStrategy(tree.blocked, tree.layout)
        ChildPreBases = " , ".join(
            preBasisTypeName(c, gridViewTypeName) for c in tree.children
        )
        return (
            "Dune::Functions::CompositePreBasis< " + IMS + " , " + ChildPreBases + " >"
        )
    elif isinstance(tree, Power):
        IMS = indexMergingStrategy(tree.blocked, tree.layout)
        ChildPreBasis = preBasisTypeName(tree.children[0], gridViewTypeName)
        return (
            "Dune::Functions::PowerPreBasis< "
            + IMS
            + " , "
            + ChildPreBasis
            + " , "
            + str(tree.exponent)
            + " >"
        )
    else:
        raise Exception("Unknown type of tree: " + repr(tree))


from dune.functions import load


def defaultGlobalBasis(gridView, tree):
    headers = [
        "powerbasis",
        "compositebasis",
        "lagrangebasis",
        "subspacebasis",
        "defaultglobalbasis",
    ]

    includes = []
    includes += ["dune/functions/functionspacebases/" + h + ".hh" for h in headers]
    includes += list(gridView.cppIncludes)
    includes += ["dune/iga/nurbsbasis.hh"]
    typeName = (
        "Dune::Functions::DefaultGlobalBasis< "
        + preBasisTypeName(tree, gridView.cppTypeName)
        + " >"
    )

    return load(includes, typeName).GlobalBasis(gridView)

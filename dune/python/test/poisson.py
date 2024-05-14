# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
# The top part of this code is taken from https://gitlab.dune-project.org/staging/dune-functions/-/blob/master/dune/python/test/poisson.py
# Thus, see also their license https://gitlab.dune-project.org/staging/dune-functions/-/blob/master/COPYING
import setpath

setpath.set_path()
import numpy as np
import math
from scipy.sparse import lil_matrix
import scipy as sp
import dune.geometry
import dune.grid
from dune.iga import reader as readeriga
from dune.iga.basis import defaultGlobalBasis, Power, Lagrange, Nurbs
from dune.iga import IGAGrid, boundaryPatch
from dune.grid import gridFunction
from dune.common import FieldVector
import os

os.environ["ALUGRID_VERBOSITY_LEVEL"] = "0"
os.environ['DUNE_LOG_LEVEL'] = 'debug'
os.environ['DUNE_SAVE_BUILD'] = 'console'

# f(x) = 1
f = lambda x: 1

# g(x) = 0
g = lambda x: 0


# boundary of the unit square
def isDirichlet(x):
    return x[0] <= 0.0 or x[0] >= 1.0 or x[1] >= 1.0 or x[1] <= 0.0


# TODO: This assembler loop is very inefficient in terms of run time and should be improved using Python vectorization.
# See discussion at https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/295 for hints and code pointers.
def localAssembler(element, localBasis):
    n = len(localBasis)

    localA = np.zeros((n, n))
    localB = np.zeros(n)

    # choose a high enough quadrature order
    quadOrder = 3

    # create a quadrature rule and integrate
    quadRule = dune.geometry.quadratureRule(element.type, quadOrder)
    for pt in quadRule:
        pos = pt.position

        # evaluate the local basis functions (this is an array!)
        phi = localBasis.evaluateFunction(pos)

        # evaluate the local basis reference Jacobians (array of arrays)
        phiRefGradient = localBasis.evaluateJacobian(pos)

        # get the transformation matrix
        jacobianInverseTransposed = element.geometry.jacobianInverseTransposed(pos)

        # |det(jacobianInverseTransposed)|
        integrationElement = element.geometry.integrationElement(pos)

        # transform the reference Jacobians to global geometry
        phiGradient = [
            np.dot(jacobianInverseTransposed, np.array(g)[0]) for g in phiRefGradient
        ]

        posGlobal = element.geometry.toGlobal(pos)

        for i in range(n):
            for j in range(n):
                # compute grad(phi_i) * grad(phi_j)
                localA[i, j] += (
                    pt.weight
                    * integrationElement
                    * np.dot(phiGradient[i], phiGradient[j])
                )

            localB[i] += pt.weight * integrationElement * phi[i] * f(posGlobal)

    return localA, localB


def globalAssembler(basis):
    N = len(basis)

    A = lil_matrix((N, N))
    b = np.zeros(N)

    # mark all Dirichlet DOFs (first sqrt(N) are fixed)
    dichichletDOFs = np.zeros(N, dtype=bool)
    dichichletDOFs[: math.floor(math.sqrt(N))] = True

    # interpolate the boundary values
    gCoeffs = np.zeros(N)
    # basis.interpolate(gCoeffs,g)

    # extract grid and localView
    localView = basis.localView()
    grid = basis.gridView

    for element in grid.elements:
        # assign the localView to the current element
        localView.bind(element)

        # set of all shape functions with support in this element
        localBasis = localView.tree().finiteElement.localBasis

        localN = len(localBasis)

        localA, localb = localAssembler(element, localBasis)

        # copy the local entries into the global matrix using the
        # index mapping given by the localView
        for i in range(localN):
            gi = localView.index(i)[0]

            if dichichletDOFs[gi]:
                A[gi, gi] = 1.0
                b[gi] = gCoeffs[gi]
            else:
                b[gi] += localb[i]
                for j in range(localN):
                    gj = localView.index(j)[0]
                    A[gi, gj] += localA[i, j]

    return A, b


############################### START ########################################

# create a trimmed grid from file
reader = (readeriga.json, "../../iga/test/auxiliaryfiles/element_trim_xb.ibra")

gridView = IGAGrid(reader, dimgrid=2, dimworld=2)
gridView.hierarchicalGrid.globalRefine(2)

neumannVertices = np.zeros(gridView.size(2), dtype=bool)

boundaryPatch = boundaryPatch(gridView, neumannVertices)

basis = defaultGlobalBasis(gridView, Nurbs())

# compute A and b
A, b = globalAssembler(basis)

# convert A to Compressed Sparse Column format
Acsc = A.tocsc()

# solve!
x = sp.sparse.linalg.spsolve(Acsc, b)

vtkWriter = gridView.trimmedVtkWriter()
discreteFunction = basis.asFunction(x)

vtkWriter.addPointData(discreteFunction, name="g")
vtkWriter.write(name="Poisson")

xTest = [
    0.0,
    0.0,
    0.0,
    0.0,
    0.60278885,
    0.70725846,
    2.04974976,
    2.00540431,
    1.86744155,
    1.45648652,
    1.31236311,
    2.04982725,
    2.01589274,
    1.91911023,
    1.81531395,
    1.79925452,
    2.04866857,
    2.01931155,
    1.9486634,
    1.88186179,
    1.85137355,
]

assert np.linalg.norm(x - xTest) < 1e-7

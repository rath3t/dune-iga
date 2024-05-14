# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.generator.algorithm import run
from dune.generator import path
from io import StringIO
from dune.generator.importclass import load
from dune.iga import NurbsPatchDataDefault

runCircularArc = """
#include <dune/iga/geometrykernel/makecirculararc.hh>
auto run( double radius = 1.0,  double startAngle = 0.0, double endAngle = 360.0,
                       const Dune::FieldVector<double, 3>& origin = {0, 0, 0},
                       const Dune::FieldVector<double, 3>& X      = {1, 0, 0},
                       const Dune::FieldVector<double, 3>& Y      = {0, 1, 0})
{
  return Dune::IGA::makeCircularArc(radius,startAngle,endAngle,origin,X,Y);
}
"""


def makeCircularArc(
    radius, startAngle=None, endAngle=None, origin=None, X=None, Y=None
):
    NurbsPatchDataDefault(dim=1, worldDim=3)
    if startAngle is None:
        return run("run", StringIO(runCircularArc), radius)
    elif endAngle is None:
        return run("run", StringIO(runCircularArc), radius, startAngle)
    elif origin is None:
        return run("run", StringIO(runCircularArc), radius, startAngle, endAngle)
    elif X is None:
        return run(
            "run", StringIO(runCircularArc), radius, startAngle, endAngle, origin
        )
    elif Y is None:
        return run(
            "run", StringIO(runCircularArc), radius, startAngle, endAngle, origin, X
        )
    else:
        return run(
            "run", StringIO(runCircularArc), radius, startAngle, endAngle, origin, X, Y
        )


runSurfaceOfRevolution = """
#include <dune/iga/geometrykernel/makesurfaceofrevolution.hh>
auto run(const Dune::IGA::NURBSPatchData<1, 3, double>& generatrix,
                               const Dune::FieldVector<double, 3>& point,
                               const Dune::FieldVector<double, 3>& revolutionaxisI,
                               const double& revolutionAngle = 360.0)
{
  return Dune::IGA::makeSurfaceOfRevolution(generatrix,point,revolutionaxisI,revolutionAngle);
}
"""


def makeSurfaceOfRevolution(generatrix, point, revolutionaxisI, revolutionAngle=None):
    if generatrix.patchDim != 1 and generatrix.dimworld != 3:
        raise ValueError(
            f"The passed nurbs patchdata has to be a curve in 3D. You passed a {generatrix.patchDim}-dimensional object in {generatrix.dimworld}D space"
        )
    from dune.common import FieldVector

    NurbsPatchDataDefault(dim=2, worldDim=3)
    pointFV = FieldVector(point)
    revolutionaxisIFV = FieldVector(revolutionaxisI)
    print(generatrix)
    if revolutionAngle is None:
        return run(
            "run",
            StringIO(runSurfaceOfRevolution),
            generatrix,
            pointFV,
            revolutionaxisIFV,
        )
    else:
        return run(
            "run",
            StringIO(runSurfaceOfRevolution),
            generatrix,
            pointFV,
            revolutionaxisIFV,
            revolutionAngle,
        )

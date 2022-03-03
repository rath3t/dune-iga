// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/common/test/testsuite.hh>
#include <dune/iga/kernel/ANurbs.h>

using namespace Dune;



void testTorusGeometry() {

  anurbs::NurbsSurfaceGeometry<3> surf(2,1,4,3, false);
  Eigen::VectorXd knotsu(5);
  knotsu <<0, 0, 7.5, 15, 15;
  surf.knots_u() = knotsu;

  Eigen::VectorXd knotsv(3);
  knotsu <<0, 10, 20;
  surf.knots_u() = knotsu;
  Eigen::MatrixX3d cps(12,3);
  cps <<
      -10.0, - 5.0, -1.0         ,
      -12.0,   3.0,  3.0         ,
      - 9.0,  11.0, -0.0701928417,
      - 5.0, - 3.0,  1.0         ,
      - 6.0,   4.0, -2.0         ,
      - 5.0,   7.0,  0.9298071583,
        0.0, - 4.0, -1.0         ,
        1.0,   6.0,  5.0         ,
        0.0,  13.0, -0.2350184214,
        4.0, - 2.0,  0.0         ,
        5.0,   4.0, -1.0         ,
        5.0,  11.0,  0.7649815786;
  surf.poles()= cps;

  Eigen::Vector3d vec =surf.point_at(12.0,5.0);
  Eigen::Vector3d vecExpected;
  vecExpected<< 1.46, 0.96, 0.9;

  TestSuite test;
  test.check(vec.isApprox(vecExpected));

}






int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);

    testTorusGeometry();

  return 0;
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
}

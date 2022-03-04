// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/iga/kernel/ANurbs.h>

using namespace Dune;

void testSurfaceGeometry() {
  anurbs::NurbsSurfaceGeometry<3> surf(2, 1, 4, 3, false);
  Eigen::VectorXd knotsu(5);
  knotsu << 0, 0, 7.5, 15, 15;
  surf.knots_u() = knotsu;

  Eigen::VectorXd knotsv(3);
  knotsv << 0, 10, 20;
  surf.knots_v() = knotsv;
  Eigen::MatrixX3d cps(12, 3);
  cps << -10.0, -5.0, -1.0, -12.0, 3.0, 3.0, -9.0, 11.0, -0.0701928417, -5.0, -3.0, 1.0, -6.0, 4.0, -2.0, -5.0, 7.0, 0.9298071583, 0.0,
      -4.0, -1.0, 1.0, 6.0, 5.0, 0.0, 13.0, -0.2350184214, 4.0, -2.0, 0.0, 5.0, 4.0, -1.0, 5.0, 11.0, 0.7649815786;
  surf.poles() = cps;

  Eigen::Vector3d vec = surf.point_at(12.0, 5.0);
  Eigen::Vector3d vecExpected;
  vecExpected << 1.46, 0.96, 0.9;

  TestSuite test;
  test.check(vec.isApprox(vecExpected), "The evaluated position is not correct!");
}

void testPointOnSurfaceProjection() {
  anurbs::Model model;

  model.load("../../../../dune/iga/test/data/point_on_surface_projection.ibra");

  auto surf       = model.template of_type<anurbs::NurbsSurfaceGeometry<3>>()[0].data();
  auto projection = anurbs::PointOnSurfaceProjection<3>(surf);
  Eigen::MatrixX2d test_data(100, 2);
  test_data << 0.2486685842, 0.1107439772, 0.9978910512, 0.875821514, 0.0830195928, 0.9259765381, 0.9496291005, 0.1502956968, 0.459878301,
      0.6225033205, 0.6025333179, 0.9861121243, 0.0790516166, 0.5709693747, 0.5655969403, 0.2574736384, 0.8305452102, 0.5624009979,
      0.2794794688, 0.8097236998, 0.3071577155, 0.3822328865, 0.99964868, 0.4057108985, 0.0189992101, 0.2765154416, 0.0056710457,
      0.0347536937, 0.7897584256, 0.8313703019, 0.7195281758, 0.1178902761, 0.4422092407, 0.0883118692, 0.6658555305, 0.4627695477,
      0.5076745905, 0.8196537363, 0.2832733436, 0.6378145738, 0.0645325822, 0.7572020105, 0.9879030846, 0.6847038878, 0.7960193063,
      0.3109410034, 0.6945010543, 0.6643853852, 0.8336644866, 0.0037464825, 0.3930552445, 0.9678212548, 0.4785460939, 0.3808379352,
      0.585713005, 0.0070576593, 0.3991026969, 0.2632817497, 0.2006969681, 0.2421369582, 0.1202849104, 0.4221172297, 0.7358854449,
      0.9698849781, 0.5188882311, 0.5018097253, 0.3654126839, 0.5201963845, 0.6403577298, 0.8467138102, 0.8726223264, 0.9882953218,
      0.0928173764, 0.1633056552, 0.155525548, 0.028770891, 0.0010793009, 0.456353415, 0.5686607461, 0.7125619164, 0.1982382905,
      0.7190583882, 0.3281335357, 0.0049134241, 0.8250263221, 0.6903648873, 0.2025081111, 0.9819264952, 0.2204824454, 0.5096595127,
      0.3834477823, 0.8546278755, 0.9842980774, 0.0390799376, 0.1078386875, 0.3153865423, 0.3720683969, 0.7303910846, 0.9586810097,
      0.5417624263, 0.6905896192, 0.2880027598, 0.5518418474, 0.1284805062, 0.7797412704, 0.4522544888, 0.9001546897, 0.7999264513,
      0.7019192277, 0.5565699151, 0.8618938363, 0.8969601173, 0.573701815, 0.3725593828, 0.5645628416, 0.6001849126, 0.9691815544,
      0.2581378097, 0.2934454783, 0.2084521582, 0.5663407671, 0.890895538, 0.8255418832, 0.1905749553, 0.9205142534, 0.4522520338,
      0.4010180167, 0.4240173988, 0.6520043102, 0.7641253019, 0.4829338433, 0.2314186642, 0.6469415439, 0.1985799857, 0.0100815608,
      0.8321075117, 0.2825933291, 0.9272068771, 0.2162422287, 0.4111952472, 0.1691770848, 0.6298445815, 0.9101820793, 0.3544982711,
      0.7672993754, 0.7465931879, 0.9772622413, 0.7692543193, 0.0181493168, 0.3759801157, 0.7127572883, 0.0072015384, 0.5138806088,
      0.9528929917, 0.8881267998, 0.633660035, 0.1209032033, 0.8468105364, 0.3423240503, 0.1053561797, 0.0060057966, 0.1696525231,
      0.8454230683, 0.0981951589, 0.1910913443, 0.8742888891, 0.9726746664, 0.9913883214, 0.714155909, 0.3804675184, 0.8380541009,
      0.3771821253, 0.4770630204, 0.6988195776, 0.1394817876, 0.4965714307, 0.2955455749, 0.7173407975, 0.3722430087, 0.3385346254,
      0.8603301383, 0.259464834, 0.0317961453, 0.6580847188, 0.0028120549, 0.9568088497, 0.4066441252, 0.1674200367, 0.6186872817,
      0.3124124684, 0.001451293, 0.5436467522, 0.6059298895, 0.5078425922, 0.6547791262, 0.0598361939, 0.1788401213, 0.1605233313,
      0.460314065, 0.8863533544;

  TestSuite test;
  for (int i = 0; auto row : test_data.rowwise()) {
    auto expected_point               = surf->point_at(row[0], row[1]);
    const auto [success, u, v, point] = projection.get(expected_point);
    test.check(success, "Projection didn't work.");
    test.check(Dune::FloatCmp::eq(u, row[0], 1e-9), "u coord doesn't match.");
    test.check(Dune::FloatCmp::eq(v, row[1], 1e-9), "v coord doesn't match.");
    test.check(point.isApprox(expected_point, 1e-10), "3d coordinate doesn't match.");
  }
}

void testPointOnCurveProjection() {
  auto curv = std::make_shared<anurbs::NurbsCurveGeometry<3>>(4, 8, false);
  Eigen::VectorXd knotsu(11);
  knotsu << 3, 3, 3, 3, 6.5, 10, 13.5, 17, 17, 17, 17;
  curv->knots() = knotsu;

  Eigen::MatrixX3d cps(8, 3);
  cps << 0, -25, -5, -15, -15, 0, 5, -5, -3, 15, -15, 3, 25, 0, 6, 15, 15, 6, -5, -5, -3, -25, 15, 4;
  curv->poles() = cps;

  TestSuite test;
  const double tol = 1e-7;
  auto projection = anurbs::PointOnCurveProjection<3>(curv, tol, 1e-8);

  projection.compute({-25.1331415843, -38.9256022819, -3.2989320128});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.39832829120118,tol));

  projection.compute({35.6464813397, 27.3703996918, -41.1153099924});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.3339477286703,tol));

  projection.compute({-40.3995502695, 45.1689836547, -1.7412051334});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({36.117074539, 44.5183648237, 47.2049699152});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.0827431544026,tol));

  projection.compute({36.8315563476, -48.7314244261, 46.3990433125});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 9.4334414007513,tol));

  projection.compute({-39.7935307537, 1.0082708909, -48.4975476742});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 16.7432603141062,tol));

  projection.compute({39.2152096095, -39.0656723124, -28.995046196});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({-11.5997738492, 6.2506795657, 41.5889377667});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 16.7208217289241,tol));

  projection.compute({-49.8732305131, -40.7106279818, 48.4922331285});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.82191619179725,tol));

  projection.compute({-0.5889005263, 15.2143434459, -2.7129801701});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 15.4833842886148,tol));

  projection.compute({48.969280533, 1.8857173398, -5.5880641358});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 11.4650237679298,tol));

  projection.compute({-8.4467404794, 45.5121414715, -45.3669887015});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 16.9062391608171,tol));

  projection.compute({30.4369597139, -1.7965056709, 48.9445074922});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 11.656194363566,tol));

  projection.compute({-44.3057219006, 33.0192715316, 47.8292196048});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-5.7458762805, -43.1324416274, 40.1634508698});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.24864425102874,tol));

  projection.compute({-40.9041742286, 3.1722395463, 4.5642140576});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-31.6555538129, -49.6355080975, -48.3113358721});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.01659661932984,tol));

  projection.compute({-7.825023475, 48.957493342, 43.3268881837});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-44.9014713033, 48.0409349306, -44.2031802117});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({22.8517943401, -29.0174949817, 12.8639449658});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 8.9331640386934,tol));
  projection.compute({27.4416171375, 22.6609359834, 15.6104371723});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.0698217030927,tol));

  projection.compute({-30.2095402406, -18.2692825646, 24.9043642426});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 4.04797179207173,tol));

  projection.compute({48.586275195, 41.7056994008, -14.2714379655});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.0251600482992,tol));

  projection.compute({18.275234135, -3.5222361579, -22.7704009846});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 7.30093236511636,tol));

  projection.compute({6.3712748496, -41.5209055373, -17.2412156906});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({1.2251402024, -22.859654237, -49.4462563188});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({-0.2169936664, 45.897933932, 7.9948189473});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-10.3139075266, 7.8029314325, -49.8249060008});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 16.1366226544753,tol));

  projection.compute({42.6518123563, -7.5629428763, -48.0427275868});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 9.59004016090874,tol));

  projection.compute({-45.737014057, -28.2994790833, -30.3337322922});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.68726766634552,tol));

  projection.compute({1.2162083533, -9.9415968917, 14.8779786028});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 5.83796853424094,tol));

  projection.compute({29.9975908268, 19.9978367751, -14.8495243233});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.1388148975356,tol));

  projection.compute({-16.2058498553, -12.1394114393, -24.9289664323});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.92139434484101,tol));

  projection.compute({3.4080482802, -48.8883231296, -43.8845983678});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({46.1560620908, -41.6277617643, 8.0593691012});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 9.38841968826203,tol));

  projection.compute({20.6848680837, 44.7835938049, -28.666853336});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.9169158964716,tol));

  projection.compute({-48.5924598754, 31.7137622655, -23.0120238722});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({42.6690408926, 21.1015188466, 39.1260346347});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 12.5371008351429,tol));

  projection.compute({-9.4118899942, 43.2968541949, -20.261988449});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-49.8420792631, 22.5401981606, -49.1368522864});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-26.8725741547, 28.7433622306, 19.6908688032});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-47.2028514823, -47.8318335478, 15.5816407248});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.61487671307568,tol));

  projection.compute({16.8426310023, 22.1283477601, 43.479231416});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.2391795455986,tol));

  projection.compute({-46.7923978329, -1.5107623076, 43.8335186307});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({46.2800507882, -1.1699410394, 23.3208604033});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 11.4617729152605,tol));

  projection.compute({-25.1075640671, 16.0016334923, -20.8414799398});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 16.8481879784866,tol));

  projection.compute({-41.8020922652, 49.4673997161, 22.9006189261});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-23.2073077342, -44.9300117301, 22.010030305});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.42934882095823,tol));

  projection.compute({37.8241089582, -17.2122999407, -26.5997939168});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 9.23371269834461,tol));

  projection.compute({2.5119125622, 24.8735006316, -33.4853518212});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 15.4768984220725,tol));

  projection.compute({42.3360173555, -22.3200439812, 37.2103834});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 10.3336146388751,tol));

  projection.compute({24.6305152656, 47.4646406236, 24.1349146581});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.4305844795599,tol));

  projection.compute({10.5149867295, -15.3445231101, 39.6555222057});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 9.6202379101399,tol));
  projection.compute({-0.6580345103, 17.6498819923, 21.9638905823});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 14.6181118336781,tol));

  projection.compute({21.9565900378, 4.854384649, -46.3083175459});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 6.94122706590022,tol));

  projection.compute({47.2674666426, 49.1388321385, 13.4482732338});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.0568860076052,tol));

  projection.compute({-25.7504245153, 24.6689192833, -43.3493452116});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 16.8806165747731,tol));

  projection.compute({-30.1640459244, 6.0843163431, 26.2018722371});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({29.3484592714, 46.1486992408, -8.1712008725});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.4953862279618,tol));

  projection.compute({48.3516445841, 45.3574198277, -48.7276976457});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.2416276364779,tol));

  projection.compute({-45.9047522377, -19.3977520193, 2.7042823158});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.94707816295092,tol));

  projection.compute({-48.2935732223, 3.1715559089, -21.2307443243});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({19.803537554, 1.7730305678, 2.7095494572});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 11.947602263324,tol));

  projection.compute({12.4297294125, -49.8548706993, 4.3646752156});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({-23.3992451212, -33.0813171263, -24.6736706582});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.31123002005976,tol));

  projection.compute({-18.7366216764, 11.0967950249, 8.6394815979});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 16.797149017792,tol));

  projection.compute({-26.0479715076, -28.0749771642, 46.442157075});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.94319396750643));

  projection.compute({4.5678507325, -22.0657207407, -2.8295629904});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({47.0529675004, -49.7124435844, 26.3974328415});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 9.43384944631851,tol));

  projection.compute({-26.2698389014, -14.3729289828, -44.3610589459});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.71783305672742,tol));

  projection.compute({20.170082464, 3.6481735081, 24.0622370383});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 12.308375171854,tol));

  projection.compute({-34.3673068957, -16.577460741, -17.3887349513});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.90023537356724,tol));

  projection.compute({48.5242923249, 23.3141597702, -0.400653505});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 12.5578436828865,tol));

  projection.compute({46.056583941, -22.4939376919, -6.313336434});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 9.83541272781669,tol));

  projection.compute({-16.0577103338, 28.7644069077, 44.1796470406});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-49.6267911045, 23.8918445883, 27.4761605437});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-12.9306530873, 47.8111545079, 23.2428442097});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-38.7378187332, 48.3537688844, -21.8518963418});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-45.7704795039, -29.7681998367, 17.3147712682});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.85009884758553,tol));

  projection.compute({-3.9792966815, -33.0479217614, 17.9478132482});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.21367325526119,tol));

  projection.compute({13.1598295938, 48.6314966803, -46.5716411344});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 14.4196487184172,tol));

  projection.compute({-20.2005061182, -4.6676250895, -2.054497065});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 16.4492271092085,tol));

  projection.compute({34.9103078642, -42.3299725132, -10.2239740362});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({-20.2698560759, 36.5470800952, -1.4485126135});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-41.7693495246, 21.7847013249, -0.4240519136});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({41.6287406309, 13.4751787239, -27.3240220627});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 12.328665054228,tol));

  projection.compute({-28.3444108806, -48.405116982, 49.6433370279});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.59899990547848,tol));

  projection.compute({10.5399498532, -40.6306048812, 29.2711783104});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({44.1834146595, -21.1933777068, 13.3290060625});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 10.1362505804196,tol));

  projection.compute({21.063313806, -25.9638172462, -35.295762953});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({-41.1744665313, -49.737137556, -16.550619419});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.38046239308991,tol));

  projection.compute({9.8580917157, 16.7146294223, -20.1967504202});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 14.7825668458449,tol));

  projection.compute({25.3265625217, -13.2317370098, -7.9272767799});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 9.13454235028398,tol));

  projection.compute({34.0210880078, -45.3797400908, -47.5821475487});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.0,tol));

  projection.compute({-44.0322639393, -31.9711322347, -10.1224126109});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.70750485144238,tol));

  projection.compute({-5.9085791725, -21.4804987756, 32.4249613483});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 4.12342850254901,tol));

  projection.compute({7.0652345927, 38.8497738581, 43.4287495881});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 13.6988253950848,tol));

  projection.compute({-29.6910216006, 41.4709048306, 32.7122103342});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 17.0,tol));

  projection.compute({-17.4587137846, -23.6425445292, 4.9992781389});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 3.74607489714158,tol));

  projection.compute({10.8703085039, 39.8229706054, -12.0107919266});
  test.check(Dune::FloatCmp::eq(projection.parameter(), 14.1342363211832,tol));
}

int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);

  testSurfaceGeometry();
  testPointOnSurfaceProjection();
  testPointOnCurveProjection();

  return 0;
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
}

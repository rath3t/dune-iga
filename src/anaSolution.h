//
// Created by Henri on 17.05.2023.
//

#pragma once

struct AnalyticalSolution {
  double E;
  double nu;
  double Tx;
  double R;
  
  double kappa{(3 - nu) / (1 + nu)};
  double G{E / (2 * (1 + nu))};

  double factor{Tx * R / (8 * G)};
  Dune::FieldVector<double, 2> offset;

  AnalyticalSolution(auto pE, auto pnu, auto pTx, auto pR, auto poffset)
      : E{pE},
        nu{pnu},
        Tx{pTx},
        R{pR},
        offset{poffset},
        kappa{(3.0 - pnu) / (1.0 + nu)},
        G{pE / (2.0 * (1 + pnu))},
        factor{pTx * pR / (8.0 * G)}

  {}

  auto toPolar(const auto& pos) -> std::pair<double, double> {
    auto x = pos[0] - offset[0];
    auto y = pos[1] - offset[1];

    auto r     = std::hypot(x, y);
    auto theta = std::atan2(y, x);

    return {r, theta};
  };

  auto stressSolutionThetaR(const auto& pos) -> Dune::FieldVector<double, 3> {
    auto [r, theta] = toPolar(pos);

    auto th    = Tx / 2;
    auto rr2   = std::pow(R, 2) / std::pow(r, 2);
    auto rr4   = std::pow(R, 4) / std::pow(r, 4);
    auto cos2t = std::cos(2 * theta);
    auto sin2t = std::sin(2 * theta);

    Dune::FieldVector<double, 3> sigma_thetar;
    sigma_thetar[0] = th * (1 - rr2) + th * (1 + 3 * rr4 - 4 * rr2) * cos2t;
    sigma_thetar[1] = th * (1 + rr2) - th * (1 + 3 * rr4) * cos2t;
    sigma_thetar[2] = -th * (1 - 3 * rr4 + 2 * rr2) * sin2t;

    return sigma_thetar;
  };

  static auto coordinateTransform(const Dune::FieldVector<double, 3>& stress, double theta)
      -> Dune::FieldVector<double, 3> {
    auto cos2t = std::cos(2 * theta);
    auto sin2t = std::sin(2 * theta);

    const auto s_x  = stress[0];
    const auto s_y  = stress[1];
    const auto s_xy = stress[2];

    auto hpl = 0.5 * (s_x + s_y);
    auto hmi = 0.5 * (s_x - s_y);

    Dune::FieldVector<double, 3> sigma_transformed;
    sigma_transformed[0] = hpl + hmi * cos2t + s_xy * sin2t;
    sigma_transformed[1] = hpl - hmi * cos2t - s_xy * sin2t;
    sigma_transformed[2] = -hmi * sin2t + s_xy * cos2t;

    return sigma_transformed;
  };

  auto stressSolution(const auto& pos) -> Dune::FieldVector<double, 3> {
    auto sigma_thetar = stressSolutionThetaR(pos);

    auto [r, theta] = toPolar(pos);
    theta *= -1;
    return coordinateTransform(sigma_thetar, theta);
  };

  auto displacementSolution(const auto& pos) -> Dune::FieldVector<double, 2> {
    auto [r, theta] = toPolar(pos);

    auto costh  = std::cos(theta);
    auto sinth  = std::sin(theta);
    auto cos3th = std::cos(3 * theta);
    auto sin3th = std::sin(3 * theta);

    auto ra  = r / R;
    auto ar  = R / r;
    auto ar3 = std::pow(R, 3) / std::pow(r, 3);

    Dune::FieldVector<double, 2> res;
    res[0] = factor * (ra * (kappa + 1) * costh + 2 * ar * ((1 + kappa) * costh + cos3th) - 2 * ar3 * cos3th);
    res[1] = factor * (ra * (kappa - 3) * sinth + 2 * ar * ((1 - kappa) * sinth + sin3th) - 2 * ar3 * sin3th);
    return res;
  }

  auto stressLambda() {
    return [&](const auto& pos) { return stressSolution(pos); };
  }

  auto displacementLambda() {
    return [&](const auto& pos) { return displacementSolution(pos); };
  }

  std::array<double, 5> estimateError(const auto& gridView, auto& igaStress, auto& dispGlobalFunc) {
    using namespace Dune::Functions;

    auto localStress        = localFunction(igaStress);
    auto localDisplacements = localFunction(dispGlobalFunc);

    auto localStressAnalytic = localFunction(Dune::Functions::makeAnalyticGridViewFunction(stressLambda(), gridView));
    auto localDisplacementAnalytic
        = localFunction(Dune::Functions::makeAnalyticGridViewFunction(displacementLambda(), gridView));

    std::array<double, 5> l2_error{};

    for (auto& ele : elements(gridView)) {
      localStress.bind(ele);
      localStressAnalytic.bind(ele);

      localDisplacements.bind(ele);
      localDisplacementAnalytic.bind(ele);

      const auto rule = ele.impl().getQuadratureRule();

      for (auto& gp : rule) {
        const auto stress_ana = localStressAnalytic(gp.position());
        const auto stress_fe  = localStress(gp.position());

        const auto displ_ana = localDisplacementAnalytic(gp.position());
        const auto displ_fe  = localDisplacements(gp.position());

        for (int i : std::views::iota(0, 3))
          l2_error[i] += Dune::power(stress_ana[i] - stress_fe[i], 2) * ele.geometry().integrationElement(gp.position())
                         * gp.weight();
        for (int i : std::views::iota(3, 5))
          l2_error[i] += Dune::power(displ_ana[i - 3] - displ_fe[i - 3], 2)
                         * ele.geometry().integrationElement(gp.position()) * gp.weight();
      }
    }

    return l2_error;
  }
};

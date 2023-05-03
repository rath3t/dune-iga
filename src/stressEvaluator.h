//
// Created by Henri on 30.04.2023.
//

#ifndef DUNE_IGA_STRESSEVALUATOR_H
#define DUNE_IGA_STRESSEVALUATOR_H

#include <dune/vtk/vtkwriter.hh>

enum class StressEvaluatorComponents {
  normalStress,
  shearStress,
  vonMieses
};


template <class GridView, typename ElementType, StressEvaluatorComponents comps>
class StressEvaluator2D : public Dune::VTKFunction<GridView> {
 public:
  typedef typename GridView::ctype ctype;
  constexpr static int dim = GridView::dimension;
  typedef typename GridView::template Codim<0>::Entity Entity;

  double evaluate(int comp, const Entity &e, const Dune::FieldVector<ctype, dim> &xi) const override {
    auto index = e.impl().getIndex();
    return evaluateStressComponent(index, xi, comp);
  }

  [[nodiscard]] constexpr int ncomps() const override {
    if constexpr (comps == StressEvaluatorComponents::normalStress)
      return 3;
    if constexpr (comps == StressEvaluatorComponents::shearStress or comps == StressEvaluatorComponents::vonMieses)
      return 1;

  }
  [[nodiscard]] constexpr std::string name() const override {
    if constexpr (comps == StressEvaluatorComponents::normalStress)
      return "normal stress";
    if constexpr (comps == StressEvaluatorComponents::shearStress)
      return "shear stress";
    if constexpr (comps == StressEvaluatorComponents::vonMieses)
      return "von Mieses stress";
  }
  StressEvaluator2D(auto &global_displacement_solution, double lambdaLoad, std::vector<ElementType>* fes)
      : resultRequirements_(Ikarus::ResultRequirements()
                                .insertGlobalSolution(Ikarus::FESolutions::displacement, global_displacement_solution)
                                .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                                .addResultRequest(Ikarus::ResultType::cauchyStress)),
        fes_(fes) {}

 private:
  double evaluateStressComponent(int eleID, auto &xi, int comp) const {
    fes_->at(eleID).calculateAt(resultRequirements_, {xi[0], xi[1]}, res_);
    auto sigma = res_.getResult(Ikarus::ResultType::cauchyStress);

    if constexpr (comps == StressEvaluatorComponents::normalStress)
      return sigma(comp, 0);

    if constexpr (comps == StressEvaluatorComponents::shearStress)
      return sigma(2, 0);

    if constexpr (comps == StressEvaluatorComponents::vonMieses)
      return von_mieses(sigma);
  }

  double von_mieses(const auto& sigma) const requires (comps == StressEvaluatorComponents::vonMieses)  {
    const auto s_x = sigma(0, 0);
    const auto s_y = sigma(1, 0);
    const auto s_xy = sigma(2, 0);

    return std::sqrt(std::pow(s_x, 2) + std::pow(s_y, 2) - s_x * s_y + 3 * std::pow(s_xy, 2));
  }


  Ikarus::ResultRequirements<Eigen::VectorXd, double> resultRequirements_;
  std::vector<ElementType>* fes_;
  mutable Ikarus::ResultTypeMap<double> res_;
};

#endif  // DUNE_IGA_STRESSEVALUATOR_H

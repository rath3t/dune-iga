//
// Created by Henri on 30.04.2023.
//

#ifndef DUNE_IGA_STRESSEVALUATOR_H
#define DUNE_IGA_STRESSEVALUATOR_H

#include <dune/vtk/vtkwriter.hh>
#include <ikarus/finiteElements/feRequirements.hh>

enum class StressEvaluatorComponents {
  normalStress,
  shearStress,
  vonMises,
  principalStress,
  kirchhoff_moments,
  RM_moments,
  RM_forces
};


template <class GridView, typename ElementType, StressEvaluatorComponents comps>
class StressEvaluator2D : public Dune::VTKFunction<GridView> {
 public:
  typedef typename GridView::ctype ctype;
  constexpr static int dim = GridView::dimension;
  typedef typename GridView::template Codim<0>::Entity Entity;

  double evaluate(int comp, const Entity& e, const Dune::FieldVector<ctype, dim>& xi) const override {
    auto index = indexSet.index(e);
    return evaluateStressComponent(index, xi, comp);
  }
  double evaluate(int comp, int index, const Dune::FieldVector<ctype, dim>& xi) {
    assert(index < fes_.size());
    assert(comp < ncomps());
    return evaluateStressComponent(index, xi, comp);
  }

  [[nodiscard]] constexpr int ncomps() const override {
    if constexpr (comps == StressEvaluatorComponents::normalStress or
                  comps == StressEvaluatorComponents::principalStress or
                  comps == StressEvaluatorComponents::RM_forces)
      return 2;
    if constexpr (comps == StressEvaluatorComponents::shearStress or
                  comps == StressEvaluatorComponents::vonMises)
      return 1;
    if constexpr (comps == StressEvaluatorComponents::kirchhoff_moments or
                  comps == StressEvaluatorComponents::RM_moments)
      return 3;
  }
  [[nodiscard]] constexpr std::string name() const override {
    if constexpr (comps == StressEvaluatorComponents::normalStress)
      return "normal stress";
    if constexpr (comps == StressEvaluatorComponents::shearStress)
      return "shear stress";
    if constexpr (comps == StressEvaluatorComponents::vonMises)
      return "von Mises stress";
    if constexpr (comps == StressEvaluatorComponents::principalStress)
      return "principal stress";
    if constexpr (comps == StressEvaluatorComponents::kirchhoff_moments or
                  comps == StressEvaluatorComponents::RM_moments)
      return "moments (x, y, xy)";
    if constexpr (comps == StressEvaluatorComponents::RM_forces)
      return "forces (vx, vy)";
  }
  StressEvaluator2D(GridView& gV, std::vector<ElementType>* fes, auto &global_displacement_solution, double lambdaLoad = 1)
      : indexSet(gV.indexSet()),
        resultRequirements_(Ikarus::ResultRequirements()
                                .insertGlobalSolution(Ikarus::FESolutions::displacement, global_displacement_solution)
                                .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                                .addResultRequest(Ikarus::ResultType::cauchyStress)),
        fes_(fes) {}

  // This object is not copy constructable because of the const ref to the indexSet
  StressEvaluator2D() = delete;
  StressEvaluator2D(const StressEvaluator2D& ) = delete;

 private:
  double evaluateStressComponent(int eleID, auto &xi, int comp) const {
    fes_->at(eleID).calculateAt(resultRequirements_, {xi[0], xi[1]}, res_);
    auto sigma = res_.getResult(Ikarus::ResultType::cauchyStress);

    if constexpr (comps == StressEvaluatorComponents::normalStress or
                  comps == StressEvaluatorComponents::kirchhoff_moments)
      return sigma(comp, 0);

    if constexpr (comps == StressEvaluatorComponents::shearStress)
      return sigma(2, 0);

    if constexpr (comps == StressEvaluatorComponents::vonMises)
      return von_mieses(sigma);

    if constexpr (comps == StressEvaluatorComponents::principalStress)
      return principalStress(sigma)[comp];

    if constexpr (comps == StressEvaluatorComponents::RM_moments) {
      if (comp < 2)
        return sigma(comp, comp);
      else
        return sigma(0, 1);
    }
    if constexpr (comps == StressEvaluatorComponents::RM_forces)
      return sigma(comp, 2);

  }

  double von_mieses(const auto& sigma) const requires (comps == StressEvaluatorComponents::vonMises)  {
    const auto s_x = sigma(0, 0);
    const auto s_y = sigma(1, 0);
    const auto s_xy = sigma(2, 0);

    return std::sqrt(std::pow(s_x, 2) + std::pow(s_y, 2) - s_x * s_y + 3 * std::pow(s_xy, 2));
  }

  std::array<double, 2> principalStress(const auto& sigma) const requires (comps == StressEvaluatorComponents::principalStress) {
    // ref https://www.continuummechanics.org/principalstressesandstrains.html
    const auto s_x = sigma(0, 0);
    const auto s_y = sigma(1, 0);
    const auto s_xy = sigma(2, 0);

    auto t1 = (s_x + s_y) / 2;
    auto t2 = std::sqrt(std::pow((s_x - s_y) / 2, 2) + std::pow(s_xy, 2));

    auto s_1 = t1 + t2;
    auto s_2 = t1 - t2;

    if (s_2 > s_1)
      return {s_2, s_1};
    return {s_1, s_2};
  }

  const GridView::IndexSet& indexSet;
  Ikarus::ResultRequirements<Eigen::VectorXd, double> resultRequirements_;
  std::vector<ElementType>* fes_;
  mutable Ikarus::ResultTypeMap<double> res_;


};

#endif  // DUNE_IGA_STRESSEVALUATOR_H

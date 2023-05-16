//
// Created by Henri on 30.04.2023.
//

#ifndef DUNE_IGA_STRESSEVALUATOR_H
#define DUNE_IGA_STRESSEVALUATOR_H

#include <dune/vtk/vtkwriter.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>

enum class StressEvaluatorComponents {
  normalStress,
  shearStress,
  vonMises,
  principalStress,
  kirchhoff_moments,
  RM_moments,
  RM_forces,
  user_function
};


template <class GridView, typename ElementType, StressEvaluatorComponents comps, Ikarus::ResultType ResType = Ikarus::ResultType::cauchyStress, int user_ncomps = 3>
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
    if constexpr (comps == StressEvaluatorComponents::user_function)
      return user_ncomps;
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
    if constexpr (comps == StressEvaluatorComponents::user_function)
      return user_name;
  }
  StressEvaluator2D(GridView& gV, std::vector<ElementType>* fes, auto &global_displacement_solution, double lambdaLoad = 1)
      : indexSet(gV.indexSet()),
        resultRequirements_(Ikarus::ResultRequirements()
                                .insertGlobalSolution(Ikarus::FESolutions::displacement, global_displacement_solution)
                                .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                                .addResultRequest(ResType)),
        fes_(fes) {}

  StressEvaluator2D(GridView& gV, std::vector<ElementType>* fes, auto &global_displacement_solution, auto userFunction, std::string&& user_name, double lambdaLoad = 1) requires (comps ==  StressEvaluatorComponents::user_function)
      : indexSet(gV.indexSet()),
        resultRequirements_(Ikarus::ResultRequirements()
                                .insertGlobalSolution(Ikarus::FESolutions::displacement, global_displacement_solution)
                                .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                                .addResultRequest(ResType)),
        fes_(fes), user_name(user_name), userFunction_(userFunction) {}

  // This object is not copy constructable because of the const ref to the indexSet
  StressEvaluator2D() = delete;
  StressEvaluator2D(const StressEvaluator2D& ) = delete;

 private:
  double evaluateStressComponent(int eleID, auto& xi, int comp) const {

    if (not (eleID == cachedIndex and Dune::FloatCmp::eq(xi, cachedXi))) {
      fes_->at(eleID).calculateAt(resultRequirements_, Dune::toEigen(xi), res_);
      sigma = res_.getResult(ResType);

      cachedIndex = eleID;
      cachedXi = xi;
    }

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

    if constexpr (comps == StressEvaluatorComponents::user_function)
      return userFunction_(sigma)[comp];
  }

  static double von_mieses(const auto& sigma) requires (comps == StressEvaluatorComponents::vonMises)  {
    const auto s_x = sigma(0, 0);
    const auto s_y = sigma(1, 0);
    const auto s_xy = sigma(2, 0);

    return std::sqrt(std::pow(s_x, 2) + std::pow(s_y, 2) - s_x * s_y + 3 * std::pow(s_xy, 2));
  }

  static std::array<double, 2> principalStress(const auto& sigma) requires (comps == StressEvaluatorComponents::principalStress) {
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
  std::string user_name{};
  std::function<Dune::FieldVector<double, user_ncomps>(Eigen::MatrixXd)> userFunction_{};

  mutable int cachedIndex {-1};
  mutable Dune::FieldVector<double, 2> cachedXi {-1, -1};
  mutable Eigen::MatrixXd sigma;
};

#endif  // DUNE_IGA_STRESSEVALUATOR_H

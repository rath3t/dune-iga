
#pragma once

namespace Dune::IGANEW::DefaultTrim {

  template <typename GridImp>
  struct PatchTrimDataImpl {
    using TrimmingCurve = typename GridImp::GridFamily::TrimmerTraits::TrimmingCurve;

    struct BoundaryLoop {
      void insertTrimCurve(const TrimmingCurve& curve) { curves_.push_back(curve); }

      const auto& curves() const { return curves_; }

     private:
      std::vector<TrimmingCurve> curves_;
    };

    void addLoop() { loops_.push_back({}); }

    void insertTrimCurve(const TrimmingCurve& curve, const int toLoop) {
      assert(loops_.size() > toLoop);
      loops_[toLoop].insertTrimCurve(curve);
    }

    const auto& loops() const { return loops_; }
    // const auto& allCurves() const { return curves_; }

   private:
    std::vector<BoundaryLoop> loops_;
  };

}  // namespace Dune::IGANEW::DefaultTrim

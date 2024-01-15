
#pragma once
#include <clipper2/clipper.core.h>

namespace Dune::IGANEW::DefaultTrim {

  template <typename GridImp>
  struct PatchTrimDataImpl {
    using TrimmingCurve = typename GridImp::GridFamily::TrimmerTraits::TrimmingCurve;
    using PatchData = NURBSPatchData<GridImp::GridFamily::patchDim, GridImp::GridFamily::worldDim, typename GridImp::ctype>;

    struct BoundaryLoop {
      void insertTrimCurve(const TrimmingCurve& curve) { curves_.push_back(curve); }

      const auto& curves() const { return curves_; }

    private:
      std::vector<TrimmingCurve> curves_;
    };



    class CurveManager {
    public:
      using idx_t = u_int64_t;
      explicit CurveManager(const idx_t splitter = 100) : splitter_(splitter) {}

      void addLoop(const BoundaryLoop& loop) {
        Clipper2Lib::PathD path;
        for (idx_t i = loops_.empty() ? splitter_ : loops_.back().back().z + 1;
          const auto& curve : loop.curves()) {
          if (curve.degree().front() == 1) {
            auto p1 = curve.corner(0);
            auto p2 = curve.corner(1);
            path.emplace_back(p1[0], p1[1], i);
            path.emplace_back(p2[0], p2[1], i + splitter_ - 1);
            i += splitter_;
            continue;
          }
          for (const auto v : Utilities::linspace(curve.domain()[0], splitter_)) {
            auto fV = curve.global({v});
            path.emplace_back(fV[0], fV[1], ++i);
          }
        }
        loops_.push_back(path);
        loopIndices_.push_back(loops_.back().back().z);
      }

      [[nodiscard]] auto getIndices(const idx_t val) const -> std::pair<size_t, size_t> {
        size_t loopIdx = 0;
        if (loops_.size() > 1) {
          auto it = std::ranges::find_if(loopIndices_, [&](const idx_t t) {
                return val > t;
          });
          if (it == loopIndices_.end())
            DUNE_THROW(Dune::IOError, "Invalid z-Value");
          loopIdx = std::ranges::distance(loopIndices_.begin(), it);
        }
        // \todo multiplt loop e.g. - size of all loops before
        size_t curveIdx = static_cast<size_t>(std::floor(val / splitter_) - 1);

        return std::make_pair(loopIdx, curveIdx);
      }

      [[nodiscard]] auto getSplitter() const -> idx_t {
        return splitter_;
      }

      Clipper2Lib::PathsD loops_;
    private:
      idx_t splitter_{};
      std::vector<idx_t> loopIndices_{};
    };

    void addLoop() { loops_.push_back({}); }

    void insertTrimCurve(const TrimmingCurve& curve, const int toLoop) {
      assert(loops_.size() > toLoop);
      loops_[toLoop].insertTrimCurve(curve);
    }

    const auto& loops() const { return loops_; }
    const auto& clipperLoops() const {
      if (not finished_)
        DUNE_THROW(Dune::GridError, "Call prepare() before quering for loops");
      return manager_.loops_;
    }

    auto getIndices(const typename CurveManager::idx_t val) const -> std::pair<size_t, size_t> {
      return manager_.getIndices(val);
    }
    auto getCurve(const typename CurveManager::idx_t val) const -> const TrimmingCurve& {
      auto [loopIdx, curveIdx] = manager_.getIndices(val);
      return loops_[loopIdx].curves()[curveIdx];
    }
    auto getCurve(const std::pair<size_t, size_t>& indices) const -> const TrimmingCurve& {
      return loops_[indices.first].curves()[indices.second];
    }
    auto getSplitter() const -> typename CurveManager::idx_t {
      return manager_.getSplitter();
    }

    void prepare(PatchData* patchData) {
      std::ranges::for_each(loops_, [&](const auto& loop) { manager_.addLoop(loop);});


      finished_ = true;
    }

  private:
    bool finished_ = false;
    std::vector<BoundaryLoop> loops_;
    CurveManager manager_;

  };

}  // namespace Dune::IGANEW::DefaultTrim

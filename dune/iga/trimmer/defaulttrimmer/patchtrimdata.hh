
#pragma once
#include <clipper2/clipper.core.h>

namespace Dune::IGANEW::DefaultTrim {

  template <typename GridImp>
  struct PatchTrimDataImpl {
    using TrimmingCurve = typename GridImp::GridFamily::TrimmerTraits::TrimmingCurve;
    using ParameterType = typename GridImp::GridFamily::Trimmer::ParameterType;
    using ctype = typename GridImp::GridFamily::ctype;

    struct BoundaryLoop {
      void insertTrimCurve(const TrimmingCurve& curve) { curves_.push_back(curve); }

      const auto& curves() const { return curves_; }

    private:
      std::vector<TrimmingCurve> curves_;
    };

    struct PointInPatch {
      FieldVector<ctype, 2> pt;
      size_t curveIdxI;
      size_t curveIdxJ;
    };

    class CurveManager {
      friend PatchTrimDataImpl;

    public:
      using idx_t = u_int64_t;

      void addLoop(const BoundaryLoop& loop) {
        Clipper2Lib::PathD path;
        for (idx_t i = loops_.empty() ? splitter_ : loops_.back().back().z + 1; const auto& curve : loop.curves()) {
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
            path.emplace_back(fV[0], fV[1], i++);
          }
        }
        loops_.push_back(path);
        loopIndices_.emplace_back(loops_.back().back().z, loop.curves().size());
      }

      [[nodiscard]] auto getIndices(const idx_t val) const -> std::pair<size_t, size_t> {
        size_t loopIdx = 0;
        if (loops_.size() > 1) {
          auto it = std::ranges::find_if(loopIndices_, [&](const std::pair<idx_t, size_t>& t) { return val > t.first; });
          if (it == loopIndices_.end()) DUNE_THROW(Dune::IOError, "Invalid z-Value");
          loopIdx = std::ranges::distance(loopIndices_.begin(), it) + 1;
        }

        auto curveIdx = static_cast<size_t>(std::floor(val / splitter_) - 1);
        for (const auto& [_, sizeOfLoop] : std::ranges::take_view(loopIndices_, static_cast<long>(loopIdx))) {
          curveIdx -= sizeOfLoop;
        }

        return std::make_pair(loopIdx, curveIdx);
      }

    private:
      Clipper2Lib::PathsD loops_;
      idx_t splitter_{};
      std::vector<std::pair<idx_t, size_t>> loopIndices_{};
    };

    void addLoop() { loops_.push_back({}); }

    void insertTrimCurve(const TrimmingCurve& curve, const int toLoop) {
      assert(loops_.size() > toLoop);
      loops_[toLoop].insertTrimCurve(curve);
    }

    const auto& loops() const { return loops_; }
    const auto& clipperLoops() const {
      if (not finished_) DUNE_THROW(Dune::GridError, "Call trimmer.setup() before quering for loops");
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
    auto getSplitter() const -> typename CurveManager::idx_t { return manager_.splitter_; }

    auto getPointsInPatch(size_t loopIndex) const -> const std::vector<PointInPatch>& {
      return pointsInPatch_.at(loopIndex);
    }

    void prepare(const ParameterType& parameters, const std::array<std::vector<ctype>, 2>& coordinates) {
      manager_.splitter_ = parameters.splitter;
      std::ranges::for_each(loops_, [&](const auto& loop) { manager_.addLoop(loop); });

      std::array domainU{coordinates[0].front(), coordinates[0].back()};
      std::array domainV{coordinates[1].front(), coordinates[1].back()};

      // Determine points of the outer boundary that are not on the edges
      auto isInsidePatch = [&](const FieldVector<ctype, 2>& pt) {
        const bool isInsideU = FloatCmp::gt(pt[0], domainU.front()) and FloatCmp::lt(pt[0], domainU.back());
        const bool isInsideV = FloatCmp::gt(pt[1], domainV.front()) and FloatCmp::lt(pt[1], domainV.back());

        return isInsideU and isInsideV;
      };

      pointsInPatch_.push_back({});
      for (size_t i = 0; const auto& curve : loops_.front().curves()) {
        if (const auto pt = curve.corner(1); isInsidePatch(pt)) {
          pointsInPatch_.front().emplace_back(pt, i, (i + 1) % loops_.front().curves().size());
        }
        ++i;
      }

      finished_ = true;
    }

  private:
    std::vector<std::vector<PointInPatch>> pointsInPatch_;
    bool finished_ = false;
    std::vector<BoundaryLoop> loops_;
    CurveManager manager_;
  };

}  // namespace Dune::IGANEW::DefaultTrim

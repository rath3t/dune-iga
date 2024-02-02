
#pragma once
#include <clipper2/clipper.core.h>

namespace Dune::IGANEW::DefaultTrim {

  namespace Impl {
    template<typename TrimmingCurve>
    struct BoundaryLoop {
      void insertTrimCurve(const TrimmingCurve& curve) { curves_.push_back(curve); }

      const auto& curves() const { return curves_; }

      [[nodiscard]] size_t size() const { return curves_.size(); }

    private:
      std::vector<TrimmingCurve> curves_;
    };

    template <typename ctype>
    struct PointInPatch {
      FieldVector<ctype, 2> pt;
      size_t curveIdxI;
      size_t curveIdxJ;
    };



    class CurveLoopIndexEncoder {
      static constexpr int loopBits = 4;
      static constexpr u_int64_t loopMask{(1ULL << loopBits) - 1};
    public:

      static int64_t encode(int curveIndex, int loopIndex)  {
        // the left hand side shifts the curve index value away from loop values and then inserts the curve index
        // example  curveBits = 4 and loopBits = 3;
        // curveIndex = 3 and loopIndex = 2 = 00010_b
        // curveIndex = 3 = 00011_b
        // curveIndex<< loopBits  = 11000_b
        // curveIndex<< loopBits  = 11000_b
        // (curveIndex<< loopBits) | loopIndex  = 11000_b | 00010_b = 11010_b = 26
        return (static_cast<u_int64_t>(curveIndex) << loopBits) | loopIndex;
      }

      static std::pair<int, int> decode(int64_t singleIndex)  {
        // example  curveBits = 4 and loopBits = 3;
        // loopMask = 00111
        // example  singleIndex = 11010_b = 26
        // curveIndex = 3 and loopIndex = 2 = 00010_b
        // curveIndex = 3 = 00011_b
        // loopMask =
        // (singleIndex >> loopBits)  = 00011_b = curveIndex
        int curveIndex = (singleIndex >> loopBits);
        // since the loop index is stored in the right-most bits, we only have to apply the bitwise and, which returns simply the loop index
        // 11010_b AND
        // 00111_b
        // 00010_b = 2
        int loopIndex = singleIndex & loopMask;
        return {loopIndex,curveIndex};
      }

    };
  }

  template <typename GridImp>
  struct PatchTrimDataImpl {
    using TrimmingCurve = typename GridImp::GridFamily::TrimmerTraits::TrimmingCurve;
    using ParameterType = typename GridImp::GridFamily::Trimmer::ParameterType;
    using ctype = typename GridImp::GridFamily::ctype;



    class CurveManager {
      friend PatchTrimDataImpl;

    public:
      using idx_t = u_int64_t;
      static constexpr idx_t indexOffSet = 4; //Since the first 4 values are reserved for the corners of the untrimmed square all indices start at 4

      void addLoop(const Impl::BoundaryLoop<TrimmingCurve>& loop) {
        Clipper2Lib::PathD path;
        // since the first
        for (int curveIndex=0; const auto& curve : loop.curves()) {
          auto index= getZValue(curveIndex,loops_.size());
          if (curve.affine() == 1 and loops_.empty()) {
            auto p1 = curve.corner(0);
            auto p2 = curve.corner(1);
            path.emplace_back(p1[0], p1[1], index);
            path.emplace_back(p2[0], p2[1], index );
          }else {
            for (const auto v : Utilities::linspace(curve.domain()[0], splitter_)) {
              auto fV = curve.global({v});
              path.emplace_back(fV[0], fV[1], index);
            }
          }
          ++curveIndex;
        }
        loops_.push_back(path);
      }

      [[nodiscard]] auto getIndices(const idx_t val) const -> std::pair<size_t, size_t> {
        return Impl::CurveLoopIndexEncoder::decode(val-indexOffSet);
      }

      [[nodiscard]] auto getZValue(int edgeIndex, int loopIndex )const {
        return Impl::CurveLoopIndexEncoder::encode(edgeIndex, loopIndex)+indexOffSet;
      }

    private:
      Clipper2Lib::PathsD loops_;
      idx_t splitter_{};

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

    auto getZValue(int edgeIndex, int loopIndex ) const  {
      return manager_.getZValue(edgeIndex,loopIndex);
    }

    auto getCurve(const typename CurveManager::idx_t val) const -> const TrimmingCurve& {
      auto [loopIdx, curveIdx] = manager_.getIndices(val);
      return loops_[loopIdx].curves()[curveIdx];
    }
    auto getCurve(const std::pair<size_t, size_t>& indices) const -> const TrimmingCurve& {
      return loops_[indices.first].curves()[indices.second];
    }
    auto getSplitter() const -> typename CurveManager::idx_t { return manager_.splitter_; }

    auto getPointsInPatch(size_t loopIndex) const -> const std::vector<Impl::PointInPatch<ctype>>& {
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

      for (const auto& loop : loops()) {
        pointsInPatch_.push_back({});
        for (size_t i = 0; const auto& curve : loop.curves()) {
          if (const auto pt = curve.corner(1); isInsidePatch(pt)) {
            pointsInPatch_.back().emplace_back(pt, i, (i + 1) % loop.size());
          }
          ++i;
        }
      }

      finished_ = true;
    }

  private:
    std::vector<std::vector<Impl::PointInPatch<ctype>>> pointsInPatch_;
    bool finished_ = false;
    std::vector<Impl::BoundaryLoop<TrimmingCurve>> loops_;
    CurveManager manager_;
  };

}  // namespace Dune::IGANEW::DefaultTrim

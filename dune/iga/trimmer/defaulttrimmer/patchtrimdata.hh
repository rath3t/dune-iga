
#pragma once
#include <clipper2/clipper.core.h>

namespace Dune::IGANEW::DefaultTrim {

namespace Impl {
  template <typename TrimmingCurve>
  struct BoundaryLoop
  {
    void insertTrimCurve(const TrimmingCurve& curve) {
      curves_.push_back(curve);
    }

    const auto& curves() const {
      return curves_;
    }

    [[nodiscard]] size_t size() const {
      return curves_.size();
    }

  private:
    std::vector<TrimmingCurve> curves_;
  };

  template <typename ctype>
  struct PointInPatch
  {
    FieldVector<ctype, 2> pt;
    size_t curveIdxI;
    size_t curveIdxJ;
  };

  class CurveLoopIndexEncoder
  {
    static constexpr int loopBits   = 4;
    static constexpr int curveBits  = 8;
    static constexpr int sampleBits = 12;

    static constexpr u_int64_t loopMask{(1ULL << loopBits) - 1};
    static constexpr u_int64_t curveMask{(1ULL << curveBits) - 1};
    static constexpr u_int64_t sampleMask{(1ULL << sampleBits) - 1};

  public:
    static int64_t encode(int loopIndex, int curveIndex, int sampleIndex) {
      // the left hand side shifts the curve index value away from loop values and then inserts the curve index
      // example  curveBits = 4 and loopBits = 3 and sampleBits=12;
      // curveIndex = 3 and loopIndex = 2 = 00010_b , sampleIndex= 200= 11001000_b
      // curveIndex = 3 = 00011_b
      // curveIndex<< (loopBits+ sampleBits)  = 110,000,000,000,000,000_b // shift the curve index 16 bits to the left
      // loopIndex<< sampleBits  =                   10,000,000,000,000_b // shift the loop index 12 bits to the left
      // sampleIndex =                                       11,001,000_b // shift the loop index 12 bits to the left
      // bitwise or of the three values:        110,010,000,011,001,000_b //
      // (curveIndex<< loopBits) | loopIndex  = 11000_b | 00010_b = 11010_b = 205000_d
      return (static_cast<u_int64_t>(curveIndex) << (loopBits + sampleBits)) |
             (static_cast<u_int64_t>(loopIndex) << sampleBits) | sampleIndex;
    }

    struct IndexResult
    {
      int loop;
      int curve;
      int sample;

      friend std::ostream& operator<<(std::ostream& os, const IndexResult& res) {
        os << "Loop: " << res.loop << " Curve: " << res.curve << " Sample: " << res.sample;
        return os;
      }
    };

    static IndexResult decode(int64_t singleIndex) {
      // since the sample index is stored in the right-most bits, we only have to apply the bitwise "and" with
      // sampleMask, which returns simply the sample index
      int sampleIndex = singleIndex & sampleMask;
      // shift single index to have the loop index at the rightmost position and apply the bitwise and with the loop
      // mask
      const int loopAndCurveIndex = singleIndex >> sampleBits;
      int loopIndex               = loopAndCurveIndex & loopMask;
      // shift again to have the curve index at the rightmost position and apply the bitwise and with the curve mask
      int curveIndex = (loopAndCurveIndex >> loopBits) & curveMask;
      return {loopIndex, curveIndex, sampleIndex};
    }
  };
} // namespace Impl

template <typename GridImp>
struct PatchTrimDataImpl
{
  using TrimmingCurve = typename GridImp::GridFamily::TrimmerTraits::TrimmingCurve;
  using ParameterType = typename GridImp::GridFamily::Trimmer::ParameterType;
  using ctype         = typename GridImp::GridFamily::ctype;

  class CurveManager
  {
    friend PatchTrimDataImpl;

  public:
    using idx_t                        = u_int64_t;
    static constexpr idx_t indexOffSet = 5; // Since the 0 is reserved for untouched points and the following 4
    // values are reserved for the corners of the untrimmed square all indices start at 5
    void addLoop(const Impl::BoundaryLoop<TrimmingCurve>& loop) {
      Clipper2Lib::PathD path;
      // since the first
      for (int curveIndex = 0; const auto& curve : loop.curves()) {
        if (curve.affine() == 1 and loops_.empty()) {
          auto p1 = curve.corner(0);
          auto p2 = curve.corner(1);
          path.emplace_back(p1[0], p1[1], getZValue(loops_.size(), curveIndex, 0));
          path.emplace_back(p2[0], p2[1], getZValue(loops_.size(), curveIndex, 1));
        } else {
          for (int vi = 0; const auto v : Utilities::linspace(curve.domain()[0], splitter_)) {
            auto fV = curve.global({v});
            path.emplace_back(fV[0], fV[1], getZValue(loops_.size(), curveIndex, vi++));
          }
        }
        ++curveIndex;
      }
      loops_.push_back(path);
    }

    [[nodiscard]] auto getIndices(const idx_t val) const -> Impl::CurveLoopIndexEncoder::IndexResult {
      if (val == 0) // untouched point usually created by Clipperlib as Intersection point
        return Impl::CurveLoopIndexEncoder::IndexResult{-1, -1, -1};
      assert(val >= indexOffSet &&
             "Passed index not decodable, since it is smaller than the offset and non-zero, which usually means that "
             "the point is a host vertex");
      return Impl::CurveLoopIndexEncoder::decode(val - indexOffSet);
    }

    [[nodiscard]] auto getZValue(int loopIndex, int curveIndex, int sampleIndex) const {
      return Impl::CurveLoopIndexEncoder::encode(loopIndex, curveIndex, sampleIndex) + indexOffSet;
    }

  private:
    Clipper2Lib::PathsD loops_;
    idx_t splitter_{};
  };

  void addLoop() {
    loops_.push_back({});
  }

  void insertTrimCurve(const TrimmingCurve& curve, const int toLoop) {
    assert(loops_.size() > toLoop);
    loops_[toLoop].insertTrimCurve(curve);
  }

  const auto& loops() const {
    return loops_;
  }
  const auto& clipperLoops() const {
    if (not finished_)
      DUNE_THROW(Dune::GridError, "Call trimmer.setup() before quering for loops");
    return manager_.loops_;
  }

  auto getIndices(const typename CurveManager::idx_t val) const -> Impl::CurveLoopIndexEncoder::IndexResult {
    return manager_.getIndices(val);
  }

  auto getZValue(int loopIndex, int edgeIndex, int sampleIndex) const {
    return manager_.getZValue(loopIndex, edgeIndex, sampleIndex);
  }

  auto getCurve(const typename CurveManager::idx_t val) const -> const TrimmingCurve& {
    return getCurve(manager_.getIndices(val));
  }
  auto getCurve(const Impl::CurveLoopIndexEncoder::IndexResult& indices) const -> const TrimmingCurve& {
    return loops_[indices.loop].curves()[indices.curve];
  }

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

} // namespace Dune::IGANEW::DefaultTrim

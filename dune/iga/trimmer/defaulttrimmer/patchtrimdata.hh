
#pragma once

namespace Dune {
  namespace IGANEW {
    namespace DefaultTrim {

      template <typename TrimmerType_>
      struct PatchTrimData {
        using TrimmingCurve = typename TrimmerType_::TrimmingCurve;

        void insertTrimCurve(const TrimmingCurve& curve) { curves.push_back(curve); }

        std::vector<TrimmingCurve> curves;

        //???
      };

    }  // namespace DefaultTrim
  }    // namespace IGANEW
}  // namespace Dune


#pragma once
#include "../io/ibrareader.hh"

#include <nlohmann/json.hpp>

#include <dune/grid/common/gridfactory.hh>

namespace Dune {

template <int dim_, int dimworld_, template <int, int, typename> typename TrimmerType_, typename ScalarType>
class GridFactory<IGANEW::PatchGrid<dim_, dimworld_, TrimmerType_, ScalarType>>
{
  /** @brief The grid world dimension */
  constexpr static int dimworld = dimworld_;
  constexpr static int dim      = dim_;
  using PatchGrid               = IGANEW::PatchGrid<dim, dimworld, TrimmerType_, ScalarType>;
  using TrimmerType             = typename PatchGrid::Trimmer;
  using PatchTrimData           = typename TrimmerType::PatchTrimData;

public:
  using TrimParameterType = typename TrimmerType::ParameterType;

  /** @brief Type used by the grid for coordinates */
  typedef typename PatchGrid::ctype ctype;

  /** @brief Insert a patch into the grid
      @param patchData The patch data
      @param patchTrimData Trimming data for this patch
   */
  void insertPatch(const IGANEW::NURBSPatchData<dim, dimworld, ctype>& patchData,
                   const std::optional<PatchTrimData>& patchTrimData = std::nullopt) {
    patchData_     = patchData;
    patchTrimData_ = patchTrimData;
  }

  void insertTrimmingCurve(const IGANEW::NURBSPatchData<dim - 1, dim, ctype>& curve) {
    DUNE_THROW(NotImplemented, "insertTrimmingCurve not yet implemented, stay tuned ");
  }

  void insertTrimParameters(const TrimParameterType& par) {
    parameters_ = par;
  }

  void insertJson(const std::string& filename, const bool trim = true, std::array<int, 2> preKnotRefine = {0, 0}) {
    json_                      = filename;
    auto [patchData, trimData] = IGANEW::IbraReader<dim, dimworld, PatchGrid>::read(filename, trim, preKnotRefine);
    insertPatch(patchData, trimData);
  }

  /** @brief Finalize grid creation and hand over the grid

     The receiver takes responsibility of the memory allocated for the grid
   */
  std::unique_ptr<PatchGrid> createGrid() const {
    if (patchTrimData_.has_value()) {
      auto grid = std::make_unique<PatchGrid>(patchData_, patchTrimData_, parameters_);
      return grid;
    }
    auto grid = std::make_unique<PatchGrid>(patchData_);
    return grid;
  }

  IGANEW::NURBSPatchData<dim, dimworld, ctype> patchData_;
  std::optional<PatchTrimData> patchTrimData_;
  std::string json_;
  TrimParameterType parameters_{};
};

} // namespace Dune


#pragma once

#include <nlohmann/json.hpp>
namespace Dune {

  template <int dim_, int dimworld_, template<int,typename > typename TrimmerType_, typename ScalarType>
  class GridFactory<IGANEW::PatchGrid<dim_,dimworld_,TrimmerType_,ScalarType>>
  {

    /** @brief The grid world dimension */
    constexpr static int dimworld     = dimworld_;
    constexpr static int dim          = dim_;
    using PatchGrid=IGANEW::PatchGrid<dim,dimworld,TrimmerType_,ScalarType>;
    using TrimmerType                 = typename PatchGrid::TrimmerType;
    using PatchTrimData               = typename TrimmerType::PatchTrimData;
    using TrimParameterType           = typename TrimmerType::ParameterType;

    /** @brief Type used by the grid for coordinates */
    typedef typename PatchGrid::ctype ctype;

   public:
    /** @brief Forward setting to trimmer
    @param parameter The parameters
 */
    void setupTrimmer(const TrimParameterType& parameter) { trimmer.setup(parameter); }

    /** @brief Insert a patch into the grid
        @param patchData The patch data
        @param patchTrimData Trimming data for this patch
     */
    void insertPatch(const IGANEW::NURBSPatchData<dim, dimworld, ctype>& patchData,
                     const std::optional<PatchTrimData>& patchTrimData = std::nullopt) {
      patchData_     = patchData;
      patchTrimData_ = patchTrimData;
    }

        /** @brief Insert a trimming curve into the grid
        @param patchData The patch data
        @param patchTrimData Trimming data for this patch
     */
    void insertTrimmingCurve(const IGANEW::NURBSPatchData<dim-1,dim,ctype>& curve) {
      trimCurves.push_back(curve);
    }

    /** @brief Insert a patch into the grid
    @param patchData The patch data
    @param patchTrimData Trimming data for this patch
 */
    void insertJson(const std::string& filename) { json_ = filename; }

    /** @brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    std::unique_ptr<PatchGrid> createGrid() {
      if (patchTrimData_) {
        auto grid = std::make_unique<PatchGrid>(patchData_, patchTrimData_);

        return grid;
      } else
        return std::make_unique<PatchGrid>(patchData_);
      // create Grid and setup Element trimming inf through additional (private) constructor
    }

    IGANEW::NURBSPatchData<dim, dimworld, ctype> patchData_;
    std::optional<PatchTrimData> patchTrimData_;
    std::vector<IGANEW::NURBSPatchData<dim-1,dim,ctype>> trimCurves;
    std::string json_;
    TrimmerType trimmer;
  };

}  // namespace Dune::IGANEW

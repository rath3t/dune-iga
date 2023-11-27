
#pragma once

#include <nlohmann/json.hpp>
namespace Dune::IGANEW {

  template <class GridType>
  class GridFactory {
    typedef GridFactoryInterface<GridType> Base;

    /** \brief The grid world dimension */
    constexpr static int dimworld = GridType::dimensionworld;
    constexpr static int dim      = GridType::dimension;
    using TrimmerType             = typename GridType::TrimmerType;
    using PatchTrimData           = typename TrimmerType::PatchTrimData;
    using TrimParameterType       = typename TrimmerType::ParameterType;
    using UntrimmedParameterSpaceGrid       = typename TrimmerType::UntrimmedParameterSpaceGrid;
    using ParameterSpaceGrid       = typename TrimmerType::ParameterSpaceGrid;

    /** \brief Type used by the grid for coordinates */
    typedef typename GridType::ctype ctype;

   public:
    void setupTrimmer(const TrimParameterType& parameter) { trimmer.parameter = parameter; }

    /** \brief Insert a patch into the grid
        \param patchData The patch data
        \param patchTrimData Trimming data for this patch
     */
    void insertPatch(const NURBSPatchData<dim, dimworld, ctype>& patchData,
                     const std::optional<PatchTrimData>& patchTrimData = std::nullopt) {
      patchData_     = patchData;
      patchTrimData_ = patchTrimData;
    }

    /** \brief Insert a patch into the grid
    \param patchData The patch data
    \param patchTrimData Trimming data for this patch
 */
    void insertJson(const std::string&  filename) {
      json_     = filename;
    }

    /** \brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    std::unique_ptr<GridType> createGrid() {
      if (patchTrimData_) {
        trimmer.trimElements(patchData_, patchTrimData_.value());
        auto uniqueKnotSpans= Splines::createUniqueKnotSpans(patchData_);
        //create SubGrid...
        auto hostGrid =std::make_unique<UntrimmedParameterSpaceGrid>(
              uniqueKnotSpans);

        ParameterSpaceGrid subGrid(hostGrid);
        subGrid.createBegin();
        for (auto hostEntity : elements(hostGrid.leafGridView())) {
          //if decide which elements are full or trim and add them to the subgrid
          //subGrid.insert(hostEntity);
        }
        subGrid.createEnd();

        return std::make_unique<GridType>(patchData_, patchTrimData_, std::move(trimmer));

      }
      // create Grid and setup Element trimming inf through additional (private) constructor
    }

    NURBSPatchData<dim, dimworld, ctype> patchData_;
    std::optional<PatchTrimData> patchTrimData_;
    std::string json_;
    TrimmerType trimmer;
  };

}  // namespace Dune::IGANEW

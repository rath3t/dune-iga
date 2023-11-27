
#pragma once

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

    /** \brief Type used by the grid for coordinates */
    typedef typename GridType::ctype ctype;

   public:
    void setupTrimmer(const TrimParameterType& parameter) { trimmer.parameter = parameter; }

    /** \brief Insert a patch into the grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering

        Make sure the inserted element is not inverted (this holds even
        for simplices).  There are grids that can't handle inverted tets.
     */
    void insertPatch([[maybe_unused]] const NURBSPatchData<dim, dimworld, ctype>& patchData,
                     const std::optional<PatchTrimData>& patchTrimData = std::nullopt) {
      patchData_     = patchData;
      patchTrimData_ = patchTrimData;
    }

    /** \brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    std::unique_ptr<GridType> createGrid() {
      auto elementTrimInfo = trimmer.trimElements(patchData_, patchTrimData_);
      // create Grid and setup Element trimming inf through additional (private) constructor
      return std::make_unique<GridType>(patchData_, patchTrimData_, trimmer);
    }

    NURBSPatchData<dim, dimworld, ctype> patchData_;
    std::optional<PatchTrimData> patchTrimData_;
    TrimmerType trimmer;
  };

}  // namespace Dune::IGANEW


#pragma once

namespace Dune {
  namespace IGANEW {
    namespace DefaultTrim {
// : TrimPatchEntitiy
  template <int codim_, int dim, class GridImp>
      class TrimmedParameterSpaceGridEntity {
        using ctype = typename GridImp::ctype  ;

        static constexpr int mydimension = dim;
        // [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }
        using LocalCoordinate = FieldVector<ctype, mydimension>;

    using TrimmerType= typename GridImp::TrimmerType;
    using GlobalIdSetIdType= typename TrimmerType::GlobalIdSetIdType;
    using ElementTrimData= typename TrimmerType::ElementTrimData;
    using HostParameterSpaceGridEntity= typename TrimmerType::template HostParameterSpaceGridEntity<codim_>;
    using LocalParameterSpaceGeometry= typename TrimmerType::template LocalParameterSpaceGeometry<codim_>;
  public:
    TrimmedParameterSpaceGridEntity(const GridImp* grid,const HostParameterSpaceGridEntity& untrimmedElement, GlobalIdSetIdType id,
      const std::optional<std::reference_wrapper<const ElementTrimData>>& trimData=std::nullopt) :grid_{grid},hostEntity_{untrimmedElement},id_{id},trimData_{trimData},localId_{0}{
    }

    auto & id() const  {
      return id_;
    }


  private:
    HostParameterSpaceGridEntity hostEntity_;
    GlobalIdSetIdType id_;
    std::optional<std::reference_wrapper<const ElementTrimData>> trimData_;
    size_t localId_;
  public:
    [[nodiscard]] bool equals(const TrimmedParameterSpaceGridEntity& other) const { return hostEntity_== other.hostEntity_; }

    //! returns true if father entity exists
    [[nodiscard]] bool hasFather() const {
      //@todo Trim this is crasy
      return hostEntity_.hasFather();
    }

    //! Create EntitySeed
    [[nodiscard]] auto seed() const {
      return hostEntity_.seed();
    }

    //! Level of this element
    [[nodiscard]] int level() const {
      return hostEntity_.level();
    }

    /** @brief The partition type for parallel computing */
    [[nodiscard]] PartitionType partitionType() const {
      //@todo Trim this is crasy

      return hostEntity_.partitionType();
    }

    //! Geometry of this entity
    [[nodiscard]] decltype(auto) geometry() const {
      //@todo Trim this is crasy
      // if(trimData_)
      //   return trimData_.template geometry<codim_>(localId_);
      return hostEntity_.geometry();
    }

    /** @brief Return the number of subEntities of codimension codim.
     */
    [[nodiscard]] unsigned int subEntities(unsigned int codim) const {
      //@todo Trim this is crasy
      // if(trimData_)
      //   return trimData_. subEntities(codim,localId_);
      return hostEntity_.subEntities(codim);
    }

    /** @brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template <int cc> requires (codim_==0)
    [[nodiscard]] decltype(auto) subEntity(int i) const {
      // if(trimData_)
      //   return trimData_.template subEntity<codim_,cc>(i,localId_);
      auto id = grid_->entityContainer().subId(id_,i,cc);
      return TrimmedParameterSpaceGridEntity<cc,mydimension,GridImp>(grid_,hostEntity_.template subEntity<cc>(i),id);
    }

    //! First level intersection
    template <typename =void> requires (codim_==0)
    [[nodiscard]] decltype(auto) ilevelbegin()  const {
      // if(trimData_)
      //   return trimData_.template ilevelbegin<codim_>(localId_);
      return hostEntity_.ilevelbegin();
    }

    //! Reference to one past the last neighbor
    template <typename=void > requires (codim_==0)
    decltype(auto) ilevelend() const {
      return hostEntity_.ilevelend();
    }

    //! First leaf intersection
    template <typename=void > requires (codim_==0)
    decltype(auto) ileafbegin() const {
      return hostEntity_.ileafbegin();
    }

    //! Reference to one past the last leaf intersection
    template <typename=void > requires (codim_==0)
    decltype(auto) ileafend() const {
      return hostEntity_.ileafend();
    }

    //! returns true if Entity has NO children
    template <typename=void > requires (codim_==0)
    bool isLeaf() const {
      return hostEntity_.isLeaf();
    }

    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    template <typename=void > requires (codim_==0)
    decltype(auto) father() const {
      return hostEntity_.father();
    }

    /** @brief Location of this element relative to the reference element element of the father.
     * This is sufficient to interpolate all dofs in conforming case.
     * Nonconforming may require access to neighbors of father and
     * computations with local coordinates.
     * On the fly case is somewhat inefficient since dofs  are visited several times.
     * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
     * implementation of numerical algorithms is only done for simple discretizations.
     * Assumes that meshes are nested.
     */
    template <typename=void > requires (codim_==0)
    decltype(auto) geometryInFather() const {
      return hostEntity_.geometryInFather();
    }

    /** @brief Inter-level access to son elements on higher levels<=maxlevel.
     * This is provided for sparsely stored nested unstructured meshes.
     * Returns iterator to first son.
     */
    template <typename=void > requires (codim_==0)
    decltype(auto) hbegin(int maxLevel) const {
      return hostEntity_.hbegin(maxLevel);
    }

    //! Returns iterator to one past the last son
    template <typename=void > requires (codim_==0)
    decltype(auto) hend(int maxLevel) const {
      return hostEntity_.hend(maxLevel);
    }

    //! @todo Please doc me !
    template <typename=void > requires (codim_==0)
    bool wasRefined() const {
      return hostEntity_.wasRefined();
    }

    //! @todo Please doc me !
    template <typename=void > requires (codim_==0)

    bool mightBeCoarsened() const {
      return hostEntity_.mightBeCoarsened();
    }

    const auto& hostEntity()const  {
      return hostEntity_;
    }

    const GridImp* grid_;


      };

    }  // namespace DefaultTrim
  }    // namespace IGANEW
}  // namespace Dune

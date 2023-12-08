
#pragma once

#include <variant>

namespace Dune {
  namespace IGANEW {
    namespace DefaultTrim {

      template <class GridImp>
      class PatchGridGlobalIdSet
          : public IdSet<GridImp, PatchGridGlobalIdSet<GridImp>,
                         typename std::remove_const<GridImp>::type::Traits::GlobalIdSet::IdType> {
        typedef typename std::remove_const<GridImp>::type::ParameterSpaceGrid ParameterSpaceGrid;

        using Trimmer = typename GridImp::Trimmer;
        // friend class GridImp::TrimmerType;

        // using TrimmingCurve= typename Trimmer::TrimmingCurve;
        using UntrimmedParameterSpaceGrid = typename Trimmer::UntrimmedParameterSpaceGrid;
        // using HostIdType= typename Trimmer::HostIdType;

       public:
        //! constructor stores reference to a grid
        PatchGridGlobalIdSet() = default;
        PatchGridGlobalIdSet(const GridImp& g) : grid_(&g) {}
        // PatchGridGlobalIdSet(const GridImp& g, const std::vector<TrimmingCurve>& trimmingCurves) : grid_(&g) {}

        //! define the type used for persistent indices
        using IdType = typename Trimmer::TrimmerTraits::GlobalIdSetId;

        //! get id of an entity
        /*
           We use the remove_const to extract the Type from the mutable class,
           because the const class is not instantiated yet.
         */
        template <int cd>
        IdType id(const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const {
          // Return id of the host entity
          return e.impl().getHostEntity().id();
        }

        //! get id of subEntity
        /*

         */
        IdType subId(const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i,
                     int codim) const {
          // @todo Trim, the sub indeces are wrong!!!
          //  Return sub id of the host entity
          return e.impl().getHostEntity().subId(i, codim);
        }

        /** @todo Should be private */
        void update() {}

       public:
        // using IndexVariant = IdType;

        // auto getStableIndex(IdType thirdPartyIndex) -> IdType {
        //   auto it = myIndexMapping.find(thirdPartyIndex);
        //   if (it == myIndexMapping.end()) {
        //     // If not found, generate a new index
        //     IdType newIndex{myIndexMapping.size()};
        //     myIndexMapping[thirdPartyIndex] = newIndex;
        //     return myIndexMapping[thirdPartyIndex];
        //   }
        //   // If found, return the existing index
        //   return it->second;
        // };

        IdType newFreeIndex() { return nextFreeIndex_++; }

        // std::map<IdType, IdType> myIndexMapping;
        IdType nextFreeIndex_;

        // store

        const GridImp* grid_;
      };
    }  // namespace DefaultTrim
  }    // namespace IGANEW
}  // namespace Dune

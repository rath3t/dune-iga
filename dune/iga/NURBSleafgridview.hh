// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/common/tuplevector.hh>
#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSgridentity.hh>
#include <dune/iga/NURBSgridleafiterator.hh>
#include <dune/iga/NURBSgridindexsets.hh>


namespace Dune::IGA {


    /** \brief Collect several types associated to OneDGrid LeafGridViews */
    template<class GridImp>
    struct NurbsLeafGridViewTraits {
        typedef NurbsLeafGridViewTraits<GridImp> GridViewImp;

        /** \brief type of the grid */
        typedef typename std::remove_const<GridImp>::type Grid;

//          /** \brief type of the index set */
//          typedef typename GridImpl :: Traits :: LeafIndexSet IndexSet;
//
//          /** \brief type of the intersection */
//          typedef typename GridImpl :: Traits :: LeafIntersection Intersection;
//
//          /** \brief type of the intersection iterator */
//          typedef typename GridImpl :: Traits :: LeafIntersectionIterator IntersectionIterator;
//
//          /** \brief type of the collective communication */
//          typedef typename GridImpl :: Traits :: CollectiveCommunication CollectiveCommunication;

        template<int cd>
        struct Codim {
            typedef typename Grid::Traits ::template Codim<cd>::template Partition<All_Partition>::LeafIterator
                    Iterator;

            typedef typename Grid::Traits::template Codim<cd>::Entity Entity;

            typedef typename Grid::Traits::template Codim<cd>::Geometry Geometry;
            typedef typename Grid::Traits::template Codim<cd>::Geometry LocalGeometry;

            /** \brief Define types needed to iterate over entities of a given partition type */
            template<PartitionIteratorType pit>
            struct Partition {
                /** \brief iterator over a given codim and partition type */
                typedef typename Grid::Traits::template Codim<cd>
                ::template Partition<pit>::LeafIterator
                        Iterator;
            };
        };

        enum {
            conforming = true
        };
    };


    template<typename GridView, int... codim>
    std::tuple<std::vector<NURBSGridEntity<codim, GridView>>...> gridEntityTupleGenerator(
            std::integer_sequence<int, codim...>);

    template<typename GridImpl>
    auto elements(const NURBSLeafGridView<GridImpl> &gridLeafView);

    /** \brief NURBS grid manager */
    template<typename GridImpl>
    class NURBSLeafGridView {
    public:

        using NurbsGridLinearAlgebraTraits = typename GridImpl::NurbsGridLinearAlgebraTraits;
        using GlobalCoordinateType = typename GridImpl::GlobalCoordinateType;
        using LocalCoordinateType = typename GridImpl::LocalCoordinateType;
        using JacobianTransposedType = typename GridImpl::JacobianTransposedType;
        using JacobianInverseTransposed = typename GridImpl::JacobianInverseTransposed;

        using ControlPointNetType = typename GridImpl::ControlPointNetType;

        template<int codim, class GridViewImp>
        friend
        class NURBSGridEntity;

        using Traits = NurbsLeafGridViewTraits<GridImpl>;

        using ctype = double;
        static constexpr int dimension = GridImpl::dimension;
        static constexpr int dimensionworld = GridImpl::dimensionworld;

        template<typename NURBSEntity>
        friend
        class NURBSGridLeafIterator;

        using Grid = typename Traits::Grid;
        typedef NURBSLeafGridView<GridImpl> NURBSGridView;
        typedef NURBSGeometry<dimension, dimensionworld,typename GridImpl::NurbsGridLinearAlgebraTraits> Geometry;
        typedef NURBSGridLeafIndexSet<NURBSGridView> IndexSet;

        template<int cd>
        struct Codim : public Traits::template Codim<cd> {
        };

        /** \brief  constructor
         *
         *  \param[in] knotSpans vector of knotSpans for each dimension
         *  \param[in] controlPoints a n-dimensional net of control points
         *  \param[in] order order of the B-Spline structure for each dimension
         */

        NURBSLeafGridView(const NURBSPatchData<dimension, dimensionworld,NurbsGridLinearAlgebraTraits> &patchData, const Grid &grid)
                : NURBSLeafGridView(patchData.getKnots(), patchData.getControlPoints(),
                                    patchData.getOrder(), grid) {
        }

        NURBSLeafGridView(const std::array<std::vector<double>, dimension> &knotSpans,
                          const ControlPointNetType& controlPoints,
                          const std::array<int, dimension> order, const Grid &grid)
                : NURBSpatch_(
                std::make_shared<NURBSPatch<dimension, dimensionworld,NurbsGridLinearAlgebraTraits>>(knotSpans, controlPoints, order)),
                  grid_{&grid}, entityVector_{std::make_shared<decltype(gridEntityTupleGenerator<NURBSLeafGridView>(
                        std::make_integer_sequence<int, dimension + 1>()))>()}
                        ,indexSet_{*this}
                        {
            int elementSize = NURBSpatch_->size(0);

            Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<1>()), [&](const auto i) {
                std::get<i>(*entityVector_.get()).reserve(elementSize);
                for (unsigned int j = 0; j < NURBSpatch_->size(i); ++j)

                    std::get<i>(*entityVector_.get()).emplace_back(*this, j);

            });
        }

        template<int codim>
        typename Codim<codim>::Entity &getEntity(unsigned int directIndex) const {
            //need to be rewrite for other codims
            if (codim == 0) {
                return *(entityVector_.get()[Indices::_0]->at(directIndex));
            }

        }

        /** \brief obtain collective communication object */
        const auto &comm() const {
            return grid().comm();
        }

        /** \brief obtain collective communication object */
        const auto &grid() const {
            return *grid_;
        }

        template<int cd>
        typename Codim<cd>::Iterator begin() const {
            return typename Codim<cd>::Iterator(std::get<cd>(*entityVector_.get()).begin());
        }

        template<int cd>
        typename Codim<cd>::Iterator end() const {
            return typename Codim<cd>::Iterator(std::get<cd>(*entityVector_.get()).end());
        }


        template<int cd,PartitionIteratorType piType>
        typename Codim<cd>::template Partition<piType>::Iterator begin() const {
            if (piType!=Ghost_Partition)
            return typename Codim<cd>::template Partition<piType>::Iterator(this->template begin<cd>());
            else
                return typename Codim<cd>::template Partition<piType>::Iterator(this->template end<cd>());
        }

        template<int cd,PartitionIteratorType piType>
        typename Codim<cd>::template Partition<piType>::Iterator end() const {
            return typename Codim<cd>::template Partition<piType>::Iterator(this->template end<cd>());
        }

        const IndexSet& indexSet() const {
            return indexSet_;
        }

        auto size(int codim) const {
            if (codim == 0)
                return std::get<0>(*entityVector_.get()).size();
            else if (codim == 1)
                return std::get<1>(*entityVector_.get()).size();
            if constexpr(dimension > 1)
                if (codim == 2)
                    return std::get<2>(*entityVector_.get()).size();
            if constexpr(dimension > 2)
                if (codim == 3)
                    return std::get<3>(*entityVector_.get()).size();
        }

    private:
        friend auto elements<GridImpl>(const NURBSLeafGridView<GridImpl> &gridLeafView);
        std::shared_ptr<NURBSPatch<dimension, dimensionworld,NurbsGridLinearAlgebraTraits>> NURBSpatch_;
        NURBSGridLeafIndexSet<NURBSGridView> indexSet_;
        const Grid * grid_;
        using EntityVectorType = decltype(gridEntityTupleGenerator<NURBSLeafGridView>(
                std::make_integer_sequence<int, dimension + 1>()));
        std::shared_ptr<EntityVectorType> entityVector_{};
    };

    template<typename GridImpl>
    auto elements(const NURBSLeafGridView<GridImpl> &gridLeafView) {
        return std::get<0>(*gridLeafView.entityVector_.get());
    }
}



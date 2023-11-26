// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <ranges>

#include "hierachicpatchgridlocalgeometry.hh"
namespace Dune {
  namespace Geo {

    /** \class TrimmedReferenceElement
     *  \ingroup GeometryTrimmedReferenceElements
     *  \brief This class provides access to geometric and topological
     *  properties of a reference element.
     *
     *  This includes its type,
     *  the number of subentities, the volume, and a method for checking
     *  if a point is contained in the reference element.
     *  The embedding of each subentity into the reference element is also
     *  provided.
     *
     *  This class has value semantics, i.e. you can (and should) pass it
     *  around by value and not by reference and store a copy of it.
     *
     *  Instances of this object for a given geometry type can be retrieved
     *  from the TrimmedReferenceElements class.
     *
     */
    template <typename GridentityImpl>
    class TrimmedReferenceElement {
     public:
      //! The dimension of the reference element.
      static constexpr int dimension = GridentityImpl::dimension;

      /** \brief Collection of types depending on the codimension */
      template <int codim>
      struct Codim {
        //! type of geometry embedding a subentity into the reference element
        using TrimmedLocalGeometry
            = IGANEW::TrimmedPatchGridLocalGeometry<dimension - codim, dimension,
                                                    typename GridentityImpl::Grid::Implementation>;
        using UnTrimmedLocalGeometry = IGANEW::PatchGridLocalGeometry<dimension - codim, dimension,
                                                                      typename GridentityImpl::Grid::Implementation>;
        using Geometry               = std::conditional_t<codim == 0, TrimmedLocalGeometry, UnTrimmedLocalGeometry>;
      };

      //! The coordinate field type.
      using ctype = typename GridentityImpl::Implementation::ctype;

      //! The coordinate field type.
      using CoordinateField = ctype;

      //! The coordinate type.
      using Coordinate = Dune::FieldVector<ctype, dimension>;

      /** \brief Type used for volume */
      typedef ctype Volume;

      /** \brief number of subentities of codimension c
       *
       *  \param[in]  c  codimension whose size is desired
       */
      int size(int c) const {
        // TODO Trim
        return {};
      }

      /** \brief number of subentities of codimension cc of subentity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the number of subentities
       *  of codimension cc of the current reference element, that are also
       *  a subentity of E. If cc<c this number is zero.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E (0 <= c <= dim)
       *  \param[in]  cc  codimension whose size is desired (0 <= cc <= dim)
       */
      int size(int i, int c, int cc) const {
        // TODO Trim
        return {};
      }

      /** \brief obtain number of ii-th subentity with codim cc of (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. And denote by S the ii-th subentity of codimension
       *  (cc-c) of E. Then, S is a also a subentity of codimension cc of the current
       *  reference element. This method returns the number of S with respect
       *  to the current reference element.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       *  \param[in]  ii  number of subentity S (with respect to E)
       *  \param[in]  cc  codimension of subentity S (c <= cc <= dim)
       */
      int subEntity(int i, int c, int ii, int cc) const {
        // TODO Trim
        return {};
      }

      /** \brief Obtain the range of numbers of subentities with codim cc of (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns a range of numbers of
       *  all subentities of E with codimension cc. Notice that the sub-subentity
       *  codimension as well as the numbers in the returned range are
       *  given with respect to the reference element itself and not with
       *  respect to E. For 0<=cc<c this will return an empty range.
       *  The returned range r provide the methods r.begin(), r.end(),
       *  r.contains(std::size_t) and r.size() mimicking an immutable
       *  iterable set.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       *  \param[in]  cc  codimension of subentity S (0 <= cc <= dim)
       *
       *  \returns An iterable range of numbers of the sub-subentities.
       */
      auto subEntities(int i, int c, int cc) const {
        // TODO TRIM this should return something usefull
        return std::ranges::iota_view(0, 10);
      }

      /** \brief obtain the type of subentity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the GeometryType of E.
       *
       *  \param[in]  i      number of subentity E (0 <= i < size( c ))
       *  \param[in]  c      codimension of subentity E
       */
      GeometryType type(int i, int c) const {
        // TODO This method makes only sense for the 2D trimming case, since for 3D the facets subentities could also be
        // none
        if (c == 0)
          return GeometryTypes::none(dimension);
        else
          return GeometryTypes::cube(dimension);
      }

      /** \brief obtain the type of this reference element
      Since it is a trimmed element we basically only return none here as the most general case
       */
      GeometryType type() const {
        // TODO for some cases we could also return triangle or something else, but im not sure if
        //  this is too complicated and also unnecessary
        return GeometryTypes::none(dimension);
      }

      /** \brief position of the barycenter of entity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the coordinates of
       *  the center of gravity of E within the current reference element.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       */
      Coordinate position(int i, int c) const {
        // TODO this functions could be a bit complicated
        //  we have to implement https://en.wikipedia.org/wiki/Center_of_mass#A_continuous_volume
        //  M as the volume and rho(R)=1
        return {};
      }

      /** \brief check if a coordinate is in the reference element
       *
       *  This method returns true if the given local coordinate is within this
       *  reference element.
       *
       *  \param[in]  local  coordinates of the point
       */
      bool checkInside(const Coordinate& local) const {
        // TODO this functions could be a bit complicated basically we have to make sure the point lies inside the outer
        // boundary loop, but outside the inner loops thus we have to implement something as
        // https://en.wikipedia.org/wiki/Point_in_polygon#:~:text=One%20simple%20way%20of%20finding,an%20even%20number%20of%20times.
        //  maybe what we are searching for is already existing in Clipperlib
        //  https://angusj.com/clipper2/Docs/Units/Clipper/Functions/PointInPolygon.htm looks promising

        return {};
      }

      /** \brief obtain the embedding of subentity (i,codim) into the reference
       *         element
       *
       *  Denote by E the i-th subentity of codimension codim of the current
       *  reference element. This method returns a \ref Dune::AffineGeometry
       *  that maps the reference element of E into the current reference element.
       *
       *  \tparam     codim  codimension of subentity E
       *
       *  \param[in]  i      number of subentity E (0 <= i < size( codim ))
       */
      template <int codim>
      typename Codim<codim>::Geometry geometry(int i) const {
        // TODO trim returns the reference element in geometry space
        //  return _impl->template geometry<codim>(i);
        //  if constexpr (codim==0)
        //    return IGANEW::TrimmedPatchGridLocalGeometry();
        //  else
        //    return IGANEW::PatchGridLocalGeometry
      }

      /** \brief obtain the volume of the reference element */
      CoordinateField volume() const {
        // TODO trim, integrate on the trimmed patch
        return {};
      }

      /** \brief obtain the integration outer normal of the reference element
       *
       *  The integration outer normal is the outer normal whose length coincides
       *  with the face's integration element.
       *
       *  \param[in]  face  index of the face, whose normal is desired
       */
      Coordinate integrationOuterNormal(int face) const {
        // TODO compute tangent of the curve and compute by cross-product the outword normal, only 2D
        return {};
      }

      /** \brief Constructs an empty reference element.
       *
       * This constructor creates an empty (invalid) reference element. This element may not be
       * used in any way except for assigning other reference elements to it. After
       * assigning a valid reference element (obtained from TrimmedReferenceElements), it may
       * be used without restrictions.
       */
      TrimmedReferenceElement() : entity{nullptr} {}

      explicit TrimmedReferenceElement(const GridentityImpl& ent) : entity{&ent} {}

      //! Compares for equality with another reference element.
      bool operator==(const TrimmedReferenceElement& r) const {
        // TODO, just check if the entities cooincide
        return entity == r.entity;
      }

      //! Compares for inequality with another reference element.
      bool operator!=(const TrimmedReferenceElement& r) const { return not(*this == r); }

      //! Yields a hash value suitable for storing the reference element a in hash table
      friend std::size_t hash_value(const TrimmedReferenceElement& r) {
        // TODO, this is not needed maybe
        return {};
      }

     private:
      // The implementation must be a friend to construct a wrapper around itself.
      const GridentityImpl* entity;
      // TODO mayby store here all the trimming information anyway?
      // But this should have value semantics and therefore it should be cheap to copy, thus maybe store it at the
      // entity
    };
  }  // namespace Geo
}  // namespace Dune

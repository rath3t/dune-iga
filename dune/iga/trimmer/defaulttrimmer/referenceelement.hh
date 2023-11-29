// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <ranges>

#include "trimmedlocalgeometry.hh"
#include "trimmer.hh"

namespace Dune {
  namespace IGANEW {
    namespace DefaultTrim {

      template <int mydim_, typename ScalarType>
      struct ElementTrimData;

      /** \class DefaultTrimmedReferenceElement
       *  \ingroup GeometryTrimmedReferenceElements
       *  @brief This class provides access to geometric and topological
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
      template <int dim, typename ct>
      class DefaultTrimmedReferenceElement {
       public:
        //! The dimension of the reference element.
        static constexpr int mydimension = dim;
        //! The coordinate field type.
        using ctype                         = ct;
        using ParameterSpaceGrid            = YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>;
        using TrimDataType                  = ElementTrimData<mydimension, ctype>;
        using TrimDataTypeOptionalReference = std::optional<std::reference_wrapper<const TrimDataType>>;

        /** @brief Collection of types depending on the codimension */
        template <int codim>
        struct Codim {
          //! type of geometry embedding a subentity into the reference element

          using UnTrimmedGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;
          using Geometry          = TrimmedPatchLocalGeometry<mydimension - codim, mydimension, ctype>;
        };

        //! The coordinate field type.
        using CoordinateField = ctype;

        //! The coordinate type.
        using Coordinate = Dune::FieldVector<ctype, mydimension>;

        /** @brief Type used for volume */
        typedef ctype Volume;

        /** @brief number of subentities of codimension c
         *
         *  @param[in]  c  codimension whose size is desired
         */
        int size(int c) const {
          // TODO Trim
          if (not trimData_) return cubeGeometry.size(c);
          assert(false);
        }

        /** @brief number of subentities of codimension cc of subentity (i,c)
         *
         *  Denote by E the i-th subentity of codimension c of the current
         *  reference element. This method returns the number of subentities
         *  of codimension cc of the current reference element, that are also
         *  a subentity of E. If cc<c this number is zero.
         *
         *  @param[in]  i   number of subentity E (0 <= i < size( c ))
         *  @param[in]  c   codimension of subentity E (0 <= c <= dim)
         *  @param[in]  cc  codimension whose size is desired (0 <= cc <= dim)
         */
        int size(int i, int c, int cc) const {
          // TODO Trim
          if (not trimData_) return cubeGeometry.size(i, c, cc);
          assert(false);
        }

        /** @brief obtain number of ii-th subentity with codim cc of (i,c)
         *
         *  Denote by E the i-th subentity of codimension c of the current
         *  reference element. And denote by S the ii-th subentity of codimension
         *  (cc-c) of E. Then, S is a also a subentity of codimension cc of the current
         *  reference element. This method returns the number of S with respect
         *  to the current reference element.
         *
         *  @param[in]  i   number of subentity E (0 <= i < size( c ))
         *  @param[in]  c   codimension of subentity E
         *  @param[in]  ii  number of subentity S (with respect to E)
         *  @param[in]  cc  codimension of subentity S (c <= cc <= dim)
         */
        int subEntity(int i, int c, int ii, int cc) const {
          // TODO Trim
          if (not trimData_) return cubeGeometry.subEntity(i, c, ii, cc);
          assert(false);
        }

        /** @brief Obtain the range of numbers of subentities with codim cc of (i,c)
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
         *  @param[in]  i   number of subentity E (0 <= i < size( c ))
         *  @param[in]  c   codimension of subentity E
         *  @param[in]  cc  codimension of subentity S (0 <= cc <= dim)
         *
         *  \returns An iterable range of numbers of the sub-subentities.
         */
        auto subEntities(int i, int c, int cc) const {
          // TODO TRIM this should return something usefull
          if (not trimData_) return cubeGeometry.subEntities(i, c, cc);
          assert(false);
        }

        /** @brief obtain the type of subentity (i,c)
         *
         *  Denote by E the i-th subentity of codimension c of the current
         *  reference element. This method returns the GeometryType of E.
         *
         *  @param[in]  i      number of subentity E (0 <= i < size( c ))
         *  @param[in]  c      codimension of subentity E
         */
        GeometryType type(int i, int c) const {
          // TODO This method makes only sense for the 2D trimming case, since for 3D the facets subentities could also
          // be none
          if (c == 0 and not trimData_) return GeometryTypes::none(mydimension);

          if (trimData_) assert(false);

          return GeometryTypes::cube(mydimension);
        }

        /** @brief obtain the type of this reference element
        Since it is a trimmed element we basically only return none here as the most general case
         */
        GeometryType type() const {
          // TODO for some cases we could also return triangle or something else, but im not sure if
          //  this is too complicated and also unnecessary
          if (not trimData_) return GeometryTypes::cube(mydimension);

          assert(false);
          return GeometryTypes::none(mydimension);
        }

        /** @brief position of the barycenter of entity (i,c)
         *
         *  Denote by E the i-th subentity of codimension c of the current
         *  reference element. This method returns the coordinates of
         *  the center of gravity of E within the current reference element.
         *
         *  @param[in]  i   number of subentity E (0 <= i < size( c ))
         *  @param[in]  c   codimension of subentity E
         */
        Coordinate position(int i, int c) const {
          // TODO this functions could be a bit complicated
          //  we have to implement https://en.wikipedia.org/wiki/Center_of_mass#A_continuous_volume
          //  M as the volume and rho(R)=1
          if (not trimData_) return cubeGeometry.position(i, c);
          assert(false);
        }

        /** @brief check if a coordinate is in the reference element
         *
         *  This method returns true if the given local coordinate is within this
         *  reference element.
         *
         *  @param[in]  local  coordinates of the point
         */
        bool checkInside(const Coordinate& local) const {
          if (not trimData_) cubeGeometry.checkInside(local);
          return trimData_->checkInside(local);
        }

        /** @brief obtain the embedding of subentity (i,codim) into the reference
         *         element
         *
         *  Denote by E the i-th subentity of codimension codim of the current
         *  reference element. This method returns a \ref Dune::AffineGeometry
         *  that maps the reference element of E into the current reference element.
         *
         *  @tparam     codim  codimension of subentity E
         *
         *  @param[in]  i      number of subentity E (0 <= i < size( codim ))
         */
        template <int codim>
        typename Codim<codim>::Geometry geometry(int i) const {
          // TODO trim returns the reference element in geometry space
          //  return _impl->template geometry<codim>(i);
          //  if constexpr (codim==0)
          //    return IGANEW::TrimmedPatchGridLocalGeometry();
          //  else
          //    return IGANEW::PatchGridLocalGeometry
          return typename Codim<codim>::Geometry(trimData_);
        }

        /** @brief obtain the volume of the reference element */
        CoordinateField volume() const {
          // TODO trim, integrate on the trimmed patch
          if (not trimData_) return cubeGeometry.volume();
        }

        /** @brief obtain the integration outer normal of the reference element
         *
         *  The integration outer normal is the outer normal whose length coincides
         *  with the face's integration element.
         *
         *  @param[in]  face  index of the face, whose normal is desired
         */
        Coordinate integrationOuterNormal(int face) const {
          // TODO compute tangent of the curve and compute by cross-product the outword normal, only 2D
          return cubeGeometry.integrationOuterNormal(face);
        }

        /** @brief Constructs an empty reference element.
         *
         * This constructor creates an empty (invalid) reference element. This element may not be
         * used in any way except for assigning other reference elements to it. After
         * assigning a valid reference element (obtained from TrimmedReferenceElements), it may
         * be used without restrictions.
         */
        DefaultTrimmedReferenceElement() = default;

        explicit DefaultTrimmedReferenceElement(TrimDataTypeOptionalReference trimData) : trimData_{trimData} {}

        //! Compares for equality with another reference element.
        bool operator==(const DefaultTrimmedReferenceElement& r) const {
          // TODO, just check if the triangulations cooincide
          return trimData_ == r.trimData_;
        }

        //! Compares for inequality with another reference element.
        bool operator!=(const DefaultTrimmedReferenceElement& r) const { return not(*this == r); }

        //! Yields a hash value suitable for storing the reference element a in hash table
        friend std::size_t hash_value(const DefaultTrimmedReferenceElement& r) {
          // TODO, this is not needed maybe
          return hash_value(ReferenceElements<ctype, mydimension>::cube());
          return {};
        }

       private:
        TrimDataTypeOptionalReference trimData_;
        const typename ReferenceElements<ctype, mydimension>::ReferenceElement& cubeGeometry{
            ReferenceElements<ctype, mydimension>::cube()};
        // TODO mayby store here all the trimming information anyway?
        // But this should have value semantics and therefore it should be cheap to copy, thus maybe store it at the
        // entity
      };
    }  // namespace DefaultTrim
  }    // namespace IGANEW
}  // namespace Dune

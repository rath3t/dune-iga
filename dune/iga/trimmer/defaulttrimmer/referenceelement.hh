// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include "elementtrimdata.hh"
#include "trimmedlocalgeometry.hh"

#include <ranges>

#include <dune/common/reservedvector.hh>

namespace Dune {
namespace IGANEW {
  namespace DefaultTrim {
    // template <int mydim_, typename ScalarType>
    // struct ElementTrimData;

    // template <int mydim_, typename ScalarType>
    // struct ElementTrimData;

    /** \class TrimmedReferenceElement
     *  @ingroup GeometryTrimmedReferenceElements
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
    template <int dim, typename GridImp>
    class TrimmedReferenceElement
    {
    public:
      // The dimension of the reference element.
      static constexpr int mydimension = dim;
      static constexpr int dimension   = mydimension;
      // The coordinate field type.
      using Trimmer            = typename GridImp::Trimmer;
      using ctype              = typename Trimmer::ctype;
      using ParameterSpaceGrid = YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>;
      using TrimDataType       = typename GridImp::Trimmer::TrimmerTraits::ElementTrimData;

      /** @brief Collection of types depending on the codimension */
      template <int codim>
      struct Codim
      {
        // type of geometry embedding a subentity into the reference element

        using Geometry = TrimmedLocalGeometryImpl<mydimension - codim, mydimension, const GridImp,
                                                  LocalGeometryTag::InReferenceElement>;
      };

      // The coordinate field type.
      using CoordinateField = ctype;

      // The coordinate type.
      using Coordinate = Dune::FieldVector<ctype, mydimension>;

      /** @brief Type used for volume */
      typedef ctype Volume;

      /** @brief number of subentities of codimension c
       *
       *  @param[in]  c  codimension whose size is desired
       */
      int size(int codim) const {
        if (trimData_.flag() == ElementTrimFlag::full)
          return cubeGeometry.size(codim);
        if (trimData_.flag() == ElementTrimFlag::trimmed) {
          switch (codim) {
            case 0:
              return 1;
            case 1:
              return trimData_.edges().size();
            case 2:
              return trimData_.vertices().size();
            default:
              assert(false && "Wrong codim requested");
          }
        }
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
        // auto subEntityName = [](std::size_t codim)-> std::string {
        //   return codim==0 ? "element" : codim==1 ? "edges" :  "vertices" ;
        // };
        // std::cout<<"The size of  "<<subEntityName(cc)<<" with codim "<<cc<<" of the "<<i<<"-th
        // "<<subEntityName(c)<<" of codim "<<c<<" is requested."<< " i: "<<i<<" c: "<<c<<" cc: "<<cc<<std::endl;
        //
        // std::cout<<"size i: "<<i<<" c: "<<c<<" cc: "<<cc<<std::endl;
        if (trimData_.flag() == ElementTrimFlag::full)
          return cubeGeometry.size(i, c, cc);

        if (trimData_.flag() == ElementTrimFlag::trimmed) {
          if (cc < c)
            return 0;
          else if (c == cc) // Number ofs element of the element 1
          {
            return 1;
          } else if (i == 0 and c == 0 and cc == 1) // Number of edges of the element
          {
            return edgeIndices.size();
          } else if (i == 0 and c == 0 and cc == 2) // Number of vertices of the element
          {
            return trimData_.vertices().size();
          } else if (c == 1 and cc == 2) // Number of vertices of an edge
          {
            return 2;
          }
        }
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
        // auto subEntityName = [](std::size_t codim)-> std::string {
        //   return codim==0 ? "element" : codim==1 ? "edges" :  "vertex" ;
        // };
        // std::cout<<"The index of the "<<ii<<"-th "<<subEntityName(cc)<<" with codim "<<cc<<" of the "<<i<<"-th
        // "<<subEntityName(c)<<" of codim "<<c<<" is requested."<< " i: "<<i<<" c: "<<c<<" ii: "<<ii<<" cc:
        // "<<cc<<std::endl;
        //
        if (trimData_.flag() == ElementTrimFlag::full)
          return cubeGeometry.subEntity(i, c, ii, cc);
        if (trimData_.flag() == ElementTrimFlag::trimmed) {
          if (cc < c)
            return 0;
          else if (i == 0 and c == cc) // Index of the element of the element is a single 0
          {
            return 0;
          } else if (i == 0 and c == 0 and cc == 1) // Get index edge of the element
          {
            return ii;
          } else if (i == 0 and c == 0 and cc == 2) // Index of vertex of the element
          {
            return ii;
          } else if (c == 1 and cc == 2) // Index of vertex of an edge
          {
            return ii == 0 ? i : (i + 1 < edgeIndices.size() ? i + 1 : 0);
            // if we are the left vertex the index is the same as the edge index,
            // if we are the right vertex the index is the edge index +1
            // if we are the one after the last vertex, we are the first vertex (with index 0)
          } else if (c == cc) // Index of entity of the entity
          {
            return i;
          }
        }
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
        // auto subEntityName = [](std::size_t codim)-> std::string {
        //   return codim==0 ? "element" : codim==1 ? "edges" :  "vertices" ;
        // };
        // std::cout<<"The range numbers of "<<subEntityName(cc)<<" with codim "<<cc<<" of the "<<i<<"-th
        // "<<subEntityName(c)<<" of codim "<<c<<" is requested."<<" i: "<<i<<" c: "<<c<<" cc: "<<cc<<std::endl;

        if (trimData_.flag() == ElementTrimFlag::full) {
          auto range = cubeGeometry.subEntities(i, c, cc);
          return SubEntityRangeImpl(range.begin(), range.end(), true, i, c, cc);
        }
        if (trimData_.flag() == ElementTrimFlag::trimmed) {
          if (cc < c)
            return SubEntityRangeImpl(dummyArray.begin(), dummyArray.end());
          else if (i == 0 and c == 0 and
                   cc == 0) // iterate over the entity itself (only one element in storage with value zero)
          {
            SubEntityRangeImpl range(elementIndexDummy.begin(), elementIndexDummy.end());
            return range;
          } else

              if (i == 0 and c == 0 and cc == 1) // iterate over the edges of the element
          {
            SubEntityRangeImpl range(edgeIndices.begin(), edgeIndices.end());
            return range;
          } else if (c == 0 and cc == 2) // iterate over the vertices of the element
          {
            SubEntityRangeImpl range(verticesIndices.begin(), verticesIndices.end() - 1);
            return range;
          }

          else if (c == 1 and cc == 2) // iterate over the vertices of edge
          {
            SubEntityRangeImpl range(verticesIndices.begin() + i, verticesIndices.begin() + i + 2);
            return range;
          }

          else if (c == cc) // iterate over the vertices of edge
          {
            if (c == 2) {
              SubEntityRangeImpl range(verticesIndices.begin() + i, verticesIndices.begin() + i + 1);
              return range;
            } else if (cc == 1) {
              SubEntityRangeImpl range(edgeIndices.begin() + i, edgeIndices.begin() + i + 1);
              return range;
            }
          }
        }

        assert(false);
      }

      class SubEntityRangeImpl : public Dune::IteratorRange<const unsigned int*>
      {
        using Base = typename Dune::IteratorRange<const unsigned int*>;

      public:
        using iterator       = Base::iterator;
        using const_iterator = Base::const_iterator;

        SubEntityRangeImpl(const iterator& begin, const iterator& end, bool host, int i, int c, int cc)
            : Base(begin, end),
              host_(host),
              size_(end - begin),
              i_{i},
              c_{c},
              cc_{cc}

        {}

        SubEntityRangeImpl(const iterator& begin, const iterator& end)
            : Base(begin, end),
              size_(end - begin) {}

        SubEntityRangeImpl()
            : size_(0) {}

        std::size_t size() const {
          return size_;
        }

        bool contains(std::size_t i) const {
          if (host_)
            return cubeGeometry.subEntities(i_, c_, cc_).contains(i);
          else {
            if (size_ == 0)
              return false;
            else if (size_ == 2) // two vertices at one edge case
            {
              const int back  = *(end() - 1);
              const int front = *begin();
              return i == front or i == back;
            }
            return size_ == 0 ? false : i >= *begin() and i <= *(end() - 1); // general case
          }
        }

      private:
        const typename ReferenceElements<ctype, mydimension>::ReferenceElement& cubeGeometry{
            ReferenceElements<ctype, mydimension>::cube()};
        std::size_t size_;
        bool host_{};
        int i_{};
        int c_{};
        int cc_{};
      };

      /** @brief obtain the type of subentity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the GeometryType of E.
       *
       *  @param[in]  i      number of subentity E (0 <= i < size( c ))
       *  @param[in]  c      codimension of subentity E
       */
      GeometryType type(int i, int c) const {
        // @todo This method makes only sense for the 2D trimming case, since for 3D the facets subentities could also
        // be none
        if (c == 0 and not trimData_.flag() != ElementTrimFlag::full)
          return GeometryTypes::none(mydimension);
        if (c == 0 and trimData_.flag() == ElementTrimFlag::full)
          return GeometryTypes::cube(mydimension);

        return GeometryTypes::cube(mydimension);
      }

      /** @brief obtain the type of this reference element
      Since it is a trimmed element we basically only return none here as the most general case
       */
      GeometryType type() const {
        // @todo for some cases we could also return triangle or something else, but im not sure if
        //  this is too complicated and also unnecessary
        if (trimData_.flag() == ElementTrimFlag::full)
          return GeometryTypes::cube(mydimension);

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
        // @todo this functions could be a bit complicated
        //  we have to implement https://en.wikipedia.org/wiki/Center_of_mass#A_continuous_volume
        //  M as the volume and rho(R)=1
        if (not trimData_.flag() == ElementTrimFlag::full)
          return cubeGeometry.position(i, c);

        assert(false);
      }

      /** @brief check if a coordinate is in the reference element
       *
       *  This method returns true if the given local coordinate is within this
       *  reference element.
       *
       *  @param[in]  local  coordinates of the point

        // @todo this functions could be a bit complicated basically we have to make sure the point lies inside the
        // outer boundary loop, but outside the inner loops thus we have to implement something as
        //
       https://en.wikipedia.org/wiki/Point_in_polygon#:~:text=One%20simple%20way%20of%20finding,an%20even%20number%20of%20times.
        //  maybe what we are searching for is already existing in Clipperlib
        //  https://angusj.com/clipper2/Docs/Units/Clipper/Functions/PointInPolygon.htm looks promising
        return true;
*/
      bool checkInside(const Coordinate& local) const {
        if (trimData_.flag() == ElementTrimFlag::full)
          return cubeGeometry.checkInside(local);
        return trimData_.checkInside(local);
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
        // @todo trim returns the reference element in geometry space
        //  return _impl->template geometry<codim>(i);
        //  if constexpr (codim==0)
        //    return IGANEW::TrimmedPatchGridLocalGeometry();
        //  else
        //    return IGANEW::PatchGridLocalGeometry
        return typename Codim<codim>::Geometry(trimData_);
      }

      /** @brief obtain the volume of the reference element */
      CoordinateField volume() const {
        // @todo trim, integrate on the trimmed patch
        if (trimData_.flag() == ElementTrimFlag::full)
          return cubeGeometry.volume();
        return trimData_.volume();
      }

      /** @brief obtain the integration outer normal of the reference element
       *
       *  The integration outer normal is the outer normal whose length coincides
       *  with the face's integration element.
       *
       *  @param[in]  face  index of the face, whose normal is desired
       */
      Coordinate integrationOuterNormal(int face) const {
        // @todo compute tangent of the curve and compute by cross-product the outword normal, only 2D
        return cubeGeometry.integrationOuterNormal(face);
      }

      /** @brief Constructs an empty reference element.
       *
       * This constructor creates an empty (invalid) reference element. This element may not be
       * used in any way except for assigning other reference elements to it. After
       * assigning a valid reference element (obtained from TrimmedReferenceElements), it may
       * be used without restrictions.
       */
      TrimmedReferenceElement() = default;

      explicit TrimmedReferenceElement(const TrimDataType& trimData)
          : trimData_{trimData} {
        if (trimData_.flag() == ElementTrimFlag::trimmed) {
          edgeIndices.resize(trimData_.edges().size());
          for (auto eI : Dune::range(trimData_.edges().size()))
            edgeIndices[eI] = eI;

          verticesIndices.resize(trimData_.vertices().size() + 1);
          for (auto vI : Dune::range(trimData_.vertices().size())) {
            verticesIndices[vI] = vI;
          }
          // to mimic a circular vector we add the first vertex at the end
          verticesIndices.back() = 0;
        }
      }

      // Compares for equality with another reference element.
      bool operator==(const TrimmedReferenceElement& r) const {
        // @todo, just check if the triangulations cooincide
        return trimData_ == r.trimData_;
      }

      // Compares for inequality with another reference element.
      bool operator!=(const TrimmedReferenceElement& r) const {
        return not(*this == r);
      }

      // Yields a hash value suitable for storing the reference element a in hash table
      friend std::size_t hash_value(const TrimmedReferenceElement& r) {
        // @todo, this is not needed maybe
        return hash_value(ReferenceElements<ctype, mydimension>::cube());
        return {};
      }

    private:
      Dune::ReservedVector<unsigned int, 32> edgeIndices;
      Dune::ReservedVector<unsigned int, 32> verticesIndices;
      std::array<unsigned int, 1> elementIndexDummy{{0}};
      std::array<unsigned int, 2> verteEdgeIndexDummy{
          {0, 1}
      };
      std::array<unsigned int, 0> dummyArray;
      ElementTrimDataImpl<GridImp> trimData_;
      const typename ReferenceElements<ctype, mydimension>::ReferenceElement& cubeGeometry{
          ReferenceElements<ctype, mydimension>::cube()};
      // @todo mayby store here all the trimming information anyway?
      // But this should have value semantics and therefore it should be cheap to copy, thus maybe store it at the
      // entity
    };
  } // namespace DefaultTrim
} // namespace IGANEW
} // namespace Dune

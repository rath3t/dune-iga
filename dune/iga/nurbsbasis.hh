// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NURBSBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NURBSBASIS_HH

/** \file
 * \brief The B-spline global function space basis
 */

#include <array>
#include <numeric>

#include <dune/iga/bsplinealgorithms.hh>
#include <dune/iga/concepts.hh>
#include <dune/iga/dunelinearalgebratraits.hh>
#include <dune/iga/igaalgorithms.hh>


#include <dune/common/dynmatrix.hh>
#include <dune/common/diagonalmatrix.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune::Functions {

  // A maze of dependencies between the different parts of this.  We need a few forward declarations
  template <typename GV, typename R, typename MI>
  class NurbsLocalFiniteElement;

  template <typename GV, class MI, Dune::IGA::LinearAlgebra NurbsGridLinearAlgebraTraits = Dune::IGA::DuneLinearAlgebraTraits<double>>
  class NurbsPreBasis;

  /** \brief LocalBasis class in the sense of dune-localfunctions, presenting the restriction
   * of a B-spline patch to a knot span
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV Grid view that the basis is defined on
   * \tparam R Number type used for spline function values
   */
  template <class GV, class R, class MI>
  class NurbsLocalBasis {
    friend class NurbsLocalFiniteElement<GV, R, MI>;

    typedef typename GV::ctype D;
    enum { dim = GV::dimension };

  public:
    //! \brief export type traits for function signature
    using Traits = LocalBasisTraits<D, dim, FieldVector<D, dim>, R, 1, FieldVector<R, 1>, FieldMatrix<R, 1, dim>>;

    /** \brief Constructor with a given B-spline patch
     *
     * The patch object does all the work.
     */
    NurbsLocalBasis(const NurbsPreBasis<GV, MI>& preBasis, const NurbsLocalFiniteElement<GV, R, MI>& lFE)
        : preBasis_(preBasis), lFE_(lFE) {}

    /** \brief Evaluate all shape functions
     * \param in Coordinates where to evaluate the functions, in local coordinates of the current knot span
     */
    void evaluateFunction(const FieldVector<D, dim>& in, std::vector<FieldVector<R, 1>>& out) const {
      FieldVector<D, dim> globalIn = offset_;
      scaling_.umv(in, globalIn);

      preBasis_.evaluateFunction(globalIn, out, lFE_.currentKnotSpan_);
    }

    /** \brief Evaluate Jacobian of all shape functions
     * \param in Coordinates where to evaluate the Jacobian, in local coordinates of the current knot span
     */
    void evaluateJacobian(const FieldVector<D, dim>& in, std::vector<FieldMatrix<D, 1, dim>>& out) const {
      FieldVector<D, dim> globalIn = offset_;
      scaling_.umv(in, globalIn);

      preBasis_.evaluateJacobian(globalIn, out, lFE_.currentKnotSpan_);

      for (size_t i = 0; i < out.size(); i++)
        for (int j = 0; j < dim; j++)
          out[i][0][j] *= scaling_[j][j];
    }

    //! \brief Evaluate all shape functions and derivatives of any degree
    inline void partial(const typename std::array<unsigned int, dim>& order, const typename Traits::DomainType& in,
                        std::vector<typename Traits::RangeType>& out) const {
      FieldVector<D, dim> globalIn = offset_;
      scaling_.umv(in, globalIn);
      preBasis_.partial(order, globalIn, out, lFE_.currentKnotSpan_);

      for (size_t d = 0; d < dim; ++d)
        for (size_t i = 0; i < out.size(); i++)
          for (std::size_t fac = 0; fac < order[d]; ++fac)
            out[i][0] *= scaling_[d][d];
    }

    /** \brief Polynomial degree of the shape functions
     *
     * Unfortunately, the general interface of the LocalBasis class mandates that the 'degree' method
     * takes no arguments, and returns a single integer.  It therefore cannot reflect that fact that
     * a B-spline basis function can easily have different orders in the different coordinate directions.
     * We therefore take the conservative choice and return the maximum over the orders of all directions.
     */
    [[nodiscard]] unsigned int order() const {
      return *std::max_element(preBasis_.patchData_.degree.begin(), preBasis_.patchData_.degree.end());
    }

    /** \brief Return the number of basis functions on the current knot span
     */
    [[nodiscard]] std::size_t size() const { return lFE_.size(); }

  private:
    const NurbsPreBasis<GV, MI>& preBasis_;

    const NurbsLocalFiniteElement<GV, R, MI>& lFE_;

    // Coordinates in a single knot span differ from coordinates on the B-spline patch
    // by an affine transformation.  This transformation is stored in offset_ and scaling_.
    FieldVector<D, dim> offset_;
    DiagonalMatrix<D, dim> scaling_;
  };

  /** \brief Attaches a shape function to an entity
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * The attachment uses the same degree as for Qk elements.  This does *not* provide sufficient information
   * to compute global indices for the shape functions.  However, it does allow to find all degrees of freedom
   * that belong to the grid boundary, if the knot vector is open.
   *
   * \note Currently only implemented for 1d and 2d grids.  For higher dimensions you can still use
   *   the BSplineBasis, but you won't be able to find the degrees of freedom on the grid boundary.
   *
   * \tparam dim Dimension of the reference cube
   */
  template <int dim>
  class NurbsLocalCoefficients {
    // Return i as a d-digit number in the (k+1)-nary system
    std::array<unsigned int, dim> multiindex(unsigned int i) const {
      std::array<unsigned int, dim> alpha;
      for (int j = 0; j < dim; j++) {
        alpha[j] = i % sizes_[j];
        i        = i / sizes_[j];
      }
      return alpha;
    }

    /** \brief Set the 'subentity' field for each dof for a 1d element */
    void setup1d(std::vector<unsigned int>& subEntity) {
      if (sizes_[0] == 1) {
        subEntity[0] = 0;
        return;
      }

      /* edge and vertex numbering
         0----0----1
       */
      unsigned lastIndex     = 0;
      subEntity[lastIndex++] = 0;  // corner 0
      for (unsigned i = 0; i < sizes_[0] - 2; ++i)
        subEntity[lastIndex++] = 0;  // inner dofs of element (0)

      subEntity[lastIndex++] = 1;  // corner 1

      assert(size() == lastIndex);
    }

    void setup2d(std::vector<unsigned int>& subEntity) {
      unsigned lastIndex = 0;

      // LocalKey: entity number , entity codim, dof indices within each entity
      /* edge and vertex numbering
         2----3----3
         |         |
         |         |
         0         1
         |         |
         |         |
         0----2----1
       */

      // lower edge (2)
      subEntity[lastIndex++] = 0;  // corner 0
      for (unsigned i = 0; i < sizes_[0] - 2; ++i)
        subEntity[lastIndex++] = 2;  // inner dofs of lower edge (2)

      subEntity[lastIndex++] = 1;  // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < sizes_[1] - 2; ++e) {
        subEntity[lastIndex++] = 0;  // left edge (0)
        for (unsigned i = 0; i < sizes_[0] - 2; ++i)
          subEntity[lastIndex++] = 0;  // face dofs
        subEntity[lastIndex++] = 1;    // right edge (1)
      }

      // upper edge (3)
      subEntity[lastIndex++] = 2;  // corner 2
      for (unsigned i = 0; i < sizes_[0] - 2; ++i)
        subEntity[lastIndex++] = 3;  // inner dofs of upper edge (3)

      subEntity[lastIndex++] = 3;  // corner 3

      assert(size() == lastIndex);
    }

  public:
    void init(const std::array<unsigned, dim>& sizes) {
      sizes_ = sizes;

      li_.resize(size());

      // Set up array of codimension-per-dof-number
      std::vector<unsigned int> codim(li_.size());

      for (std::size_t i = 0; i < codim.size(); i++) {
        codim[i] = 0;
        // Codimension gets increased by 1 for each coordinate direction
        // where dof is on boundary
        std::array<unsigned int, dim> mIdx = multiindex(i);
        for (int j = 0; j < dim; j++)
          if (mIdx[j] == 0 or mIdx[j] == sizes[j] - 1) codim[i]++;
      }

      // Set up index vector (the index of the dof in the set of dofs of a given subentity)
      // Algorithm: the 'index' has the same ordering as the dof number 'i'.
      // To make it consecutive we interpret 'i' in the (k+1)-adic system, omit all digits
      // that correspond to axes where the dof is on the element boundary, and transform the
      // rest to the (k-1)-adic system.
      std::vector<unsigned int> index(size());

      for (std::size_t i = 0; i < index.size(); i++) {
        index[i] = 0;

        std::array<unsigned int, dim> mIdx = multiindex(i);

        for (int j = dim - 1; j >= 0; j--)
          if (mIdx[j] > 0 and mIdx[j] < sizes[j] - 1) index[i] = (sizes[j] - 1) * index[i] + (mIdx[j] - 1);
      }

      // Set up entity and dof numbers for each (supported) dimension separately
      std::vector<unsigned int> subEntity(li_.size());

      if (!subEntity.empty()) {
        if (dim == 1) {
          setup1d(subEntity);

        } else if (dim == 2 and sizes_[0] > 1 and sizes_[1] > 1) {
          setup2d(subEntity);
        }
      }

      for (size_t i = 0; i < li_.size(); i++)
        li_[i] = LocalKey(subEntity[i], codim[i], index[i]);
    }

    //! number of coefficients
    [[nodiscard]] std::size_t size() const { return std::accumulate(sizes_.begin(), sizes_.end(), 1, std::multiplies<>()); }

    //! get i'th index
    [[nodiscard]] const LocalKey& localKey(std::size_t i) const { return li_[i]; }

  private:
    // Number of shape functions on this element per coordinate direction
    std::array<unsigned, dim> sizes_;

    std::vector<LocalKey> li_;
  };

  /** \brief Local interpolation in the sense of dune-localfunctions, for the B-spline basis on tensor-product grids
   *
   * \ingroup FunctionSpaceBasesImplementations
   */
  template <int dim, class LB>
  class NurbsLocalInterpolation {
  public:
    //! \brief Local interpolation of a function
    template <typename F, typename C>
    void interpolate(const F& f, std::vector<C>& out) const {
      DUNE_THROW(NotImplemented, "NurbsLocalInterpolation::interpolate");
    }
  };

  /** \brief LocalFiniteElement in the sense of dune-localfunctions, for the B-spline basis on tensor-product grids
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * This class ties together the implementation classes NurbsLocalBasis, NurbsLocalCoefficients, and NurbsLocalInterpolation
   *
   * \tparam D Number type used for domain coordinates
   * \tparam R Number type used for spline function values
   * \tparam MI Global multi-index type.  Only here for technical reasons
   */
  template <class GV, class R, class MI>
  class NurbsLocalFiniteElement {
    typedef typename GV::ctype D;
    enum { dim = GV::dimension };
    friend class NurbsLocalBasis<GV, R, MI>;

  public:
    /** \brief Export various types related to this LocalFiniteElement
     */
    typedef LocalFiniteElementTraits<NurbsLocalBasis<GV, R, MI>, NurbsLocalCoefficients<dim>,
                                     NurbsLocalInterpolation<dim, NurbsLocalBasis<GV, R, MI>>>
        Traits;

    /** \brief Constructor with a given B-spline basis
     */
    explicit NurbsLocalFiniteElement(const NurbsPreBasis<GV, MI>& preBasis) : preBasis_(preBasis), localBasis_(preBasis, *this) {}

    /** \brief Copy constructor
     */
    NurbsLocalFiniteElement(const NurbsLocalFiniteElement& other) : preBasis_(other.preBasis_), localBasis_(preBasis_, *this) {}

    /** \brief Bind LocalFiniteElement to a specific knot span of the spline patch
     *
     * Elements are the non-empty knot spans, here we do the renumbering
     *
     * \param ijk Integer coordinates in the tensor product patch
     */
    void bind(const std::array<unsigned, dim>& elementIdx) {
      const auto& patchData = preBasis_.patchData_;
      for (size_t i = 0; i < elementIdx.size(); i++) {
        currentKnotSpan_[i] = Dune::IGA::findSpan(patchData.degree[i], *(preBasis_.uniqueKnotVector_[i].begin() + elementIdx[i]),
                                                  patchData.knotSpans[i], elementIdx[i]);

        // Compute the geometric transformation from knotspan-local to global coordinates
        localBasis_.offset_[i] = preBasis_.patchData_.knotSpans[i][currentKnotSpan_[i]];
        localBasis_.scaling_[i][i]
            = preBasis_.patchData_.knotSpans[i][currentKnotSpan_[i] + 1] - preBasis_.patchData_.knotSpans[i][currentKnotSpan_[i]];
      }

      // Set up the LocalCoefficients object
      std::array<unsigned int, dim> sizes;
      for (size_t i = 0; i < dim; i++)
        sizes[i] = size(i);
      localCoefficients_.init(sizes);
    }

    /** \brief Hand out a LocalBasis object */
    const NurbsLocalBasis<GV, R, MI>& localBasis() const { return localBasis_; }

    /** \brief Hand out a LocalCoefficients object */
    const NurbsLocalCoefficients<dim>& localCoefficients() const { return localCoefficients_; }

    /** \brief Hand out a LocalInterpolation object */
    const NurbsLocalInterpolation<dim, NurbsLocalBasis<GV, R, MI>>& localInterpolation() const { return localInterpolation_; }

    /** \brief Number of shape functions in this finite element */
    [[nodiscard]] unsigned size() const {
      std::size_t r = 1;
      for (int i = 0; i < dim; i++)
        r *= size(i);
      return r;
    }

    /** \brief Return the reference element that the local finite element is defined on (here, a hypercube)
     */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(dim); }

    /** \brief Number of degrees of freedom for one coordinate direction */
    [[nodiscard]] unsigned int size(int i) const { return preBasis_.patchData_.degree[i] + 1; }

    const NurbsPreBasis<GV, MI>& preBasis_;

    NurbsLocalBasis<GV, R, MI> localBasis_;
    NurbsLocalCoefficients<dim> localCoefficients_;
    NurbsLocalInterpolation<dim, NurbsLocalBasis<GV, R, MI>> localInterpolation_;

    // The knot span we are bound to
    std::array<int, dim> currentKnotSpan_;
  };

  template <typename GV, typename MI>
  class NurbsNode;

  /** \brief Pre-basis for B-spline basis
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV The GridView that the space is defined on
   * \tparam MI Type to be used for multi-indices
   *
   * The BSplinePreBasis can be used to embed a BSplineBasis
   * in a larger basis for the construction of product spaces.
   */
  template <typename GV, class MI, Dune::IGA::LinearAlgebra NurbsGridLinearAlgebraTraits>
  class NurbsPreBasis {
    static const auto dim      = GV::dimension;
    static const auto dimworld = GV::dimensionworld;

    /** \brief Simple dim-dimensional multi-index class */
    class MultiDigitCounter {
    public:
      /** \brief Constructs a new multi-index, and sets all digits to zero
       *  \param limits Number of different digit values for each digit, i.e., digit i counts from 0 to limits[i]-1
       */
      explicit MultiDigitCounter(const std::array<unsigned int, dim>& limits) : limits_(limits) {
        std::fill(counter_.begin(), counter_.end(), 0);
      }

      /** \brief Increment the multi-index */
      MultiDigitCounter& operator++() {
        for (int i = 0; i < dim; i++) {
          ++counter_[i];

          // no overflow?
          if (counter_[i] < limits_[i]) break;

          counter_[i] = 0;
        }
        return *this;
      }

      /** \brief Access the i-th digit of the multi-index */
      const unsigned int& operator[](int i) const { return counter_[i]; }

      /** \brief How many times can you increment this multi-index before it overflows? */
      [[nodiscard]] unsigned int cycle() const {
        unsigned int r = 1;
        for (int i = 0; i < dim; i++)
          r *= limits_[i];
        return r;
      }

    private:
      /** \brief The number of different digit values for each place */
      const std::array<unsigned int, dim> limits_;

      /** \brief The values of the multi-index.  Each array entry is one digit */
      std::array<unsigned int, dim> counter_;
    };

  public:
    /** \brief The grid view that the FE space is defined on */
    using GridView  = GV;
    using size_type = std::size_t;

    using Node = NurbsNode<GV, MI>;

    //! Type of created tree node index set. \deprecated
    using IndexSet = Impl::DefaultNodeIndexSet<NurbsPreBasis>;

    /** \brief Type used for global numbering of the basis vectors */
    using MultiIndex = MI;

    using SizePrefix = Dune::ReservedVector<size_type, 1>;

    // Type used for function values
    using R = typename NurbsGridLinearAlgebraTraits::value_type;

    explicit NurbsPreBasis(const GridView& gridView) : NurbsPreBasis(gridView, gridView.getPatchData()) {}

    NurbsPreBasis(const GridView& gridView, const Dune::IGA::NURBSPatchData<dim, dimworld>& patchData)
        : gridView_{gridView}, patchData_{patchData} {
      for (int i = 0; i < dim; ++i)
        std::ranges::unique_copy(patchData_.knotSpans[i], std::back_inserter(uniqueKnotVector_[i]),
                                 [](auto& l, auto& r) { return Dune::FloatCmp::eq(l, r); });

      std::ranges::transform(uniqueKnotVector_, elements_.begin(), [](auto& v) { return v.size() - 1; });
    }

    //! Initialize the global indices
    void initializeIndices() {}

    //! Obtain the grid view that the basis is defined on
    const GridView& gridView() const { return gridView_; }

    //! Update the stored grid view, to be called if the grid has changed
    void update(const GridView& gv) { gridView_ = gv; }

    /**
     * \brief Create tree node
     */
    Node makeNode() const { return Node{this}; }

    //! Return number of possible values for next position in multi index
    [[nodiscard]] size_type size(const SizePrefix prefix) const {
      assert(prefix.empty() || prefix.size() == 1);
      return (prefix.empty()) ? size() : 0;
    }

    //! Get the total dimension of the space spanned by this basis
    [[nodiscard]] size_type dimension() const { return size(); }

    //! Get the maximal number of DOFs associated to node for any element
    [[nodiscard]] size_type maxNodeSize() const {
      size_type result = 1;
      for (int i = 0; i < dim; i++)
        result *= patchData_.degree[i] + 1;
      return result;
    }

    //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
    template <typename It>
    It indices(const Node& node, It it) const {
      // Local degrees of freedom are arranged in a lattice.
      // We need the lattice dimensions to be able to compute lattice coordinates from a local index
      std::array<unsigned int, dim> localSizes;
      for (int i = 0; i < dim; i++)
        localSizes[i] = node.finiteElement().size(i);
      for (size_type i = 0, end = node.size(); i < end; ++i, ++it) {
        std::array<unsigned int, dim> localIJK = getIJK(i, localSizes);

        const auto currentKnotSpan = node.finiteElement().currentKnotSpan_;
        const auto order           = patchData_.degree;

        std::array<unsigned int, dim> globalIJK;
        for (int j = 0; j < dim; j++)
          globalIJK[j] = std::max((int)currentKnotSpan[j] - (int)order[j], 0) + localIJK[j];  // needs to be a signed type!

        // Make one global flat index from the globalIJK tuple
        size_type globalIdx = globalIJK[dim - 1];

        for (int j = dim - 2; j >= 0; j--)
          globalIdx = globalIdx * sizePerDirection(j) + globalIJK[j];

        *it = {{globalIdx}};
      }
      return it;
    }

    //! \brief Total number of B-spline basis functions
    [[nodiscard]] unsigned int size() const {
      unsigned int result = 1;
      for (size_t i = 0; i < dim; i++)
        result *= sizePerDirection(i);
      return result;
    }

    //! \brief Number of shape functions in one direction
    [[nodiscard]] unsigned int sizePerDirection(size_t d) const { return patchData_.knotSpans[d].size() - patchData_.degree[d] - 1; }

    /** \brief Evaluate all B-spline basis functions at a given point
     */
    void evaluateFunction(const FieldVector<typename GV::ctype, dim>& in, std::vector<FieldVector<R, 1>>& out,
                          const std::array<int, dim>& currentKnotSpan) const {
      std::array<typename GV::ctype, dim> inArray;
      std::ranges::copy(in, inArray.begin());
      const auto N = IGA::Nurbs<dim, NurbsGridLinearAlgebraTraits>::basisFunctions(inArray, patchData_.knotSpans, patchData_.degree,
                                                                             extractWeights(patchData_.controlPoints), currentKnotSpan);
      out.resize(N.directSize());
      std::ranges::copy(N.directGetAll(), out.begin());
    }

    /** \brief Evaluate Jacobian of all B-spline basis functions
     *
     * In theory this is easy: just look up the formula in a B-spline text of your choice.
     * The challenge is compute only the values needed for the current knot span.
     */
    void evaluateJacobian(const FieldVector<typename GV::ctype, dim>& in, std::vector<FieldMatrix<R, 1, dim>>& out,
                          const std::array<int, dim>& currentKnotSpan) const {
      std::array<typename GV::ctype, dim> inArray;
      std::ranges::copy(in, inArray.begin());
      const auto dN = IGA::Nurbs<dim, NurbsGridLinearAlgebraTraits>::basisFunctionDerivatives(
          inArray, patchData_.knotSpans, patchData_.degree, extractWeights(patchData_.controlPoints), 1, false, currentKnotSpan);
      out.resize(dN.get(std::array<int, dim>{}).directSize());
      for (int j = 0; j < dim; ++j) {
        std::array<int, dim> multiIndex{};
        multiIndex[j]  = 1;
        const auto& dNcurrent = dN.get(multiIndex);
        for (int i = 0; i < dNcurrent.directSize(); ++i)
          out[i][0][j] = dNcurrent.directGet(i);
      }
    }

    //! \brief Evaluate Derivatives of all B-spline basis functions

    void partial(const std::array<int, dim>& order, const FieldVector<typename GV::ctype, dim>& in, std::vector<FieldVector<R, 1>>& out,
                 const std::array<unsigned, dim>& currentKnotSpan) const {
      std::array<typename GV::ctype, dim> inArray;
      std::ranges::copy(in, inArray.begin());

      const auto dN = IGA::Nurbs<dim, NurbsGridLinearAlgebraTraits>::basisFunctionDerivatives(inArray, patchData_.knotSpans, patchData_.degree,
                                                                                        extractWeights(patchData_.controlPoints), order);
      out     = dN.get(order).directGetAll();
    }

    /** \brief Compute integer element coordinates from the element index
     * \warning This method makes strong assumptions about the grid, namely that it is
     *   structured, and that indices are given in a x-fastest fashion.
     */
    static std::array<unsigned int, dim> getIJK(typename GridView::IndexSet::IndexType idx, std::array<unsigned int, dim> elements) {
      std::array<unsigned, dim> result;
      for (int i = 0; i < dim; i++) {
        result[i] = idx % elements[i];
        idx /= elements[i];
      }
      return result;
    }

    std::array<std::vector<double>, dim> uniqueKnotVector_;

    /** \brief Order of the B-spline for each space dimension */
    Dune::IGA::NURBSPatchData<dim, dimworld> patchData_;

    /** \brief Number of grid elements in the different coordinate directions */
    std::array<unsigned, dim> elements_;

    GridView gridView_;
  };

  template <typename GV, typename MI>
  class NurbsNode : public LeafBasisNode {
    static const int dim = GV::dimension;

  public:
    using size_type     = std::size_t;
    using Element       = typename GV::template Codim<0>::Entity;
    using FiniteElement = NurbsLocalFiniteElement<GV, double, MI>;

    NurbsNode(const NurbsPreBasis<GV, MI>* preBasis) : preBasis_(preBasis), finiteElement_(*preBasis) {}

    //! Return current element, throw if unbound
    const Element& element() const { return element_; }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
     */
    const FiniteElement& finiteElement() const { return finiteElement_; }

    //! Bind to element.
    void bind(const Element& e) {
      element_          = e;
      auto elementIndex = preBasis_->gridView().indexSet().index(e);
      finiteElement_.bind(preBasis_->getIJK(elementIndex, preBasis_->elements_));
      this->setSize(finiteElement_.size());
    }

  protected:
    const NurbsPreBasis<GV, MI>* preBasis_;

    FiniteElement finiteElement_;
    Element element_;
  };

  namespace BasisFactory {

    namespace Imp {

      template <std::integral auto dim, std::integral auto dimworld>
      class NurbsPreBasisFactory {
      public:
        static const std::size_t requiredMultiIndexSize = 1;

        explicit NurbsPreBasisFactory(const Dune::IGA::NURBSPatchData<dim, dimworld>& patchData) : patchData_(patchData) {}

        template <class MultiIndex, class GridView>
        auto makePreBasis(const GridView& gridView) const {
          return NurbsPreBasis<GridView, MultiIndex>(gridView, patchData_);
        }

      private:
        const Dune::IGA::NURBSPatchData<dim, dimworld>& patchData_;
      };

    }  // namespace Imp

    /**
     * \brief Create a pre-basis factory that can create a B-spline pre-basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     */
    template <std::integral auto dim, std::integral auto dimworld>
    auto nurbs(const Dune::IGA::NURBSPatchData<dim, dimworld>& data) {
      return Imp::NurbsPreBasisFactory<dim, dimworld>(data);
    }

  }  // end namespace BasisFactory

  // *****************************************************************************
  // This is the actual global basis implementation based on the reusable parts.
  // *****************************************************************************

  /** \brief A global B-spline basis
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV The GridView that the space is defined on
   */
  template <typename GV>
  using NurbsBasis = DefaultGlobalBasis<NurbsPreBasis<GV, FlatMultiIndex<std::size_t>>>;

}  // namespace Dune::Functions

#endif  // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NURBSBASIS_HH

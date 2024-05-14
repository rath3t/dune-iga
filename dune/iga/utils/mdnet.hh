// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <concepts>
#include <functional>
#include <iterator>
#include <numeric>
#include <ranges>

#include <dune/common/fvector.hh>
#include <dune/iga/hierarchicpatch/concepts.hh>

namespace Dune::IGA {

namespace Impl {
  template <std::integral auto netdim, typename ValueType, bool isConst = false>
  class HyperSurfaceIterator;

} // namespace Impl

/**
 * @brief This is a condensed version of the proposal
 * https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2022/p1684r2.html It is a container with multidimensional
 * accessors. \details It has less features and only allows runtime size values but the rank is given at compile time.
 * Furthermore, there are specific constructors to have a convenient construction for NURBS controlpoints
 * It also has several overloaded operators to do linear algebra with it such as operator*
 * @tparam ElementType The type of the stored values
 * @tparam netdim The rank of the array (How many dimensions it has)
 * @tparam Container The type of storage
 */
template <int netdim, class ElementType, class Container = std::vector<ElementType>>
class MultiDimensionalNet
{
public:
  static constexpr std::size_t netDim = netdim;

  constexpr MultiDimensionalNet()                                      = default;
  constexpr MultiDimensionalNet(const MultiDimensionalNet&)            = default;
  constexpr MultiDimensionalNet(MultiDimensionalNet&&)                 = default;
  constexpr MultiDimensionalNet& operator=(const MultiDimensionalNet&) = default;
  constexpr MultiDimensionalNet& operator=(MultiDimensionalNet&&)      = default;
  ~MultiDimensionalNet()                                               = default;

  using container_type  = Container;
  using element_type    = ElementType;
  using value_type      = std::remove_cv_t<element_type>;
  using index_type      = typename Container::size_type;
  using size_type       = index_type;
  using rank_type       = int;
  using pointer         = typename container_type::pointer;
  using reference       = typename container_type::reference;
  using const_pointer   = typename container_type::const_pointer;
  using const_reference = typename container_type::const_reference;

  MultiDimensionalNet(std::initializer_list<std::initializer_list<element_type>> values)
  requires(netdim == 2)
  {
    std::vector<std::vector<element_type>> vals;
    for (auto&& val : values)
      vals.push_back(val);
    std::array<int, netdim> dimsize = {static_cast<int>(vals.size()), static_cast<int>(values.begin()->size())};
    *this                           = MultiDimensionalNet{dimsize, vals};
  }

  explicit MultiDimensionalNet(const std::vector<std::vector<element_type>>& vals)
  requires(netdim == 2)
  {
    std::array<int, netdim> dimsize = {static_cast<int>(vals.size()), static_cast<int>(vals.begin()->size())};
    *this                           = MultiDimensionalNet{dimsize, vals};
  }

  explicit MultiDimensionalNet(const std::vector<element_type>& vals)
  requires(netdim == 1)
  {
    std::array<int, 1> dimsize = {static_cast<int>(vals.size())};
    *this                      = MultiDimensionalNet{dimsize, vals};
  }

  /** @brief constructor for a net of a certain strideSizes with values unknown.
   *
   *  @param[in] dimSize array of the strideSizes of each dimension
   */
  explicit MultiDimensionalNet(const std::array<int, netdim>& dimSize)
      : dimSize_(dimSize) {
    int size_ = 1;
    for (auto ds : dimSize)
      size_ *= ds;

    values_.resize(size_);
  }

  template <typename... Args>
  explicit MultiDimensionalNet(int dimSize0, Args&&... dimSize)
      : dimSize_({dimSize0, std::forward<Args>(dimSize)...}) {
    int size = 1;
    for (auto ds : dimSize_)
      size *= ds;

    values_.resize(size);
  }

  explicit MultiDimensionalNet(const FieldVector<int, netdim>& dimSize) {
    std::ranges::copy(dimSize, dimSize_.begin());
    int size = 1;
    for (auto ds : dimSize_)
      size *= ds;

    values_.resize(size);
  }

  /** @brief constructor for a net from outerproduct of vectors
   *
   *  @param[in] values netdim vectors of values
   */
  template <std::ranges::sized_range V>
  explicit MultiDimensionalNet(const std::array<V, netdim>& values) {
    for (int i = 0; i < netdim; ++i)
      dimSize_[i] = values[i].size();

    const int localSize = std::accumulate(dimSize_.begin(), dimSize_.end(), 1, std::multiplies{});
    values_.resize(localSize);
    for (int i = 0; i < localSize; i += dimSize_[0])
      std::copy(values[0].begin(), values[0].end(), values_.begin() + i);
    // https://godbolt.org/z/TbTPszGqT
    for (int curSize = dimSize_[0], k = 1; k < netdim; curSize *= (dimSize_[k]), ++k)
      for (int ij = 0; ij < localSize; ij += curSize * (dimSize_[k]))
        for (int j = 0; auto& v2 : values[k]) {
          auto countedView = std::views::counted(values_.begin() + ij + (j++) * curSize, curSize);
          std::ranges::transform(countedView, countedView.begin(), [&v2](auto& v) { return v * v2; });
        }
  }

  /** @brief constructor intended for the 1-D, or if the values are already in a vector
   *  \note can also be used if the values are already mapped
   *
   *  @param[in] dimSize array of the strideSizes of each dimension
   *  @param[in] values vector with values
   */
  MultiDimensionalNet(std::array<int, netdim> dimSize, const container_type& values)
      : values_(values),
        dimSize_{dimSize} {}

  MultiDimensionalNet(std::array<int, netdim> dimSize, std::ranges::range auto values)
      : values_(values.begin(), values.end()),
        dimSize_{dimSize} {}

  /** @brief constructor intended for the 2-D if the values are already in a matrix
   *
   *  @param[in] dimSize array of the strideSizes of each dimension
   *  @param[in] values matrix with values
   */
  MultiDimensionalNet(std::array<int, netdim> dimSize, const std::vector<std::vector<value_type>>& values)
  requires(netdim == 2)
      : dimSize_(dimSize) {
    values_.resize(values.size() * values[0].size());
    for (int i = 0; i < values.size(); ++i)
      for (int j = 0; j < values[0].size(); ++j)
        this->set({i, j}, values[i][j]);
  }

  /** @brief constructor intended for the 3-D
   *
   *  @param[in] dimSize array of the strideSizes of each dimension
   *  @param[in] values matrix with values
   */
  MultiDimensionalNet(std::array<int, netdim> dimSize, const std::vector<std::vector<std::vector<value_type>>>& values)
  requires(netdim == 3)
      : dimSize_(dimSize) {
    values_.resize(values.size() * values[0].size() * values[0][0].size());

    for (int i = 0; i < values.size(); ++i)
      for (int j = 0; j < values[0].size(); ++j)
        for (int k = 0; k < values[0][0].size(); ++k)
          this->set({i, j, k}, values[i][j][k]);
  }

  /** @brief constructor for a grid of the same value
   *
   *  @param[in] dimSize array of the strideSizes of each dimension
   *  @param[in] value a common value to fill the grid
   */
  MultiDimensionalNet(std::array<int, netdim> dimSize, const value_type& value)
      : dimSize_(dimSize) {
    int size = 1;
    for (auto ds : dimSize_)
      size *= ds;

    values_.resize(size);
    std::fill(values_.begin(), values_.end(), value);
  }

  /** @brief sets a value at the multiindex */
  void set(std::array<int, netdim> multiIndex, const value_type& value) {
    int index      = this->index(multiIndex);
    values_[index] = value;
  }

  /** @brief sets a value at the multiindex */
  void directSet(int index, const value_type& value) {
    values_[index] = value;
  }

  void directSet(int index, value_type&& value) {
    values_[index] = std::move(value);
  }

  template <typename... Args>
  auto& operator()(const Args... args) {
    return get({args...});
  }

  template <typename... Args>
  const auto& operator()(const Args... args) const {
    return get({args...});
  }

  /** @brief returns the value at the multiindex */
  template <typename ArrayType = std::array<int, netdim>>
  value_type& get(const ArrayType& multiIndex) {
    return values_[this->index(multiIndex)];
  }

  /** @brief returns the value at the multiindex */
  template <typename ArrayType = std::array<int, netdim>>
  const value_type& get(const ArrayType& multiIndex) const {
    return values_[this->index(multiIndex)];
  }

  auto subNet(const std::array<int, netdim>& start, const std::array<int, netdim>& size) const {
    std::vector<value_type> subValues;
    subValues.reserve(std::accumulate(size.begin(), size.end(), 1, std::multiplies{}));
    if constexpr (netdim == 1)
      for (int i = 0; i < size[0]; ++i)
        subValues.push_back(get({start[0] + i}));
    else if constexpr (netdim == 2)
      for (int j = 0; j < size[1]; ++j)
        for (int i = 0; i < size[0]; ++i)
          subValues.push_back(get({start[0] + i, start[1] + j}));
    else if constexpr (netdim == 3)
      for (int k = 0; k < size[2]; ++k) // @todo generalize
        for (int j = 0; j < size[1]; ++j)
          for (int i = 0; i < size[0]; ++i)
            subValues.push_back(get({start[0] + i, start[1] + j, start[2] + k}));

    return MultiDimensionalNet<netdim, value_type>(size, subValues);
  }

  /** @brief returns a value at an index (unmapped)
   * \note only to be used when the mapping is known
   */
  value_type& directGet(const int index) {
    return values_[index];
  }

  const value_type& directGet(const int index) const {
    return values_[index];
  }

  auto& directGetAll() {
    return values_;
  }

  const auto& directGetAll() const {
    return values_;
  }

  /** @brief returns a multiindex for a scalar index */
  template <typename ReturnType = std::array<int, netdim>>
  ReturnType directToMultiIndex(const int index) const {
    return directToMultiIndex<ReturnType>(dimSize_, index);
  }

  template <typename ReturnType = std::array<int, netdim>>
  static ReturnType directToMultiIndex(const std::array<int, netdim>& dimSize, const int index) {
    ReturnType multiIndex;

    int help = index;
    for (int i = 0; i < netdim; ++i) {
      int temp      = help % (dimSize[i]);
      multiIndex[i] = temp;
      help -= temp;
      help = help / dimSize[i];
    }
    return multiIndex;
  }

  /** @brief returns an array with the strideSizes of each dimension */
  std::array<int, netdim> strideSizes() const {
    return dimSize_;
  }

  /** @brief returns an array with the strideSizes of each dimension */
  template <std::integral T>
  std::array<T, netdim> sizeAsT() const {
    std::array<T, netdim> sizeUI;
    for (int i = 0; i < netdim; ++i)
      sizeUI[i] = static_cast<T>(dimSize_[i]);
    return sizeUI;
  }

  [[nodiscard]] std::size_t size() const {
    return values_.size();
  }

  void resize(std::array<int, netdim> dimSize) {
    dimSize_  = dimSize;
    int sizeT = 1;
    for (auto ds : dimSize_)
      sizeT *= ds;

    values_.resize(sizeT);
  }

  template <typename rValueType>
  requires Concept::MultiplyAssignAble<value_type, rValueType>
  MultiDimensionalNet<netdim, value_type>& operator*=(const MultiDimensionalNet<netdim, rValueType>& rnet) {
    assert(this->strideSizes() == rnet.strideSizes() && "The net dimensions need to match in each direction!");
    std::ranges::transform(values_, rnet.directGetAll(), values_.begin(), std::multiplies{});
    return *this;
  }

  template <typename rValueType>
  requires Concept::DivideAssignAble<value_type, rValueType>
  MultiDimensionalNet<netdim, value_type>& operator/=(const rValueType& div) {
    std::ranges::transform(values_, values_.begin(), [&div](auto& val) { return val / div; });
    return *this;
  }

  template <typename rValueType>
  requires Concept::DivideAssignAble<value_type, rValueType>
  MultiDimensionalNet<netdim, value_type>& operator/=(const MultiDimensionalNet<netdim, rValueType>& rnet) {
    assert(this->strideSizes() == rnet.strideSizes() && "The net dimensions need to match in each direction!");
    std::ranges::transform(values_, rnet.directGetAll(), values_.begin(),
                           [](const auto& lval, const auto& rval) { return lval / rval; });
    return *this;
  }

  template <typename rValueType>
  requires Concept::MultiplyAssignAble<value_type, rValueType>
  MultiDimensionalNet<netdim, value_type>& operator*=(const rValueType& fac) {
    std::ranges::transform(values_, values_.begin(), [&fac](const auto& val) { return val * fac; });
    return *this;
  }

  template <typename rValueType>
  requires Concept::AddAble<value_type, rValueType>
  MultiDimensionalNet<netdim, value_type>& operator+=(const MultiDimensionalNet<netdim, rValueType>& rnet) {
    assert(this->strideSizes() == rnet.strideSizes() && "The net dimensions need to match in each direction!");
    std::ranges::transform(values_, rnet.directGetAll(), values_.begin(), std::plus{});
    return *this;
  }

  template <typename rValueType>
  requires Concept::SubstractAble<value_type, rValueType>
  MultiDimensionalNet<netdim, value_type>& operator-=(const MultiDimensionalNet<netdim, rValueType>& rnet) {
    assert(this->strideSizes() == rnet.strideSizes() && "The net dimensions need to match in each direction!");
    std::ranges::transform(values_, rnet.directGetAll(), values_.begin(), std::minus{});
    return *this;
  }

  template <typename ArrayType = std::array<int, netdim>>
  size_t index(const ArrayType& multiIndex) const {
    assert(!std::ranges::any_of(multiIndex, [](int i) { return i < 0; }) &&
           "The passed multiIndex has negative values");
    assert(!std::ranges::any_of(multiIndex, [id = 0, this](int i) mutable { return i > dimSize_[id++] - 1; }) &&
           "The passed multiIndex has too large values");

    size_t index = 0;
    for (int i = 0; i < netdim; ++i) {
      int help = 1;
      for (int j = i - 1; j > -1; --j)
        help *= dimSize_[j];

      index += help * static_cast<size_t>(multiIndex[i]);
    }
    return index;
  }

  template <typename ArrayType = std::array<int, netdim>>
  bool isValid(const ArrayType& multiIndex) const {
    for (int i = 0; i < netdim; ++i) {
      if (multiIndex[i] > dimSize_[i] - 1 || (multiIndex[i] < 0))
        return false;
    }
    return true;
  }

  Impl::HyperSurfaceIterator<netdim, value_type> hyperSurfBegin(
      const std::array<int, (std::size_t)(netdim - 1)>& direction) {
    return Impl::HyperSurfaceIterator<netdim, value_type, false>(*this, direction, 0);
  }

  auto hyperSurfFront(const std::array<int, (std::size_t)(netdim - 1)>& direction) {
    return *this->hyperSurfBegin(direction);
  }

  Impl::HyperSurfaceIterator<netdim, value_type> hyperSurfEnd(
      const std::array<int, (std::size_t)(netdim - 1)>& direction) {
    int directionEnd;
    if constexpr (netdim != 0) {
      for (int dirI = 0, i = 0; i < netdim; ++i) {
        if (dirI < direction.size() && i == direction.at(dirI++))
          continue;
        directionEnd = this->strideSizes()[i];
        break;
      }
    } else
      directionEnd = this->strideSizes()[0];

    return Impl::HyperSurfaceIterator<netdim, value_type, false>(*this, direction, directionEnd);
  }

  auto hyperSurfBack(const std::array<int, (std::size_t)(netdim - 1)>& direction) {
    int directionEnd;
    if constexpr (netdim != 0) {
      for (int dirI = 0, i = 0; i < netdim; ++i) {
        if (dirI < direction.size() && i == direction.at(dirI++))
          continue;
        directionEnd = this->strideSizes()[i] - 1;
        break;
      }
    } else
      directionEnd = this->strideSizes()[0] - 1;

    return *Impl::HyperSurfaceIterator<netdim, value_type, false>(*this, direction, directionEnd);
  }

  Impl::HyperSurfaceIterator<netdim, value_type, true> hyperSurfBegin(
      const std::array<int, (std::size_t)(netdim - 1)>& direction) const {
    return Impl::HyperSurfaceIterator<netdim, value_type, true>(*this, direction, 0);
  }

  auto hyperSurfFront(const std::array<int, (std::size_t)(netdim - 1)>& direction) const {
    return *this->hyperSurfBegin(direction);
  }
  Impl::HyperSurfaceIterator<netdim, value_type, true> hyperSurfEnd(
      const std::array<int, (std::size_t)(netdim - 1)>& direction) const {
    int directionEnd;
    if constexpr (netdim != 0) {
      for (int dirI = 0, i = 0; i < netdim; ++i) {
        if (dirI < direction.size() && i == direction.at(dirI++))
          continue;
        directionEnd = this->strideSizes()[i];
        break;
      }
    } else
      directionEnd = this->strideSizes()[0];

    return Impl::HyperSurfaceIterator<netdim, value_type, true>(*this, direction, directionEnd);
  }

  auto hyperSurfBack(const std::array<int, (std::size_t)(netdim - 1)>& direction) const {
    int directionBack;
    if constexpr (netdim != 0) {
      for (int dirI = 0, i = 0; i < netdim; ++i) {
        if (dirI < direction.size() && i == direction.at(dirI++))
          continue;
        directionBack = this->strideSizes()[i] - 1;
        break;
      }
    } else
      directionBack = this->strideSizes()[0] - 1;

    return *Impl::HyperSurfaceIterator<netdim, value_type, true>(*this, direction, directionBack);
  }

private:
  std::array<int, netdim> dimSize_;
  Container values_;
};

template <size_t netdim, typename value_type>
MultiDimensionalNet(std::array<int, netdim> dimSize, const std::vector<std::vector<std::vector<value_type>>>& values)
    -> MultiDimensionalNet<netdim, value_type>;

namespace Impl {
  template <std::integral auto netdim, typename ValueType, bool isConst>
  class HyperSurfaceIterator
  {
  public:
    using MDNetType = std::conditional_t<isConst, const MultiDimensionalNet<netdim, ValueType>,
                                         MultiDimensionalNet<netdim, ValueType>>;
    HyperSurfaceIterator(MDNetType& net, const std::array<int, netdim - 1>& direction, int at)
        : net_{&net},
          direction_{direction},
          at_{at} {
      std::array<int, netdim - 1> indicesSurface;
      for (int i = 0; i < indicesSurface.size(); ++i)
        indicesSurface[i] = net.strideSizes().at(direction[i]);

      const int cpSize = std::accumulate(indicesSurface.begin(), indicesSurface.end(), 1, std::multiplies{});
      std::function<std::array<int, netdim - 1>(int)> func = [indicesSurface](auto i) {
        return MultiDimensionalNet<netdim - 1, int>::template directToMultiIndex<std::array<int, netdim - 1>>(
            indicesSurface, i);
      };
      viewOverIndices_ = std::ranges::iota_view(0, cpSize) | std::views::transform(func);
    }

    std::ranges::transform_view<std::ranges::iota_view<int, int>, std::function<std::array<int, netdim - 1>(int)>>
        viewOverIndices_;

    HyperSurfaceIterator& operator++() {
      ++at_;
      return *this;
    }

    HyperSurfaceIterator& operator--() {
      --at_;
      return *this;
    }

    HyperSurfaceIterator operator[](int index) {
      return HyperSurfaceIterator<netdim, ValueType>(*(this->net_), this->direction_, index);
    }

    HyperSurfaceIterator& operator-=(int index) {
      at_ -= index;
      return *this;
    }

    bool operator==(const HyperSurfaceIterator& r) const {
      assert(this->net_->strideSizes() == r.net_->strideSizes());
      return this->at_ == r.at_;
    }

    bool operator!=(const HyperSurfaceIterator& r) const {
      assert(this->net_->strideSizes() == r.net_->strideSizes());
      return this->at_ != r.at_;
    }

    auto operator*() {
      std::array<int, netdim> multiIndex;
      if constexpr (netdim != 1)
        for (int dirI = 0, i = 0; i < netdim; ++i) {
          if (dirI < direction_.size() && i == direction_.at(dirI++))
            continue;
          multiIndex[i] = at_;
          break;
        }
      else
        multiIndex[0] = at_;

      auto objectExtractor = [multiIndex, *this](auto mI) mutable -> auto& {
        if constexpr (netdim != 1)
          for (int i = 0; i < direction_.size(); ++i)
            multiIndex[direction_[i]] = mI[i];

        return net_->get(multiIndex);
      };
      return std::ranges::transform_view(viewOverIndices_, objectExtractor);
    }

    auto operator->() {
      return *this;
    }

    friend HyperSurfaceIterator operator+(const HyperSurfaceIterator& l, int inc) {
      return HyperSurfaceIterator(*(l.net_), l.direction_, l.at_ + inc);
    }

    friend HyperSurfaceIterator operator-(const HyperSurfaceIterator& l, int inc) {
      return l + (-inc);
    }

  private:
    MDNetType* net_;
    std::array<int, netdim - 1> direction_;
    int at_;
  };

} // namespace Impl

template <std::integral auto netdim, typename lValueType, typename rValueType>
requires Concept::MultiplyAble<lValueType, rValueType>
auto dot(const MultiDimensionalNet<netdim, lValueType>& lnet, const MultiDimensionalNet<netdim, rValueType>& rnet) {
  using ResultType = decltype(lnet.directGetAll()[0] * rnet.directGetAll()[0]);
  assert(lnet.strideSizes() == rnet.strideSizes() && "The net dimensions need to match in each direction!");
  return std::inner_product(lnet.directGetAll().begin(), lnet.directGetAll().end(), rnet.directGetAll().begin(),
                            ResultType(0.0));
}

template <std::integral auto netdim, typename lValueType, typename rValueType>
requires Concept::MultiplyAble<lValueType, rValueType>
auto operator-(const MultiDimensionalNet<netdim, lValueType>& lnet,
               const MultiDimensionalNet<netdim, rValueType>& rnet) {
  assert(lnet.strideSizes() == rnet.strideSizes() && "The net dimensions need to match in each direction!");
  MultiDimensionalNet<netdim, lValueType> res = lnet;
  std::ranges::transform(res.directGetAll(), rnet.directGetAll(), res.directGetAll().begin(), std::minus{});
  return res;
}

template <std::integral auto netdim, typename lValueType, typename rValueType>
requires Concept::DivideAble<lValueType, rValueType>
MultiDimensionalNet<netdim, lValueType> operator/(const MultiDimensionalNet<netdim, lValueType>& lnet,
                                                  const rValueType& div) {
  MultiDimensionalNet<netdim, lValueType> res = lnet;
  std::ranges::transform(res.directGetAll(), res.directGetAll().begin(), [&div](auto& val) { return val / div; });
  return res;
}

template <std::integral auto netdim, typename lValueType, typename rValueType>
requires Concept::DivideAble<lValueType, rValueType>
MultiDimensionalNet<netdim, lValueType> operator*(const MultiDimensionalNet<netdim, lValueType>& lnet,
                                                  const rValueType& fac) {
  MultiDimensionalNet<netdim, lValueType> res = lnet;
  std::ranges::transform(res.directGetAll(), res.directGetAll().begin(), [&fac](auto& val) { return val * fac; });
  return res;
}

template <std::integral auto netdim, typename lValueType, typename rValueType>
requires Concept::DivideAble<lValueType, rValueType>
MultiDimensionalNet<netdim, lValueType> operator*(const rValueType& fac,
                                                  const MultiDimensionalNet<netdim, lValueType>& lnet) {
  return lnet * fac;
}

/** @brief class holds a n-dim net */
template <std::size_t netdim>
class MultiDimensionalNetIndex // FIXME merge with Net
{
public:
  explicit MultiDimensionalNetIndex(const FieldVector<int, netdim>& dimSize) {
    std::ranges::copy(dimSize, dimSize_.begin());
    size_ = 1;
    for (auto ds : dimSize_)
      size_ *= ds;
  }
  /** @brief returns a multiindex for a scalar index */
  template <typename ReturnType = std::array<int, netdim>>
  ReturnType directToMultiIndex(const int index) const {
    return directToMultiIndex<ReturnType>(dimSize_, index);
  }

  template <typename ReturnType = std::array<int, netdim>>
  static ReturnType directToMultiIndex(const std::array<int, netdim>& dimSize, const int index) {
    ReturnType multiIndex;

    int help = index;
    int temp;
    for (int i = 0; i < netdim; ++i) {
      temp          = help % (dimSize[i]);
      multiIndex[i] = temp;
      help -= temp;
      help = help / dimSize[i];
    }
    return multiIndex;
  }

  template <typename ArrayType = std::array<int, netdim>>
  int index(const ArrayType& multiIndex) const {
    int index{}, help;
    assert(!std::ranges::any_of(multiIndex, [](int i) { return i < 0; }) &&
           "The passed multiIndex has negative values");
    assert(!std::ranges::any_of(multiIndex, [id = 0, this](int i) mutable { return i > dimSize_[id++] - 1; }) &&
           "The passed multiIndex has too large values");

    for (int i = 0; i < netdim; ++i) {
      help = 1;
      for (int j = i - 1; j > -1; --j)
        help *= dimSize_[j];

      index += help * multiIndex[i];
    }
    return index;
  }

  int size() const {
    return size_;
  }

private:
  std::array<int, netdim> dimSize_;
  int size_;
};

} // namespace Dune::IGA

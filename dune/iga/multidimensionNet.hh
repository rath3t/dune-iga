//
// Created by lex on 21.10.21.
//

#pragma once

#include <ranges>

namespace Dune::IGA {
  /** \brief class holds a n-dim net */
  template <int netdim, typename ValueType>
  class MultiDimensionNet {
    using value_type = ValueType;

  public:
    /** \brief constructor for a net of a certain size with values unknown.
     *
     *  \param[in] dimSize array of the size of each dimension
     */
    explicit MultiDimensionNet(std::array<unsigned int, netdim> dimSize) : dimSize_(dimSize) {
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];

      values_.resize(size);
    }

    /** \brief constructor intended for the 1-D if the values are already in a vector
     *  \note can also be used if the values are already mapped
     *
     *  \param[in] dimSize array of the size of each dimension
     *  \param[in] values vector with values
     */
    MultiDimensionNet(std::array<unsigned int, netdim> dimSize, const std::vector<ValueType> values)
        : values_(values), dimSize_{dimSize} {}

    /** \brief constructor intended for the 2-D if the values are already in a matrix
     *
     *  \param[in] dimSize array of the size of each dimension
     *  \param[in] values matrix with values
     */
    MultiDimensionNet(std::array<unsigned int, netdim> dimSize, const std::vector<std::vector<ValueType>> values)
        : dimSize_(dimSize) {
      values_.resize(values.size() * values[0].size());

      for (unsigned int i = 0; i < values.size(); ++i) {
        for (unsigned int j = 0; j < values[0].size(); ++j) {
          std::array<unsigned int, netdim> multiIndex = {i, j};
          this->set(multiIndex, values[i][j]);
        }
      }
    }

    /** \brief constructor for a grid of the same value
     *
     *  \param[in] dimSize array of the size of each dimension
     *  \param[in] value a common value to fill the grid
     */
    MultiDimensionNet(std::array<unsigned int, netdim> dimSize, const ValueType& value) : dimSize_(dimSize) {
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];
      values_.resize(size);
      std::fill(values_.begin(), values_.end(), value);
    }

    /** \brief sets a value at the multiindex */
    void set(std::array<unsigned int, netdim> multiIndex, const ValueType& value) {
      int index      = this->index(multiIndex);
      values_[index] = value;
    }

    /** \brief sets a value at the multiindex */
    void directSet(unsigned int index, ValueType value) { values_[index] = value; }

    template <typename... Args>
    auto& operator()(const Args... args) {
      return get({args...});
    }

    template <typename... Args>
    const auto& operator()(const Args... args) const {
      return get({args...});
    }

    /** \brief returns the value at the multiindex */
    ValueType& get(const std::array<unsigned int, netdim>& multiIndex) {
      int index = this->index(multiIndex);
      return values_[index];
    }

    /** \brief returns the value at the multiindex */
    const ValueType& get(const std::array<unsigned int, netdim>& multiIndex) const {
      int index = this->index(multiIndex);
      return values_[index];
    }

    /** \brief returns a value at an index (unmapped)
     * \note only to be used when the mapping is known
     */
    ValueType& directGet(const int index)  { return values_[index]; }

    const ValueType& directGet(const int index) const { return values_[index]; }

    auto& directGetAll() { return values_; }

    /** \brief returns a multiindex for a scalar index */
    std::array<unsigned int, netdim> directToMultiIndex(const unsigned int index) const {
      std::array<unsigned int, netdim> multiIndex;

      unsigned int help = index;
      int temp;
      for (int i = 0; i < netdim; ++i) {
        temp                         = help % (dimSize_[netdim - (i + 1)]);
        multiIndex[netdim - (i + 1)] = temp;
        help -= temp;
        help = help / dimSize_[netdim - (i + 1)];
      }

      return multiIndex;
    }

    /** \brief returns an array with the size of each dimension */
    std::array<unsigned int, netdim> size() const { return dimSize_; }

    [[nodiscard]] unsigned int directSize() const { return values_.size(); }

    void resize(std::array<unsigned int, netdim> dimSize) {
      dimSize_ = dimSize;
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];

      values_.resize(size);
    }

  private:
    int index(const std::array<unsigned int, netdim>& multiIndex) const {
      int index, help;
      index = 0;
      for (int i = 0; i < netdim; ++i) {
        help = 1;
        for (int j = (i + 1); j < netdim; ++j)
          help *= dimSize_[j];

        index += help * multiIndex[i];
      }
      return index;
    }

  private:
    std::array<unsigned int, netdim> dimSize_;
    std::vector<ValueType> values_;
  };

//  template <typename... Args>
//  struct At {
//    std::array< int, sizeof...(Args)> args;
//  };
//
//  template <typename... Args>
//  auto at(Args&&... args) {
//    return At<Args&&...>{std::forward<Args>(args)...};
//  }

  template <int netdim, typename ValueType>
  auto line(MultiDimensionNet<netdim, ValueType>& net, const int direction, const int at) {

      std::array<unsigned int, netdim> multiIndex;
        for ( int argCounter = 0, i = 0; i < netdim; ++i) {
          if (i == direction && netdim>1) continue;
          multiIndex[i] = at;
        }
      int indicesEnd = (netdim==1 ) ? 1 : static_cast<int>(net.size()[direction]);
      auto indices = std::ranges::iota_view{0, indicesEnd};

      auto objectExtractor = [ multiIndex, &net, direction ](auto i) mutable -> auto& {
        if constexpr (netdim!=1)
          multiIndex[direction] = i;
        return net.get(multiIndex);
      };
      return std::ranges::transform_view(indices, objectExtractor);
  }

  template <int netdim, typename ValueType>
  auto line(MultiDimensionNet<netdim, ValueType> const& net, const int direction, const int at) {
      std::array<unsigned int, netdim> multiIndex;

        for ( int argCounter = 0, i = 0; i < netdim-1; ++i) {
          if (i == direction && netdim > 1) continue;
          multiIndex[i] = at;
        }

      int indicesEnd = (netdim==1  ) ? 1 : static_cast<int>(net.size()[direction]);
      auto indices = std::ranges::iota_view{0, indicesEnd};

      auto objectExtractor = [ multiIndex, &net, direction ](const auto i) mutable -> const auto& {
        if constexpr (netdim!=1)
          multiIndex[direction] = i;
        return net.get(multiIndex);
      };
      return std::ranges::transform_view(indices, objectExtractor);
  }

  template <typename ValueType>
  auto line(MultiDimensionNet<1, ValueType>& net) {
    return line(net, 0, {});
  }

  template <typename ValueType>
  auto line(MultiDimensionNet<1, ValueType> const& net) {
    return line(net, 0, {});
  }

}  // namespace Dune::IGA

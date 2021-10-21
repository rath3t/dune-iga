//
// Created by lex on 21.10.21.
//

#pragma once

namespace Dune::IGA {
  /** \brief class holds a n-dim net */
  template <int netdim, int dimworld, typename ValueType = FieldVector<double, dimworld>>
  class MultiDimensionNet {
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
    MultiDimensionNet(std::array<unsigned int, netdim> dimSize, const ValueType &value) : dimSize_(dimSize) {
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];
      values_.resize(size);
      std::fill(values_.begin(), values_.end(), value);
    }

    /** \brief sets a value at the multiindex */
    void set(std::array<unsigned int, netdim> multiIndex, const ValueType &value) {
      int index      = this->index(multiIndex);
      values_[index] = value;
    }

    /** \brief sets a value at the multiindex */
    void directSet(unsigned int index, FieldVector<double, dimworld> value) { values_[index] = value; }

    /** \brief returns the value at the multiindex */
    FieldVector<double, dimworld> get(std::array<unsigned int, netdim> multiIndex) const {
      int index = this->index(multiIndex);
      return values_[index];
    }

    /** \brief returns a value at an index (unmapped)
     * \note only to be used when the mapping is known
     */
    FieldVector<double, dimworld> directGet(int index) const { return values_[index]; }

    /** \brief returns a multiindex for a scalar index */
    std::array<unsigned int, netdim> directToMultiIndex(unsigned int index) const {
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

  private:
    int index(std::array<unsigned int, netdim> multiIndex) const {
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
}
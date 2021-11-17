//
// Created by lex on 21.10.21.
//

#pragma once

#include <ranges>
#include <concepts>
#include <numeric>
#include <dune/iga/concepts.hh>

namespace Dune::IGA {
  /** \brief class holds a n-dim net */
  template <std::integral auto netdim, typename ValueType>
  class MultiDimensionNet {
    using value_type = ValueType;

  public:
    MultiDimensionNet() = default;


    MultiDimensionNet(std::initializer_list<std::initializer_list<ValueType>>  values)
{
      std::vector<std::vector<ValueType>> vals;
      for( auto& val : values)
        vals.push_back(val);
        std::array< int, netdim> dimsize
            = {static_cast< int>(vals.size()), static_cast< int>(values.begin()->size())};
      *this = MultiDimensionNet{dimsize,vals};
    }

    /** \brief constructor for a net of a certain size with values unknown.
     *
     *  \param[in] dimSize array of the size of each dimension
     */
    explicit MultiDimensionNet(const std::array< int, netdim>& dimSize) : dimSize_(dimSize) {
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];

      values_.resize(size);
    }

    template<typename ...Args>
    explicit MultiDimensionNet(int dimSize0, Args&&... dimSize)  : dimSize_({dimSize0,std::forward<Args>(dimSize)...}) {
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];

      values_.resize(size);
    }

    explicit MultiDimensionNet(const FieldVector< int, netdim>& dimSize) {
      std::ranges::copy(dimSize,dimSize_.begin());
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];

      values_.resize(size);
    }

    explicit MultiDimensionNet( FieldVector< int, netdim>&& dimSize)  {
      std::ranges::copy(dimSize,dimSize_.begin());
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];

      values_.resize(size);
    }

    /** \brief constructor for a net from outerproduct of vectors
     *
     *  \param[in] values netdim vectors of values
     */
    explicit MultiDimensionNet(const std::array<std::vector<ValueType>, netdim>& values)
    {
      for (int i = 0; i < netdim; ++i)
        dimSize_[i] =values[i].size();

      const int localSize = std::accumulate(dimSize_.begin(), dimSize_.end(), 1,std::multiplies{});
      values_.resize(localSize);
      for (int i = 0; i < localSize; i += dimSize_[0])
        std::copy(values[0].begin(), values[0].end(), values_.begin() + i);
      //https://godbolt.org/z/TbTPszGqT
      for (int curSize=dimSize_[0], k=1; k< netdim; curSize*=(dimSize_[k]),++k )
        for (int ij = 0; ij < localSize; ij += curSize * (dimSize_[k]))
          for (int j = 0; auto& v2 : values[k]) {
            auto countedView = std::views::counted(values_.begin() + ij + (j++) * curSize, curSize);
            std::ranges::transform(countedView, countedView.begin(), [&v2](auto& v) { return v * v2; });
          }
    }

    /** \brief constructor intended for the 1-D if the values are already in a vector
     *  \note can also be used if the values are already mapped
     *
     *  \param[in] dimSize array of the size of each dimension
     *  \param[in] values vector with values
     */
    MultiDimensionNet(std::array<int, netdim> dimSize, const std::vector<ValueType>& values)
        : values_(values), dimSize_{dimSize} {}


    MultiDimensionNet(std::array<int, netdim> dimSize, std::ranges::range auto values)
        : values_(values.begin(),values.end()), dimSize_{dimSize} {}

    /** \brief constructor intended for the 2-D if the values are already in a matrix
     *
     *  \param[in] dimSize array of the size of each dimension
     *  \param[in] values matrix with values
     */
    MultiDimensionNet(std::array<int, netdim> dimSize, const std::vector<std::vector<ValueType>> values)
        : dimSize_(dimSize) {
      values_.resize(values.size() * values[0].size());

      for ( int i = 0; i < values.size(); ++i) {
        for ( int j = 0; j < values[0].size(); ++j) {
          std::array< int, netdim> multiIndex = {i, j};
          this->set(multiIndex, values[i][j]);
        }
      }
    }

    /** \brief constructor for a grid of the same value
     *
     *  \param[in] dimSize array of the size of each dimension
     *  \param[in] value a common value to fill the grid
     */
    MultiDimensionNet(std::array<int, netdim> dimSize, const ValueType& value) : dimSize_(dimSize) {
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];
      values_.resize(size);
      std::fill(values_.begin(), values_.end(), value);
    }

    /** \brief sets a value at the multiindex */
    void set(std::array<int, netdim> multiIndex, const ValueType& value) {
      int index      = this->index(multiIndex);
      values_[index] = value;
    }

    /** \brief sets a value at the multiindex */
    void directSet(int index, const ValueType& value) { values_[index] = value; }

    void directSet(int index, ValueType&& value) { values_[index] = std::move(value); }

    template <typename... Args>
    auto& operator()(const Args... args) {
      return get({args...});
    }



    template <typename... Args>
    const auto& operator()(const Args... args) const {
      return get({args...});
    }

    /** \brief returns the value at the multiindex */
    template<typename ArrayType= std::array<int, netdim>>
    ValueType& get(const ArrayType& multiIndex) {
      int index = this->index(multiIndex);
      return values_[index];
    }

    /** \brief returns the value at the multiindex */
    template<typename ArrayType= std::array<int, netdim>>
    const ValueType& get(const ArrayType& multiIndex) const {
      int index = this->index(multiIndex);
      return values_[index];
    }

    auto subNet(const std::array<int, netdim>& start, const std::array<int, netdim>& size) const
    {
      std::vector<ValueType>  subValues;
      subValues.reserve(std::accumulate(size.begin(),size.end(),1,std::multiplies{}));
      if constexpr (netdim==1)
        for (int i = 0; i < size[0]; ++i)
          subValues.push_back(get({start[0]+i}));
      else if constexpr (netdim==2)
          for (int j = 0; j < size[1]; ++j)
        for (int i = 0; i < size[0]; ++i)
            subValues.push_back(get({start[0]+i,start[1]+j}));
      else if constexpr (netdim==3)
            for (int k = 0; k < size[2]; ++k) //TODO generalize
          for (int j = 0; j < size[1]; ++j)
        for (int i = 0; i < size[0]; ++i)
              subValues.push_back(get({start[0]+i,start[1]+j,start[2]+k}));

      return MultiDimensionNet<netdim,ValueType>(size,subValues);
    }


    /** \brief returns a value at an index (unmapped)
     * \note only to be used when the mapping is known
     */
    ValueType& directGet(const int index)  { return values_[index]; }

    const ValueType& directGet(const int index) const { return values_[index]; }

    auto& directGetAll() { return values_; }

    const auto& directGetAll() const { return values_; }

    /** \brief returns a multiindex for a scalar index */
    template<typename ReturnType=std::array< int, netdim>>
    ReturnType directToMultiIndex(const int index) const {
      ReturnType multiIndex;

      int help = index;
      int temp;
      for (int i = 0; i < netdim; ++i) {
        temp                         = help % (dimSize_[i]);
        multiIndex[i] = temp;
        help -= temp;
        help = help / dimSize_[i];
      }
      return multiIndex;
    }

    /** \brief returns an array with the size of each dimension */
    std::array< int, netdim> size() const { return dimSize_; }

    [[nodiscard]]  std::size_t directSize() const { return values_.size(); }

    void resize(std::array< int, netdim> dimSize) {
      dimSize_ = dimSize;
      int size = 1;
      for (int i = 0; i < netdim; ++i)
        size *= dimSize_[i];

      values_.resize(size);
    }


    template<typename rValueType> requires MultiplyAssignAble<ValueType,rValueType>
    MultiDimensionNet<netdim,ValueType>& operator*=(const MultiDimensionNet<netdim,rValueType>& rnet)
    {
      assert(this->size()== rnet.size() && "The net dimensions need to match in each direction!");
      std::ranges::transform(values_,rnet.directGetAll(),values_.begin(),std::multiplies{});
      return *this;
    }

    template<typename rValueType> requires DivideAble<ValueType,rValueType>
        MultiDimensionNet<netdim,ValueType>& operator/=(const rValueType& div)
    {
      std::ranges::transform(values_,values_.begin(),[&div](auto& val){return val/div;});
      return *this;
    }

    template<typename rValueType> requires DivideAble<ValueType,rValueType>
        MultiDimensionNet<netdim,ValueType>& operator/=(const MultiDimensionNet<netdim,rValueType>& rnet)
    {
      assert(this->size()== rnet.size() && "The net dimensions need to match in each direction!");
      std::ranges::transform(values_,rnet.directGetAll(),values_.begin(),[](auto& lval, auto& rval){return lval/rval;});
      return *this;
    }

    template<typename rValueType> requires MultiplyAssignAble<ValueType,rValueType>
        MultiDimensionNet<netdim,ValueType>& operator*=(const rValueType& fac)
    {
      std::ranges::transform(values_,values_.begin(),[&fac](auto& val){return val*fac;});
      return *this;
    }

    template<typename rValueType> requires AddAble<ValueType,rValueType>
        MultiDimensionNet<netdim,ValueType>& operator+=(const MultiDimensionNet<netdim,rValueType>& rnet)
    {
      assert(this->size()== rnet.size() && "The net dimensions need to match in each direction!");
      std::ranges::transform(values_,rnet.directGetAll(),values_.begin(),std::plus{});
      return *this;
    }

    template<typename rValueType> requires SubstractAble<ValueType,rValueType>
        MultiDimensionNet<netdim,ValueType>& operator-=(const MultiDimensionNet<netdim,rValueType>& rnet)
    {
      assert(this->size()== rnet.size() && "The net dimensions need to match in each direction!");
      std::ranges::transform(values_,rnet.directGetAll(),values_.begin(),std::minus{});
      return *this;
    }

    template<typename ArrayType=std::array< int, netdim>>
    int index(const ArrayType& multiIndex) const {
      int index, help;
      index = 0;
      for (int i = 0; i < netdim; ++i) {
        help = 1;
        for (int j = i-1; j > -1; --j)
          help *= dimSize_[j];

        index += help * multiIndex[i];
      }
      return index;
    }


  private:
    std::array<int, netdim> dimSize_;
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

  template <std::integral auto  netdim, typename ValueType>
  auto line(MultiDimensionNet<netdim, ValueType>& net, const int direction, const int at) {

      std::array<int, netdim> multiIndex;
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

  template <std::integral auto  netdim, typename ValueType>
  auto line(MultiDimensionNet<netdim, ValueType> const& net, const int direction, const int at) {
      std::array<int, netdim> multiIndex;

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

//  template<int netdim,typename lValueType,typename rValueType> requires MultiplyAble<lValueType,rValueType>
//  auto operator *(const MultiDimensionNet<netdim,lValueType>& lnet,const MultiDimensionNet<netdim,rValueType>& rnet)
//  {
//    MultiDimensionNet<netdim,lValueType> res (lnet.size());
//    std::ranges::transform(lnet.directGetAll(),rnet.directGetAll(),res.directGetAll().begin(),std::multiplies{});
//    return res;
//  }

  template<std::integral auto netdim,typename lValueType,typename rValueType> requires MultiplyAble<lValueType,rValueType>
  auto dot(const MultiDimensionNet<netdim,lValueType>& lnet,const MultiDimensionNet<netdim,rValueType>& rnet)
  {
//    using ResultType = decltype(operator*(std::declval<lValueType>(),std::declval<rValueType>()));
    using ResultType = decltype(lnet.directGetAll()[0]*rnet.directGetAll()[0]);
    assert(lnet.size()== rnet.size() && "The net dimensions need to match in each direction!");
    return std::inner_product(lnet.directGetAll().begin(),lnet.directGetAll().end(),rnet.directGetAll().begin(),ResultType(0.0));
  }

  template<std::integral auto  netdim,typename lValueType,typename rValueType> requires MultiplyAble<lValueType,rValueType>
  auto operator-(const MultiDimensionNet<netdim,lValueType>& lnet,const MultiDimensionNet<netdim,rValueType>& rnet)
  {
    assert(lnet.size()== rnet.size() && "The net dimensions need to match in each direction!");
    MultiDimensionNet<netdim,lValueType> res = lnet;
    std::ranges::transform(res.directGetAll(),rnet.directGetAll(),res.directGetAll().begin(),std::minus{});
    return res;
  }

  template<std::integral auto  netdim,typename lValueType,typename rValueType> requires DivideAble<lValueType,rValueType>
      MultiDimensionNet<netdim,lValueType> operator/(const MultiDimensionNet<netdim,lValueType>& lnet,const rValueType& div)
  {
    MultiDimensionNet<netdim,lValueType> res = lnet;
    std::ranges::transform(res.directGetAll(),res.directGetAll().begin(),[&div](auto& val){return val/div;});
    return res;
  }

  template<std::integral auto  netdim,typename lValueType,typename rValueType> requires DivideAble<lValueType,rValueType>
      MultiDimensionNet<netdim,lValueType> operator*(const MultiDimensionNet<netdim,lValueType>& lnet,const rValueType& fac)
  {
    MultiDimensionNet<netdim,lValueType> res = lnet;
    std::ranges::transform(res.directGetAll(),res.directGetAll().begin(),[&fac](auto& val){return val*fac;});
    return res;
  }




  template<std::integral auto  netdim,typename lValueType,typename rValueType> requires DivideAble<lValueType,rValueType>
      MultiDimensionNet<netdim,lValueType> operator*(const rValueType& fac,const MultiDimensionNet<netdim,lValueType>& lnet)
  {
    return lnet*fac;
  }


//  //! calculate the binomial coefficient n over k for netdim directions and multiplies them and store the in the net
//  template<std::size_t netdim>
//  auto createBinomialNet(const std::array<int, netdim>& n)
//  {
//    MultiDimensionNet<netdim,int> binomNet(n);
//
//    for (int i = 0; i < binomNet.directSize(); ++i) {
//      auto multiIndex = binomNet.directToMultiIndex(i);
//      int binomialfac = 1;
//      for (int j=0 ; j< netdim ; ++j)
//        binomialfac*= Dune::binomial(n[j],multiIndex[j]);
//      binomNet.directSet(i,binomialfac);
//
//    }
//    return binomNet;
//  }



}  // namespace Dune::IGA

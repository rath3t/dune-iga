#pragma once



#include <nlohmann/json.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>


namespace anurbs {

using Index = std::ptrdiff_t;

template <typename T>
inline Index length(const T& container) noexcept
{
    return static_cast<Index>(std::size(container));
}


const double Infinity = std::numeric_limits<double>::infinity();
const double QuietNaN = std::numeric_limits<double>::quiet_NaN();


template <typename T, typename... TArgs>
std::unique_ptr<T> new_(TArgs&&... args)
{
    return std::make_unique<T>(std::forward<TArgs>(args)...);
}

template <size_t I,typename T> 
struct tuple_n
{
    template <typename... Args>
    using type = typename tuple_n<I - 1, T>::template type<T, Args...>;
};

template <typename T> 
struct tuple_n<0, T> {
    template <typename... Args>
    using type = std::tuple<Args...>;   
};

template <size_t I, typename T>
using tuple_of = typename tuple_n<I, T>::template type<>;

} // namespace anurbs
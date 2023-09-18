#ifndef PYTHONIC_NUMPY_ALL_HPP
#define PYTHONIC_NUMPY_ALL_HPP

#include "pythonic/include/numpy/all.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/numpy/multiply.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  bool _all(E begin, E end, utils::int_<1>)
  {
    return std::all_of(begin, end,
                       [](typename std::iterator_traits<E>::value_type e)
                           -> bool { return e; });
  }

  template <class E, size_t N>
  bool _all(E begin, E end, utils::int_<N>)
  {
    for (; begin != end; ++begin)
      if (!_all((*begin).begin(), (*begin).end(), utils::int_<N - 1>()))
        return false;
    return true;
  }

  template <class E>
  typename std::enable_if<types::is_numexpr_arg<E>::value, bool>::type
  all(E const &expr, types::none_type)
  {
    return _all(expr.begin(), expr.end(), utils::int_<E::value>());
  }

  template <class E>
  typename std::enable_if<
      std::is_scalar<E>::value || types::is_complex<E>::value, bool>::type
  all(E const &expr, types::none_type)
  {
    return expr;
  }

  template <class E>
  auto all(E const &array, long axis) ->
      typename std::enable_if<std::is_scalar<E>::value ||
                                  types::is_complex<E>::value,
                              decltype(all(array))>::type
  {
    if (axis != 0)
      throw types::ValueError("axis out of bounds");
    return all(array);
  }

  template <class E>
  auto all(E const &array, long axis) ->
      typename std::enable_if<E::value == 1, decltype(all(array))>::type
  {
    if (axis != 0)
      throw types::ValueError("axis out of bounds");
    return all(array);
  }

  template <class E>
  typename std::enable_if<
      E::value != 1,
      types::ndarray<typename E::dtype, types::array<long, E::value - 1>>>::type
  all(E const &array, long axis)
  {
    constexpr long N = E::value;
    typedef typename E::dtype T;
    if (axis < 0 || axis >= long(N))
      throw types::ValueError("axis out of bounds");
    if (axis == 0) {
      types::array<long, N - 1> shp;
      sutils::copy_shape<0, 1>(shp, array, utils::make_index_sequence<N - 1>());
      types::ndarray<bool, types::array<long, N - 1>> out(shp, true);
      return std::accumulate(array.begin(), array.end(), out,
                             functor::multiply());
    } else {
      types::array<long, N - 1> shp;
      sutils::copy_shape<0, 0>(shp, array, utils::make_index_sequence<N - 1>());
      types::ndarray<bool, types::array<long, N - 1>> ally(shp, builtins::None);
      std::transform(
          array.begin(), array.end(), ally.begin(),
          [=](types::ndarray<T, types::array<long, N - 1>> const &other) {
            return all(other, axis - 1);
          });
      return ally;
    }
  }
}
PYTHONIC_NS_END

#endif

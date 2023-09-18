#ifndef PYTHONIC_INCLUDE_NUMPY_PARTIAL_SUM_HPP
#define PYTHONIC_INCLUDE_NUMPY_PARTIAL_SUM_HPP

#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class Op, class E>
  using result_dtype = types::dtype_t<decltype(std::declval<Op>()(
      std::declval<typename std::remove_reference<E>::type::dtype>(),
      std::declval<typename std::remove_reference<E>::type::dtype>()))>;

  template <class Op, class E, class dtype = result_dtype<Op, E>>
  types::ndarray<typename dtype::type, types::pshape<long>>
  partial_sum(E const &expr, dtype d = dtype());

  template <class Op, class E, class dtype = result_dtype<Op, E>>
  auto partial_sum(E const &expr, long axis, dtype d = dtype()) ->
      typename std::enable_if<E::value == 1,
                              decltype(partial_sum<Op, E, dtype>(expr))>::type;

  template <class Op, class E, class dtype = result_dtype<Op, E>>
  using partial_sum_type =
      types::ndarray<typename dtype::type, types::array<long, E::value>>;
  template <class Op, class E, class dtype = result_dtype<Op, E>>
  using partial_sum_type2 =
      types::ndarray<typename dtype::type, types::array<long, E::value - 1>>;

  template <class Op, class E, class dtype = result_dtype<Op, E>>
  typename std::enable_if<E::value != 1, partial_sum_type<Op, E, dtype>>::type
  partial_sum(E const &expr, long axis, dtype d = dtype());
}
PYTHONIC_NS_END

#endif

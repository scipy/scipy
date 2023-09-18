#ifndef PYTHONIC_INCLUDE_NUMPY_FROMITER_HPP
#define PYTHONIC_INCLUDE_NUMPY_FROMITER_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class Iterable, class dtype = functor::float64>
  types::ndarray<typename std::remove_cv<typename std::remove_reference<
                     Iterable>::type>::type::value_type,
                 types::pshape<long>>
  fromiter(Iterable &&iterable, dtype d = dtype(), long count = -1);

  DEFINE_FUNCTOR(pythonic::numpy, fromiter);
}
PYTHONIC_NS_END

#endif

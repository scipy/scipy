#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/nested_container.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class pS, class dtype = functor::float64>
  types::ndarray<typename dtype::type, sutils::shape_t<pS>>
  ndarray(pS const &shape, dtype d = dtype());

  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<long>>
  ndarray(long size, dtype d = dtype());

  template <long N, class dtype = functor::float64>
  types::ndarray<typename dtype::type,
                 types::pshape<std::integral_constant<long, N>>>
  ndarray(std::integral_constant<long, N>, dtype d = dtype());

  DEFINE_FUNCTOR(pythonic::numpy, ndarray);
}
PYTHONIC_NS_END

#endif

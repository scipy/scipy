#ifndef PYTHONIC_INCLUDE_NUMPY_INDICES_HPP
#define PYTHONIC_INCLUDE_NUMPY_INDICES_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/numpy/int64.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class pS, class dtype = functor::int64>
  types::ndarray<
      typename dtype::type,
      sutils::push_front_t<
          pS, std::integral_constant<long, std::tuple_size<pS>::value>>>
  indices(pS const &shape, dtype d = dtype());

  DEFINE_FUNCTOR(pythonic::numpy, indices);
}
PYTHONIC_NS_END

#endif

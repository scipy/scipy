#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_RESHAPE_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_RESHAPE_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace ndarray
  {
    template <class T, class pS, class NpS>
    typename std::enable_if<!std::is_integral<NpS>::value,
                            types::ndarray<T, NpS>>::type
    reshape(types::ndarray<T, pS> const &expr, NpS const &new_shape);
    template <class T, class pS, class NpS>
    typename std::enable_if<std::is_integral<NpS>::value,
                            types::ndarray<T, types::pshape<long>>>::type
    reshape(types::ndarray<T, pS> const &expr, NpS const &new_shape);

    template <class T, class pS, class S0, class S1, class... S>
    auto reshape(types::ndarray<T, pS> const &expr, S0 i0, S1 i1,
                 S const &... indices)
        -> decltype(reshape(expr,
                            types::pshape<S0, S1, S...>{i0, i1, indices...}));

    NUMPY_EXPR_TO_NDARRAY0_DECL(reshape);

    DEFINE_FUNCTOR(pythonic::numpy::ndarray, reshape);
  }
}
PYTHONIC_NS_END

#endif

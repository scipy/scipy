#ifndef PYTHONIC_INCLUDE_NUMPY_DELETE_HPP
#define PYTHONIC_INCLUDE_NUMPY_DELETE_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::pshape<long>>
  delete_(types::ndarray<T, pS> const &a, long index,
          types::none_type axis = builtins::None);

  template <class T, class pS, class I>
  typename std::enable_if<!std::is_scalar<I>::value,
                          types::ndarray<T, types::pshape<long>>>::type
  delete_(types::ndarray<T, pS> const &in, I const &indices,
          types::none_type axis = builtins::None);

  NUMPY_EXPR_TO_NDARRAY0_DECL(delete_);
  DEFINE_FUNCTOR(pythonic::numpy, delete_);
}
PYTHONIC_NS_END

#endif

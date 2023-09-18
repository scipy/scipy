#ifndef PYTHONIC_INCLUDE_NUMPY_FLOAT128_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLOAT128_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {

    long double float128();
    template <class V>
    long double float128(V v);
  }

#define NUMPY_NARY_FUNC_NAME float128
#define NUMPY_NARY_FUNC_SYM details::float128
#define NUMPY_NARY_EXTRA_METHOD using type = long double;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

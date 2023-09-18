#ifndef PYTHONIC_INCLUDE_NUMPY_FLOAT64_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLOAT64_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {

    double float64();
    template <class V>
    double float64(V v);
  }

#define NUMPY_NARY_FUNC_NAME float64
#define NUMPY_NARY_FUNC_SYM details::float64
#define NUMPY_NARY_EXTRA_METHOD using type = double;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

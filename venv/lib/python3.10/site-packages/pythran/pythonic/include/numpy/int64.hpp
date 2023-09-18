#ifndef PYTHONIC_INCLUDE_NUMPY_INT64_HPP
#define PYTHONIC_INCLUDE_NUMPY_INT64_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    int64_t int64();
    template <class V>
    int64_t int64(V v);
  }

#define NUMPY_NARY_FUNC_NAME int64
#define NUMPY_NARY_FUNC_SYM details::int64
#define NUMPY_NARY_EXTRA_METHOD using type = int64_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

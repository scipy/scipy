#ifndef PYTHONIC_INCLUDE_NUMPY_INT32_HPP
#define PYTHONIC_INCLUDE_NUMPY_INT32_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    int32_t int32();
    template <class V>
    int32_t int32(V v);
  }

#define NUMPY_NARY_FUNC_NAME int32
#define NUMPY_NARY_FUNC_SYM details::int32
#define NUMPY_NARY_EXTRA_METHOD using type = int32_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

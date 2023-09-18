#ifndef PYTHONIC_INCLUDE_NUMPY_FLOAT32_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLOAT32_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    float float32();
    template <class V>
    float float32(V v);
  }

#define NUMPY_NARY_FUNC_NAME float32
#define NUMPY_NARY_FUNC_SYM details::float32
#define NUMPY_NARY_EXTRA_METHOD using type = float;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

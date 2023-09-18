#ifndef PYTHONIC_NUMPY_FLOAT32_HPP
#define PYTHONIC_NUMPY_FLOAT32_HPP

#include "pythonic/include/numpy/float32.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline float float32()
    {
      return float();
    }

    template <class V>
    float float32(V v)
    {
      return static_cast<float>(v);
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME float32
#define NUMPY_NARY_FUNC_SYM details::float32
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

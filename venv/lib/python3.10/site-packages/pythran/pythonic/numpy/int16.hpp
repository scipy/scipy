#ifndef PYTHONIC_NUMPY_INT16_HPP
#define PYTHONIC_NUMPY_INT16_HPP

#include "pythonic/include/numpy/int16.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline int16_t int16()
    {
      return int16_t();
    }

    template <class V>
    int16_t int16(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME int16
#define NUMPY_NARY_FUNC_SYM details::int16
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

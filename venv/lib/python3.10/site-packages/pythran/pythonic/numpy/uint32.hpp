#ifndef PYTHONIC_NUMPY_UINT32_HPP
#define PYTHONIC_NUMPY_UINT32_HPP

#include "pythonic/include/numpy/uint32.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline uint32_t uint32()
    {
      return uint32_t();
    }

    template <class V>
    uint32_t uint32(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME uint32
#define NUMPY_NARY_FUNC_SYM details::uint32
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

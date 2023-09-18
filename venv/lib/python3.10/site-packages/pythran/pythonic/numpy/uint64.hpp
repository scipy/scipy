#ifndef PYTHONIC_NUMPY_UINT64_HPP
#define PYTHONIC_NUMPY_UINT64_HPP

#include "pythonic/include/numpy/uint64.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline uint64_t uint64()
    {
      return uint64_t();
    }

    template <class V>
    uint64_t uint64(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME uint64
#define NUMPY_NARY_FUNC_SYM details::uint64
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

#ifndef PYTHONIC_INCLUDE_NUMPY_UINT64_HPP
#define PYTHONIC_INCLUDE_NUMPY_UINT64_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    uint64_t uint64();
    template <class V>
    uint64_t uint64(V v);
  }

#define NUMPY_NARY_FUNC_NAME uint64
#define NUMPY_NARY_FUNC_SYM details::uint64
#define NUMPY_NARY_EXTRA_METHOD using type = uint64_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

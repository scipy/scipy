#ifndef PYTHONIC_INCLUDE_NUMPY_UINT32_HPP
#define PYTHONIC_INCLUDE_NUMPY_UINT32_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    uint32_t uint32();
    template <class V>
    uint32_t uint32(V v);
  }

#define NUMPY_NARY_FUNC_NAME uint32
#define NUMPY_NARY_FUNC_SYM details::uint32
#define NUMPY_NARY_EXTRA_METHOD using type = uint32_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

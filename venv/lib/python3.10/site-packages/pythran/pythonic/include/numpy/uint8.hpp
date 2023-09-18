#ifndef PYTHONIC_INCLUDE_NUMPY_UINT8_HPP
#define PYTHONIC_INCLUDE_NUMPY_UINT8_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    uint8_t uint8();
    template <class V>
    uint8_t uint8(V v);
  }

#define NUMPY_NARY_FUNC_NAME uint8
#define NUMPY_NARY_FUNC_SYM details::uint8
#define NUMPY_NARY_EXTRA_METHOD using type = uint8_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

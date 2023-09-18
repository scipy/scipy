#ifndef PYTHONIC_INCLUDE_NUMPY_UINT16_HPP
#define PYTHONIC_INCLUDE_NUMPY_UINT16_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    uint16_t uint16();
    template <class V>
    uint16_t uint16(V v);
  }

#define NUMPY_NARY_FUNC_NAME uint16
#define NUMPY_NARY_FUNC_SYM details::uint16
#define NUMPY_NARY_EXTRA_METHOD using type = uint16_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

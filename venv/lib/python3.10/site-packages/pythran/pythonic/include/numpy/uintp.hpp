#ifndef PYTHONIC_INCLUDE_NUMPY_UINTP_HPP
#define PYTHONIC_INCLUDE_NUMPY_UINTP_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    uintptr_t uintp();
    template <class V>
    uintptr_t uintp(V v);
  }

#define NUMPY_NARY_FUNC_NAME uintp
#define NUMPY_NARY_FUNC_SYM details::uintp
#define NUMPY_NARY_EXTRA_METHOD using type = uintptr_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

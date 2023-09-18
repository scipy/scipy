#ifndef PYTHONIC_INCLUDE_NUMPY_INTP_HPP
#define PYTHONIC_INCLUDE_NUMPY_INTP_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    intptr_t intp();
    template <class V>
    intptr_t intp(V v);
  }

#define NUMPY_NARY_FUNC_NAME intp
#define NUMPY_NARY_FUNC_SYM details::intp
#define NUMPY_NARY_EXTRA_METHOD using type = intptr_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

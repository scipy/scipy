#ifndef PYTHONIC_INCLUDE_NUMPY_LONGLONG_HPP
#define PYTHONIC_INCLUDE_NUMPY_LONGLONG_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    long long longlong();
    template <class V>
    long long longlong(V v);
  }

#define NUMPY_NARY_FUNC_NAME longlong
#define NUMPY_NARY_FUNC_SYM details::longlong
#define NUMPY_NARY_EXTRA_METHOD using type = long long;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

#ifndef PYTHONIC_INCLUDE_NUMPY_BOOL_HPP
#define PYTHONIC_INCLUDE_NUMPY_BOOL_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    bool bool_();
    template <class V>
    bool bool_(V v);
  }

#define NUMPY_NARY_FUNC_NAME bool_
#define NUMPY_NARY_FUNC_SYM details::bool_
#define NUMPY_NARY_EXTRA_METHOD using type = bool;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

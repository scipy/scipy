#ifndef PYTHONIC_INCLUDE_NUMPY_BYTE_HPP
#define PYTHONIC_INCLUDE_NUMPY_BYTE_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    char byte();
    template <class V>
    char byte(V v);
  }

#define NUMPY_NARY_FUNC_NAME byte
#define NUMPY_NARY_FUNC_SYM details::byte
#define NUMPY_NARY_EXTRA_METHOD using type = char;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

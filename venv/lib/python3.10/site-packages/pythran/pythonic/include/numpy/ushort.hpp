#ifndef PYTHONIC_INCLUDE_NUMPY_USHORT_HPP
#define PYTHONIC_INCLUDE_NUMPY_USHORT_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    unsigned short ushort();
    template <class V>
    unsigned short ushort(V v);
  }

#define NUMPY_NARY_FUNC_NAME ushort
#define NUMPY_NARY_FUNC_SYM details::ushort
#define NUMPY_NARY_EXTRA_METHOD using type = unsigned short;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

#ifndef PYTHONIC_NUMPY_BYTE_HPP
#define PYTHONIC_NUMPY_BYTE_HPP

#include "pythonic/include/numpy/byte.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline char byte()
    {
      return {};
    }

    template <class V>
    char byte(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME byte
#define NUMPY_NARY_FUNC_SYM details::byte
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

#ifndef PYTHONIC_NUMPY_USHORT_HPP
#define PYTHONIC_NUMPY_USHORT_HPP

#include "pythonic/include/numpy/ushort.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline unsigned short ushort()
    {
      return {};
    }

    template <class V>
    unsigned short ushort(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME ushort
#define NUMPY_NARY_FUNC_SYM details::ushort
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

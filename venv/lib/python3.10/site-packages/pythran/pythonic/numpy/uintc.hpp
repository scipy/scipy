#ifndef PYTHONIC_NUMPY_UINTC_HPP
#define PYTHONIC_NUMPY_UINTC_HPP

#include "pythonic/include/numpy/uintc.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline unsigned uintc()
    {
      return {};
    }

    template <class V>
    unsigned uintc(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME uintc
#define NUMPY_NARY_FUNC_SYM details::uintc
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

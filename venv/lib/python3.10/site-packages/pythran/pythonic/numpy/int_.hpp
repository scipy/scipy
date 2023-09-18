#ifndef PYTHONIC_NUMPY_INT__HPP
#define PYTHONIC_NUMPY_INT__HPP

#include "pythonic/include/numpy/int_.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline long int_()
    {
      return {};
    }

    template <class V>
    long int_(V v)
    {
      return static_cast<long>(v);
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME int_
#define NUMPY_NARY_FUNC_SYM details::int_
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

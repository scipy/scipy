#ifndef PYTHONIC_NUMPY_FLOAT64_HPP
#define PYTHONIC_NUMPY_FLOAT64_HPP

#include "pythonic/include/numpy/float64.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {

    inline double float64()
    {
      return double();
    }

    template <class V>
    double float64(V v)
    {
      return static_cast<double>(v);
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME float64
#define NUMPY_NARY_FUNC_SYM details::float64
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

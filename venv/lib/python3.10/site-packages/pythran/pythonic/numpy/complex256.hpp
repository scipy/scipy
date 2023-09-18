#ifndef PYTHONIC_NUMPY_COMPLEX256_HPP
#define PYTHONIC_NUMPY_COMPLEX256_HPP

#include "pythonic/include/numpy/complex256.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {

    inline std::complex<long double> complex256()
    {
      return {};
    }

    template <class V>
    std::complex<long double> complex256(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME complex256
#define NUMPY_NARY_FUNC_SYM details::complex256
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif

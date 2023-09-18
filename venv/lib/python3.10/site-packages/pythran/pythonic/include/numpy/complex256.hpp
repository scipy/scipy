#ifndef PYTHONIC_INCLUDE_NUMPY_COMPLEX256_HPP
#define PYTHONIC_INCLUDE_NUMPY_COMPLEX256_HPP

#include "pythonic/include/types/complex.hpp"

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    std::complex<long double> complex256();
    template <class V>
    std::complex<long double> complex256(V v);
  }

#define NUMPY_NARY_FUNC_NAME complex256
#define NUMPY_NARY_FUNC_SYM details::complex256
#define NUMPY_NARY_EXTRA_METHOD using type = std::complex<long double>;
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

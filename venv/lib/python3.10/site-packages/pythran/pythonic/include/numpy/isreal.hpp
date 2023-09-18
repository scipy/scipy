#ifndef PYTHONIC_INCLUDE_NUMPY_ISREAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISREAL_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/types/traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace wrapper
  {
    template <class I>
    typename std::enable_if<types::is_complex<I>::value, bool>::type
    isreal(I const &a);

    template <class I>
    typename std::enable_if<!types::is_complex<I>::value, bool>::type
    isreal(I const &a);
  }

#define NUMPY_NARY_FUNC_NAME isreal
#define NUMPY_NARY_FUNC_SYM wrapper::isreal
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

#ifndef PYTHONIC_NUMPY_ISCOMPLEX_HPP
#define PYTHONIC_NUMPY_ISCOMPLEX_HPP

#include "pythonic/include/numpy/iscomplex.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/numpy_traits.hpp"
#include "pythonic/types/traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace wrapper
  {
    template <class I>
    typename std::enable_if<types::is_complex<I>::value, bool>::type
    iscomplex(I const &a)
    {
      return a.imag() != 0.;
    }

    template <class I>
    constexpr typename std::enable_if<!types::is_complex<I>::value, bool>::type
    iscomplex(I const &a)
    {
      return false;
    }
  }

#define NUMPY_NARY_FUNC_NAME iscomplex
#define NUMPY_NARY_FUNC_SYM wrapper::iscomplex
#include "pythonic/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif

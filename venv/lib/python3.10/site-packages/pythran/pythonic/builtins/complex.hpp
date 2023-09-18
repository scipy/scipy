#ifndef PYTHONIC_BUILTIN_COMPLEX_HPP
#define PYTHONIC_BUILTIN_COMPLEX_HPP

#include "pythonic/include/builtins/complex.hpp"

#include "pythonic/types/complex.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {
    inline complex::type complex::operator()(double v0, double v1) const
    {
      return {v0, v1};
    }
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#endif

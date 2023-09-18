#ifndef PYTHONIC_INCLUDE_BUILTIN_COMPLEX_CONJUGATE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_COMPLEX_CONJUGATE_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/numpy/conjugate.hpp"

PYTHONIC_NS_BEGIN
namespace builtins
{
  namespace complex
  {
    USING_FUNCTOR(conjugate, numpy::functor::conjugate);
  }
}
PYTHONIC_NS_END
#endif

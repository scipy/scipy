#ifndef PYTHONIC_INCLUDE_BUILTIN_ABS_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ABS_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/numpy/abs.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  // FIXME np.abs accept any iterator while builtins.abs only accept
  // numeric types && numpy.array
  USING_FUNCTOR(abs, numpy::functor::abs);
}
PYTHONIC_NS_END

#endif

#ifndef PYTHONIC_INCLUDE_BUILTIN_ALL_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ALL_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class Iterable>
  bool all(Iterable &&s);

  DEFINE_FUNCTOR(pythonic::builtins, all);
}
PYTHONIC_NS_END

#endif

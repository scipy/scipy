#ifndef PYTHONIC_INCLUDE_BUILTIN_HEX_HPP
#define PYTHONIC_INCLUDE_BUILTIN_HEX_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/str.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class T>
  types::str hex(T const &v);

  DEFINE_FUNCTOR(pythonic::builtins, hex);
}
PYTHONIC_NS_END

#endif

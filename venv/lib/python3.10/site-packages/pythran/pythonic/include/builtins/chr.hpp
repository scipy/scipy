#ifndef PYTHONIC_INCLUDE_BUILTIN_CHR_HPP
#define PYTHONIC_INCLUDE_BUILTIN_CHR_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/str.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  types::str chr(T const &v);

  DEFINE_FUNCTOR(pythonic::builtins, chr);
}
PYTHONIC_NS_END

#endif

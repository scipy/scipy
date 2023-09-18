#ifndef PYTHONIC_INCLUDE_BUILTIN_BIN_HPP
#define PYTHONIC_INCLUDE_BUILTIN_BIN_HPP

#include "pythonic/include/utils/functor.hpp"

#include "pythonic/include/types/str.hpp"

#include <type_traits>

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  typename std::enable_if<std::is_scalar<T>::value, types::str>::type
  bin(T const &v);

  DEFINE_FUNCTOR(pythonic::builtins, bin);
}
PYTHONIC_NS_END

#endif

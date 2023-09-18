#ifndef PYTHONIC_BUILTIN_PYTHRAN_STATICIFBREAK_HPP
#define PYTHONIC_BUILTIN_PYTHRAN_STATICIFBREAK_HPP

#include "pythonic/include/builtins/pythran/StaticIfBreak.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/types/static_if.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    template <class T>
    types::StaticIfBreak<T> StaticIfBreak(T const &arg)
    {
      return {arg};
    }
  }
}
PYTHONIC_NS_END

#endif

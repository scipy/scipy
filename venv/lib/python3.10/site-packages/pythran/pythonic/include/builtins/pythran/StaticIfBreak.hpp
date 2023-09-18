#ifndef PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATICIFBREAK_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATICIFBREAK_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/static_if.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    template <class T>
    types::StaticIfBreak<T> StaticIfBreak(T const &arg);

    DEFINE_FUNCTOR(pythonic::builtins::pythran, StaticIfBreak);
  }
}

PYTHONIC_NS_END

#endif

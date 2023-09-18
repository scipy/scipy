#ifndef PYTHONIC_BUILTIN_PYTHRAN_STATICIFNORETURN_HPP
#define PYTHONIC_BUILTIN_PYTHRAN_STATICIFNORETURN_HPP

#include "pythonic/include/builtins/pythran/StaticIfNoReturn.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/types/static_if.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    template <class T>
    types::StaticIfNoReturn<T> StaticIfNoReturn(T const &arg)
    {
      return {arg};
    }
  }
}
PYTHONIC_NS_END

#endif

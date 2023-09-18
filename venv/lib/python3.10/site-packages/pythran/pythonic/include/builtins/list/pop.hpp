#ifndef PYTHONIC_INCLUDE_BUILTIN_LIST_POP_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LIST_POP_HPP

#include "pythonic/include/__dispatch__/pop.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {
    USING_FUNCTOR(pop, pythonic::__dispatch__::functor::pop);
  }
}
PYTHONIC_NS_END
#endif

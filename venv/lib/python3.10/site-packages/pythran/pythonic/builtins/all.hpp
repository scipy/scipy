#ifndef PYTHONIC_BUILTIN_ALL_HPP
#define PYTHONIC_BUILTIN_ALL_HPP

#include "pythonic/utils/functor.hpp"
#include "pythonic/include/builtins/all.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class Iterable>
  bool all(Iterable &&s)
  {
    auto iend = s.end();
    for (auto iter = s.begin(); iter != iend; ++iter)
      if (!*iter)
        return false;
    return true;
  }
}
PYTHONIC_NS_END

#endif

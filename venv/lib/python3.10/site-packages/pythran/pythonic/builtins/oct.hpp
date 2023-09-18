#ifndef PYTHONIC_BUILTIN_OCT_HPP
#define PYTHONIC_BUILTIN_OCT_HPP

#include "pythonic/include/builtins/oct.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

#include <sstream>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class T>
  types::str oct(T const &v)
  {
    std::ostringstream oss;
    oss <<
#if defined(__PYTHRAN__) && __PYTHRAN__ == 3
        "0o"
#else
        '0'
#endif
        << std::oct << v;
    return oss.str();
  }
}
PYTHONIC_NS_END

#endif

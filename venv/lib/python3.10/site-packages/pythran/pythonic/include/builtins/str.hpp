#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    template <class T>
    types::str str(T const &t);

    inline types::str str();
    inline types::str str(bool b);
    inline types::str str(long value);
    inline types::str str(double l);
  }

  DEFINE_FUNCTOR(pythonic::builtins::anonymous, str);
}
PYTHONIC_NS_END

#endif

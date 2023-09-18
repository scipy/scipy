#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_ISALPHA_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_ISALPHA_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/str.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    bool isalpha(types::str const &s);

    DEFINE_FUNCTOR(pythonic::builtins::str, isalpha);
  }
}
PYTHONIC_NS_END
#endif

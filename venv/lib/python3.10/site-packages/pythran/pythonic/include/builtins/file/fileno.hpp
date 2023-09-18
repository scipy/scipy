#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_FILENO_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_FILENO_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    long fileno(types::file const &f);

    DEFINE_FUNCTOR(pythonic::builtins::file, fileno);
  }
}
PYTHONIC_NS_END
#endif

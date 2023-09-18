#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_FLUSH_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_FLUSH_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    void flush(types::file &f);
    void flush(types::file &&f);

    DEFINE_FUNCTOR(pythonic::builtins::file, flush);
  }
}
PYTHONIC_NS_END
#endif

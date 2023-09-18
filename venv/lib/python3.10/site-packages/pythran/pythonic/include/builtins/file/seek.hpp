#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_SEEK_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_SEEK_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/file.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    void seek(types::file &f, long offset);
    void seek(types::file &&f, long offset);
    void seek(types::file &f, long offset, long whence);
    void seek(types::file &&f, long offset, long whence);

    DEFINE_FUNCTOR(pythonic::builtins::file, seek);
  }
}
PYTHONIC_NS_END
#endif

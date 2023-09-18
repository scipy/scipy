#ifndef PYTHONIC_INCLUDE_BUILTIN_LIST_INSERT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LIST_INSERT_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/NoneType.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {

    template <class T, class F>
    types::none_type insert(types::list<T> &seq, long n, F &&value);

    DEFINE_FUNCTOR(pythonic::builtins::list, insert);
  }
}
PYTHONIC_NS_END
#endif

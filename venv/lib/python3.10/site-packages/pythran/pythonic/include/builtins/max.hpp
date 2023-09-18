#ifndef PYTHONIC_INCLUDE_BUILTIN_MAX_HPP
#define PYTHONIC_INCLUDE_BUILTIN_MAX_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/operator_/gt.hpp"
#include "pythonic/include/builtins/minmax.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class... Types>
  auto max(Types &&... values)
      -> decltype(details::minmax(operator_::functor::lt{},
                                  std::forward<Types>(values)...));

  DEFINE_FUNCTOR(pythonic::builtins, max);
}
PYTHONIC_NS_END

#endif

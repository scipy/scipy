#ifndef PYTHONIC_OPERATOR_IS_HPP
#define PYTHONIC_OPERATOR_IS_HPP

#include "pythonic/include/operator_/is_.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/builtins/id.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto is_(A &&a, B &&b) -> decltype(builtins::id(std::forward<A>(a)) ==
                                     builtins::id(std::forward<B>(b)))
  {
    return builtins::id(std::forward<A>(a)) == builtins::id(std::forward<B>(b));
  }
}
PYTHONIC_NS_END

#endif

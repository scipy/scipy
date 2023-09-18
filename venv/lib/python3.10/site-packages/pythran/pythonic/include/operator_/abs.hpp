#ifndef PYTHONIC_INCLUDE_OPERATOR_ABS_HPP
#define PYTHONIC_INCLUDE_OPERATOR_ABS_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/builtins/abs.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  USING_FUNCTOR(abs, builtins::functor::abs);
}
PYTHONIC_NS_END

#endif

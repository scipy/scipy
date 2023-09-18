#ifndef PYTHONIC_INCLUDE_OPERATOR_FLOORDIV_HPP
#define PYTHONIC_INCLUDE_OPERATOR_FLOORDIV_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/numpy/floor_divide.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  USING_FUNCTOR(floordiv, numpy::functor::floor_divide);
}
PYTHONIC_NS_END

#endif

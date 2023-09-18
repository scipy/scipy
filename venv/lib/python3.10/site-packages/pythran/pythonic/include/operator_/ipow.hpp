#ifndef PYTHONIC_INCLUDE_OPERATOR_IPOW_HPP
#define PYTHONIC_INCLUDE_OPERATOR_IPOW_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/builtins/pow.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  A ipow(A const &a, B &&b);
  template <class A, class B>
  A &ipow(A &a, B &&b);

  DEFINE_FUNCTOR(pythonic::operator_, ipow);
}
PYTHONIC_NS_END

#endif

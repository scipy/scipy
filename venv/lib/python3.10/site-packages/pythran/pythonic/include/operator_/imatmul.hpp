#ifndef PYTHONIC_INCLUDE_OPERATOR_IMATMUL_HPP
#define PYTHONIC_INCLUDE_OPERATOR_IMATMUL_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/numpy/dot.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  A imatmul(A const &a, B &&b);
  template <class A, class B>
  A imatmul(A &a, B &&b);

  DEFINE_FUNCTOR(pythonic::operator_, imatmul);
}
PYTHONIC_NS_END

#endif

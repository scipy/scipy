#ifndef PYTHONIC_INCLUDE_OPERATOR_MUL_HPP
#define PYTHONIC_INCLUDE_OPERATOR_MUL_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/operator_/overloads.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto mul(A &&a, B &&b) -> decltype(std::forward<A>(a) * std::forward<B>(b));

  DEFINE_ALL_OPERATOR_OVERLOADS_DECL(mul, *)

  DEFINE_FUNCTOR(pythonic::operator_, mul);
}
PYTHONIC_NS_END

#endif

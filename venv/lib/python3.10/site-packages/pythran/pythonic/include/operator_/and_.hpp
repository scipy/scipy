#ifndef PYTHONIC_INCLUDE_OPERATOR_AND_HPP
#define PYTHONIC_INCLUDE_OPERATOR_AND_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/operator_/overloads.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto and_(A &&a, B &&b) -> decltype(std::forward<A>(a) & std::forward<B>(b));

  DEFINE_ALL_OPERATOR_OVERLOADS_DECL(and_, &)

  DEFINE_FUNCTOR(pythonic::operator_, and_);
}
PYTHONIC_NS_END

#endif

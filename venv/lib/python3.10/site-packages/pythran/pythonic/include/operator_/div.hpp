#ifndef PYTHONIC_INCLUDE_OPERATOR_DIV_HPP
#define PYTHONIC_INCLUDE_OPERATOR_DIV_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/operator_/overloads.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto div(A &&a, B &&b) // for ndarrays
      -> typename std::enable_if<
          !std::is_fundamental<typename std::decay<A>::type>::value ||
              !std::is_fundamental<typename std::decay<B>::type>::value,
          decltype(std::forward<A>(a) / std::forward<B>(b))>::type;

  double div(double a, double b);

  DEFINE_FUNCTOR(pythonic::operator_, div);
}
PYTHONIC_NS_END

#endif

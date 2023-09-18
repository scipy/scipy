#ifndef PYTHONIC_OPERATOR_DIV_HPP
#define PYTHONIC_OPERATOR_DIV_HPP

#include "pythonic/include/operator_/div.hpp"

#include "pythonic/operator_/overloads.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto div(A &&a, B &&b) // for ndarrays
      -> typename std::enable_if<
          !std::is_fundamental<typename std::decay<A>::type>::value ||
              !std::is_fundamental<typename std::decay<B>::type>::value,
          decltype(std::forward<A>(a) / std::forward<B>(b))>::type
  {
    return std::forward<A>(a) / std::forward<B>(b);
  }

  inline double div(double a, double b)
  {
    assert(b != 0 && "divide by zero");
    return a / b;
  }
} // namespace operator_
PYTHONIC_NS_END

#endif

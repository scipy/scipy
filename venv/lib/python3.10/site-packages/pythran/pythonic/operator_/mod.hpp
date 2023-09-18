#ifndef PYTHONIC_OPERATOR_MOD_HPP
#define PYTHONIC_OPERATOR_MOD_HPP

#include "pythonic/include/operator_/mod.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto mod(A &&a, B &&b) -> typename std::enable_if<
      std::is_fundamental<typename std::decay<A>::type>::value &&
          std::is_fundamental<typename std::decay<B>::type>::value,
      decltype(std::forward<A>(a) % std::forward<B>(b))>::type
  {
    auto t = std::forward<A>(a) % b;
    return t < 0 ? (t + b) : t;
  }

  inline double mod(double a, long b)
  {
    auto t = std::fmod(a, double(b));
    return t < 0 ? (t + b) : t;
  }

  inline double mod(double a, double b)
  {
    auto t = std::fmod(a, b);
    return t < 0 ? (t + b) : t;
  }

  template <class A, class B>
  auto mod(A &&a, B &&b) // for ndarrays
      -> typename std::enable_if<
          !std::is_fundamental<typename std::decay<A>::type>::value ||
              !std::is_fundamental<typename std::decay<B>::type>::value,
          decltype(std::forward<A>(a) % std::forward<B>(b))>::type
  {
    return std::forward<A>(a) % std::forward<B>(b);
  }
}
PYTHONIC_NS_END

#endif

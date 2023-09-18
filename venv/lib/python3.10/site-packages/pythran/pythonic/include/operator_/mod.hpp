#ifndef PYTHONIC_INCLUDE_OPERATOR_MOD_HPP
#define PYTHONIC_INCLUDE_OPERATOR_MOD_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto mod(A &&a, B &&b) -> typename std::enable_if<
      std::is_fundamental<typename std::decay<A>::type>::value &&
          std::is_fundamental<typename std::decay<B>::type>::value,
      decltype(std::forward<A>(a) % std::forward<B>(b))>::type;

  inline double mod(double a, long b);

  inline double mod(double a, double b);

  template <class A, class B>
  auto mod(A &&a, B &&b) // for ndarrays
      -> typename std::enable_if<
          !std::is_fundamental<typename std::decay<A>::type>::value ||
              !std::is_fundamental<typename std::decay<B>::type>::value,
          decltype(std::forward<A>(a) % std::forward<B>(b))>::type;

  DEFINE_FUNCTOR(pythonic::operator_, mod);
}
PYTHONIC_NS_END

#endif

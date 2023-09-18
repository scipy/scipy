#ifndef PYTHONIC_INCLUDE_OPERATOR_IMAX_HPP
#define PYTHONIC_INCLUDE_OPERATOR_IMAX_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/numpy/maximum.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto imax(A &&a, B &&b) -> typename std::enable_if<
      std::is_const<A>::value || !std::is_assignable<A, B>::value,
      decltype(numpy::functor::maximum{}(std::forward<A>(a),
                                         std::forward<B>(b)))>::type;

  template <class A, class B>
  auto imax(A &&a, B &&b) -> typename std::enable_if<
      !std::is_const<A>::value && std::is_assignable<A, B>::value,
      decltype(a = numpy::functor::maximum{}(std::forward<A>(a),
                                             std::forward<B>(b)))>::type;

  DEFINE_FUNCTOR(pythonic::operator_, imax);
}
PYTHONIC_NS_END

#endif

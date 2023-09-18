#ifndef PYTHONIC_INCLUDE_BUILTIN_LEN_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LEN_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/yield.hpp"

#include <tuple>
#include <iterator>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class... Types>
  long len(std::tuple<Types...> const &);

  template <class T>
  typename std::enable_if<types::has_size<T>::value, long>::type
  len(T const &t);

  DEFINE_FUNCTOR(pythonic::builtins, len);
}
PYTHONIC_NS_END
#endif

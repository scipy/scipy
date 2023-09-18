#ifndef PYTHONIC_BUILTIN_LEN_HPP
#define PYTHONIC_BUILTIN_LEN_HPP

#include "pythonic/include/builtins/len.hpp"

#include "pythonic/types/traits.hpp"
#include "pythonic/utils/functor.hpp"

#include <tuple>
#include <iterator>

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class... Types>
  long len(std::tuple<Types...> const &)
  {
    return sizeof...(Types);
  }

  template <class T>
  typename std::enable_if<types::has_size<T>::value, long>::type len(T const &t)
  {
    return t.size();
  }
}
PYTHONIC_NS_END
#endif

#ifndef PYTHONIC_BUILTIN_DICT_HPP
#define PYTHONIC_BUILTIN_DICT_HPP

#include "pythonic/include/builtins/dict.hpp"

#include "pythonic/types/dict.hpp"
#include "pythonic/utils/functor.hpp"

#include <tuple>
#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    inline types::empty_dict dict()
    {
      return types::empty_dict();
    }

    template <class K, class V>
    types::dict<K, V> dict(types::dict<K, V> const &other)
    {
      return other.copy();
    }

    template <class Iterable>
    auto dict(Iterable &&iterable) -> types::dict<
        typename std::decay<decltype(std::get<0>(*iterable.begin()))>::type,
        typename std::decay<decltype(std::get<1>(*iterable.begin()))>::type>
    {
      types::dict<
          typename std::decay<decltype(std::get<0>(*iterable.begin()))>::type,
          typename std::decay<decltype(std::get<1>(*iterable.begin()))>::type>
          out = types::empty_dict();
      for (auto const &i : iterable)
        out[std::get<0>(i)] = std::get<1>(i);
      return out;
    }
  } // namespace anonymous
} // namespace builtins
PYTHONIC_NS_END

#endif

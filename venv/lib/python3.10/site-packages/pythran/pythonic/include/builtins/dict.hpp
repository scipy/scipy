#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_HPP

#include "pythonic/include/types/dict.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <tuple>
#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    types::empty_dict dict();

    template <class K, class V>
    types::dict<K, V> dict(types::dict<K, V> const &);

    template <class Iterable>
    auto dict(Iterable &&iterable) -> types::dict<
        typename std::decay<decltype(std::get<0>(*iterable.begin()))>::type,
        typename std::decay<decltype(std::get<1>(*iterable.begin()))>::type>;
  }

  DEFINE_FUNCTOR(pythonic::builtins::anonymous, dict);
}
PYTHONIC_NS_END

#endif

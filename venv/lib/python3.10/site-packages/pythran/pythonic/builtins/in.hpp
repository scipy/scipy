#ifndef PYTHONIC_BUILTIN_IN_HPP
#define PYTHONIC_BUILTIN_IN_HPP

#include "pythonic/include/builtins/in.hpp"

#include "pythonic/types/traits.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace details
{
  template <bool has_contains>
  struct in {
    template <class T, class V>
    bool operator()(T &&t, V const &v) const;
  };

  template <>
  template <class T, class V>
  bool in<false>::operator()(T &&t, V const &v) const
  {
    return std::find(t.begin(), t.end(), v) != t.end();
  }

  template <>
  template <class T, class V>
  bool in<true>::operator()(T &&t, V const &v) const
  {
    return t.contains(v);
  }
}

template <class T, class V>
bool in(T &&t, V const &v)
{
  using RT =
      typename std::remove_cv<typename std::remove_reference<T>::type>::type;
  static bool constexpr has_contains = types::has_contains<RT, V>::value;
  return details::in<has_contains>()(std::forward<T>(t), v);
}
PYTHONIC_NS_END

#endif

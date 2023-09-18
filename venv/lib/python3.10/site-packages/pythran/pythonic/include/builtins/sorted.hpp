#ifndef PYTHONIC_INCLUDE_BUILTIN_SORTED_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SORTED_HPP

#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class Iterable>
  types::list<typename std::remove_cv<typename std::iterator_traits<
      typename std::decay<Iterable>::type::iterator>::value_type>::type>
  sorted(Iterable &&seq);

  template <class Iterable, class Key>
  types::list<typename std::remove_cv<typename std::iterator_traits<
      typename std::decay<Iterable>::type::iterator>::value_type>::type>
  sorted(Iterable &&seq, Key const &key, bool reverse = false);

  template <class Iterable>
  types::list<typename std::remove_cv<typename std::iterator_traits<
      typename std::decay<Iterable>::type::iterator>::value_type>::type>
  sorted(Iterable &&seq, types::none_type const &key, bool reverse = false);

  DEFINE_FUNCTOR(pythonic::builtins, sorted);
}
PYTHONIC_NS_END

#endif

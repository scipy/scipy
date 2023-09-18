#ifndef PYTHONIC_INCLUDE_BUILTIN_TUPLE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_TUPLE_HPP

#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/types/dynamic_tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class... Types>
  std::tuple<Types...> tuple(std::tuple<Types...> const &t);

  template <class Iterable> /* this is far from perfect, but how to cope with
                               the
                               difference between python tuples && c++ ones ? */
      typename std::enable_if <
      types::len_of<typename std::remove_cv<
          typename std::remove_reference<Iterable>::type>::type>::
          value<0, types::dynamic_tuple<typename std::iterator_traits<
                       typename std::remove_cv<typename std::remove_reference<
                           Iterable>::type>::type::iterator>::value_type>>::type
          tuple(Iterable &&i);

  template <
      class StaticIterable> /* specialization if we are capable to statically
                               compute the size of the input */
  typename std::enable_if<
      types::len_of<typename std::remove_cv<typename std::remove_reference<
          StaticIterable>::type>::type>::value >= 0,
      types::array<
          typename std::iterator_traits<
              typename std::remove_cv<typename std::remove_reference<
                  StaticIterable>::type>::type::iterator>::value_type,
          types::len_of<typename std::remove_cv<typename std::remove_reference<
              StaticIterable>::type>::type>::value>>::type
  tuple(StaticIterable &&i);

  DEFINE_FUNCTOR(pythonic::builtins, tuple);
}
PYTHONIC_NS_END

#endif

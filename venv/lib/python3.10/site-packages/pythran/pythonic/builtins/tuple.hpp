#ifndef PYTHONIC_BUILTIN_TUPLE_HPP
#define PYTHONIC_BUILTIN_TUPLE_HPP

#include "pythonic/include/builtins/tuple.hpp"

#include "pythonic/types/dynamic_tuple.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class... Types>
  std::tuple<Types...> tuple(std::tuple<Types...> const &t)
  {
    return t;
  }

  template <class Iterable>
      /* this is far from perfect, but how to cope with the
         difference between python tuples && c++ ones ? */
      typename std::enable_if <
      types::len_of<typename std::remove_cv<
          typename std::remove_reference<Iterable>::type>::type>::
          value<0, types::dynamic_tuple<typename std::iterator_traits<
                       typename std::remove_cv<typename std::remove_reference<
                           Iterable>::type>::type::iterator>::value_type>>::type
          tuple(Iterable &&i)
  {
    return {i.begin(), i.end()};
  }

  template <class StaticIterable>
  /* specialization if we are capable to statically compute the size of the
     input */
  typename std::enable_if<
      types::len_of<typename std::remove_cv<typename std::remove_reference<
          StaticIterable>::type>::type>::value >= 0,
      types::array<
          typename std::iterator_traits<
              typename std::remove_cv<typename std::remove_reference<
                  StaticIterable>::type>::type::iterator>::value_type,
          types::len_of<typename std::remove_cv<typename std::remove_reference<
              StaticIterable>::type>::type>::value>>::type
  tuple(StaticIterable &&i)
  {
    types::array<
        typename std::iterator_traits<
            typename std::remove_cv<typename std::remove_reference<
                StaticIterable>::type>::type::iterator>::value_type,
        types::len_of<typename std::remove_cv<typename std::remove_reference<
            StaticIterable>::type>::type>::value> res;
    std::copy(i.begin(), i.end(), res.begin());
    return res;
  }
}
PYTHONIC_NS_END

#endif

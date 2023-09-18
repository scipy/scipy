#ifndef PYTHONIC_INCLUDE_BUILTIN_MINMAX_HPP
#define PYTHONIC_INCLUDE_BUILTIN_MINMAX_HPP

#include <utility>
#include "pythonic/include/builtins/pythran/kwonly.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace details
  {
    template <class Op, class T>
    typename std::decay<T>::type::value_type minmax(Op const &, T &&t);

    template <class Op, class T, class F>
    typename std::decay<T>::type::value_type minmax(Op const &, T &&t,
                                                    types::kwonly, F key);

    template <class Op, class T0, class T1, class... Types>
    typename std::enable_if<!std::is_same<T1, types::kwonly>::value,
                            typename __combined<T0, T1, Types...>::type>::type
    minmax(Op const &, T0 const &, T1 const &, Types const &...);
  }
}
PYTHONIC_NS_END

#endif

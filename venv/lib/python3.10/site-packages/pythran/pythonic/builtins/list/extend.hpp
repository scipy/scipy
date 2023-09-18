#ifndef PYTHONIC_BUILTIN_LIST_EXTEND_HPP
#define PYTHONIC_BUILTIN_LIST_EXTEND_HPP

#include "pythonic/include/builtins/list/extend.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/types/NoneType.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {
    template <class T0, class T1>
    typename std::enable_if<
        !std::is_same<typename std::decay<T0>::type, types::empty_list>::value,
        types::none_type>::type
    extend(T0 &&seq, T1 const &add)
    {
      std::forward<T0>(seq) += add;
      return {};
    }

    template <class T0, class T1>
    typename std::enable_if<
        std::is_same<typename std::decay<T0>::type, types::empty_list>::value,
        types::none_type>::type
    extend(T0 &&seq, T1 const &add)
    {
      return {};
    }
  }
}
PYTHONIC_NS_END
#endif

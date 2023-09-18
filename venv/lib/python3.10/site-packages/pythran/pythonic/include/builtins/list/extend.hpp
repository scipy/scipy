#ifndef PYTHONIC_INCLUDE_BUILTIN_LIST_EXTEND_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LIST_EXTEND_HPP

#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {

    template <class T0, class T1>
    typename std::enable_if<
        !std::is_same<typename std::decay<T0>::type, types::empty_list>::value,
        types::none_type>::type
    extend(T0 &&seq, T1 const &add);

    template <class T0, class T1>
    typename std::enable_if<
        std::is_same<typename std::decay<T0>::type, types::empty_list>::value,
        types::none_type>::type
    extend(T0 &&seq, T1 const &add);

    DEFINE_FUNCTOR(pythonic::builtins::list, extend);
  }
}
PYTHONIC_NS_END
#endif

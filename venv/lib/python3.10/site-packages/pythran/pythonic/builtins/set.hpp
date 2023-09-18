#ifndef PYTHONIC_BUILTIN_SET_HPP
#define PYTHONIC_BUILTIN_SET_HPP

#include "pythonic/include/builtins/set.hpp"

#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {

    inline types::empty_set set()
    {
      return types::empty_set();
    }

    template <class Iterable>
    inline types::set<typename std::iterator_traits<
        typename std::remove_reference<Iterable>::type::iterator>::value_type>
    set(Iterable &&t)
    {
      return {t.begin(), t.end()};
    }
  }
}
PYTHONIC_NS_END
#endif

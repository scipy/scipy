#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_HPP

#include "pythonic/include/types/set.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    inline types::empty_set set();

    template <class Iterable>
    inline types::set<typename std::iterator_traits<
        typename std::remove_reference<Iterable>::type::iterator>::value_type>
    set(Iterable &&t);
  }

  DEFINE_FUNCTOR(pythonic::builtins::anonymous, set);
}
PYTHONIC_NS_END
#endif

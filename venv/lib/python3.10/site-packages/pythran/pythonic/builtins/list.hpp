#ifndef PYTHONIC_BUILTIN_LIST_HPP
#define PYTHONIC_BUILTIN_LIST_HPP

#include "pythonic/include/builtins/list.hpp"

#include "pythonic/types/list.hpp"
#include "pythonic/utils/functor.hpp"

#include <iterator>
#include <type_traits>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {

    inline types::empty_list list()
    {
      return types::empty_list();
    }

    inline types::empty_list list(types::empty_list)
    {
      return types::empty_list();
    }

    template <class Iterable>
    types::list<typename std::decay<typename std::iterator_traits<
        typename std::remove_reference<Iterable>::type::iterator>::value_type>::
                    type>
    list(Iterable &&t)
    {
      return types::list<typename std::decay<
          typename std::iterator_traits<typename std::remove_reference<
              Iterable>::type::iterator>::value_type>::type>(t.begin(),
                                                             t.end());
    }
  }
}
PYTHONIC_NS_END

#endif

#ifndef PYTHONIC_INCLUDE_BUILTIN_ITER_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ITER_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {
    template <class T>
    struct iter : T::iterator {
      using iterator = typename T::iterator;

      iterator _end;
      T data;

      iter();
      iter(T data);
      iterator &begin();
      iterator const &begin() const;
      iterator const &end() const;
    };
  }

  template <class T>
  details::iter<
      typename std::remove_cv<typename std::remove_reference<T>::type>::type>
  iter(T &&t);

  DEFINE_FUNCTOR(pythonic::builtins, iter);
}
PYTHONIC_NS_END

#endif

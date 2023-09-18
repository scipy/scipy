#ifndef PYTHONIC_INCLUDE_BUILTIN_TYPE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_TYPE_HPP

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  struct type_functor;

  template <class T>
  typename type_functor<T>::type type(T const &t);
  DEFINE_FUNCTOR(pythonic::builtins, type);
}
PYTHONIC_NS_END

#endif

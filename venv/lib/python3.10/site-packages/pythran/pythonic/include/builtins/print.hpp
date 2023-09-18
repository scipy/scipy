#ifndef PYTHONIC_INCLUDE_BUILTIN_PRINT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PRINT_HPP

#include <ostream>
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  void print_nonl();

  template <typename T, typename... Types>
  void print_nonl(T const &value, Types const &... values);

  void print();

  template <typename T, typename... Types>
  void print(T const &value, Types const &... values);
  DEFINE_FUNCTOR(pythonic::builtins, print);
}
PYTHONIC_NS_END

#endif

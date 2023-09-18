#ifndef PYTHONIC_INCLUDE_BUILTIN_BOOL_HPP
#define PYTHONIC_INCLUDE_BUILTIN_BOOL_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/tuple.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {

    struct bool_ {
      using callable = void;
      using type = bool;

      bool operator()() const;

      template <class T>
      bool operator()(T const &val) const;

      template <class... Ts>
      bool operator()(std::tuple<Ts...> const &val) const;

      template <class T, size_t N>
      bool operator()(types::array<T, N> const &val) const;

      friend std::ostream &operator<<(std::ostream &os, bool_)
      {
        return os << "bool";
      }
    };
  }
}
PYTHONIC_NS_END

#endif

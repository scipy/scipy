#ifndef PYTHONIC_INCLUDE_BUILTIN_COMPLEX_HPP
#define PYTHONIC_INCLUDE_BUILTIN_COMPLEX_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/complex.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {
    struct complex {
      using callable = void;
      using type = std::complex<double>;
      // TODO: doesn't handle string as first argument
      type operator()(double v0 = 0, double v1 = 0) const;
      friend std::ostream &operator<<(std::ostream &os, complex)
      {
        return os << "complex";
      }
    };
  }
}
PYTHONIC_NS_END

#endif

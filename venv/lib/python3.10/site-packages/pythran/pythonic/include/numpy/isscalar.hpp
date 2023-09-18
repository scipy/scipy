#ifndef PYTHONIC_INCLUDE_NUMPY_ISSCALAR_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISSCALAR_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/traits.hpp"
#include "pythonic/include/types/str.hpp"

#include <type_traits>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E>
  constexpr bool isscalar(E const &);

  DEFINE_FUNCTOR(pythonic::numpy, isscalar);
}
PYTHONIC_NS_END

#endif

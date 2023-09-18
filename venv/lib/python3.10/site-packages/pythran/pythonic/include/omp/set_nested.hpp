#ifndef PYTHONIC_INCLUDE_OMP_SET_NESTED_HPP
#define PYTHONIC_INCLUDE_OMP_SET_NESTED_HPP

#include <omp.h>
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{

  void set_nested(long val);

  DEFINE_FUNCTOR(pythonic::omp, set_nested);
}
PYTHONIC_NS_END

#endif

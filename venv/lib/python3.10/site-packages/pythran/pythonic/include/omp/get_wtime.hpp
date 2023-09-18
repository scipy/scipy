#ifndef PYTHONIC_INCLUDE_OMP_GET_WTIME_HPP
#define PYTHONIC_INCLUDE_OMP_GET_WTIME_HPP

#include <omp.h>
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{

  long get_wtime();

  DEFINE_FUNCTOR(pythonic::omp, get_wtime);
}
PYTHONIC_NS_END

#endif

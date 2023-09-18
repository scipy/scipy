#ifndef PYTHONIC_OMP_GET_WTIME_HPP
#define PYTHONIC_OMP_GET_WTIME_HPP

#include "pythonic/include/omp/get_wtime.hpp"

#include <omp.h>
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{

  long get_wtime()
  {
    return omp_get_wtime();
  }
}
PYTHONIC_NS_END

#endif

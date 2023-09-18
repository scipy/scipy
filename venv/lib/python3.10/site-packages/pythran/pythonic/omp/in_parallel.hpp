#ifndef PYTHONIC_OMP_IN_PARALLEL_HPP
#define PYTHONIC_OMP_IN_PARALLEL_HPP

#include "pythonic/include/omp/in_parallel.hpp"

#include <omp.h>
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{

  bool in_parallel()
  {
    return omp_in_parallel();
  }
}
PYTHONIC_NS_END

#endif

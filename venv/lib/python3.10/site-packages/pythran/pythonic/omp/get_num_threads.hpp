#ifndef PYTHONIC_OMP_GET_NUM_THREADS_HPP
#define PYTHONIC_OMP_GET_NUM_THREADS_HPP

#include "pythonic/include/omp/get_num_threads.hpp"

#include <omp.h>
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{
  long get_num_threads()
  {
    return omp_get_num_threads();
  }
}
PYTHONIC_NS_END

#endif

#ifndef PYTHONIC_OMP_GET_THREAD_NUM_HPP
#define PYTHONIC_OMP_GET_THREAD_NUM_HPP

#include "pythonic/include/omp/get_thread_num.hpp"

#include <omp.h>
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{

  long get_thread_num()
  {
    return omp_get_thread_num();
  }
}
PYTHONIC_NS_END

#endif

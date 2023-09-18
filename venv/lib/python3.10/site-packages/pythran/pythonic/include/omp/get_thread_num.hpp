#ifndef PYTHONIC_INCLUDE_OMP_GET_THREAD_NUM_HPP
#define PYTHONIC_INCLUDE_OMP_GET_THREAD_NUM_HPP

#include <omp.h>
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{

  long get_thread_num();

  DEFINE_FUNCTOR(pythonic::omp, get_thread_num);
}
PYTHONIC_NS_END

#endif

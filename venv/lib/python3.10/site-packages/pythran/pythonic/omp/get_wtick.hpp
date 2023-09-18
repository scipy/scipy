#ifndef PYTHONIC_OMP_GET_WTICK_HPP
#define PYTHONIC_OMP_GET_WTICK_HPP

#include "pythonic/include/omp/get_wtick.hpp"

#include <omp.h>
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{
  long get_wtick()
  {
    return omp_get_wtick();
  }
}
PYTHONIC_NS_END

#endif

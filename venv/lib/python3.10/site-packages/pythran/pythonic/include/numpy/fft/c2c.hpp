#ifndef PYTHONIC_INCLUDE_NUMPY_FFT_C2C_HPP
#define PYTHONIC_INCLUDE_NUMPY_FFT_C2C_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace fft
  {

    template <class T, class pS>
    types::ndarray<std::complex<T>,
                   types::array<long, std::tuple_size<pS>::value>>
    c2c(types::ndarray<std::complex<T>, pS> const &a, long n = -1,
        long axis = -1, types::str const &norm = {}, bool const forward = true);
  }
}
PYTHONIC_NS_END

#endif

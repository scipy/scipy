#ifndef PYTHONIC_NUMPY_FFT_IRFFT_HPP
#define PYTHONIC_NUMPY_FFT_IRFFT_HPP

#include "pythonic/include/numpy/fft/irfft.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/include/utils/array_helper.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/builtins/None.hpp"
#include "pythonic/numpy/fft/c2c.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace fft
  {

    template <class T, class pS>
    types::ndarray<T, types::array<long, std::tuple_size<pS>::value>>
    irfft(types::ndarray<std::complex<T>, pS> const &in_array,
          types::none_type n, long axis, types::str const &norm)
    {
      return c2r(in_array, -1, axis, norm, false);
    }

    template <class T, class pS>
    types::ndarray<T, types::array<long, std::tuple_size<pS>::value>>
    irfft(types::ndarray<std::complex<T>, pS> const &in_array,
          types::none_type n, long axis, types::none_type norm)
    {
      return c2r(in_array, -1, axis, "", false);
    }

    template <class T, class pS>
    types::ndarray<T, types::array<long, std::tuple_size<pS>::value>>
    irfft(types::ndarray<std::complex<T>, pS> const &in_array, long n,
          long axis, types::none_type norm)
    {
      return c2r(in_array, n, axis, "", false);
    }

    template <class T, class pS>
    types::ndarray<T, types::array<long, std::tuple_size<pS>::value>>
    irfft(types::ndarray<std::complex<T>, pS> const &in_array, long n,
          long axis, types::str const &norm)
    {
      return c2r(in_array, n, axis, norm, false);
    }

    template <class T, class pS>
    types::ndarray<typename std::enable_if<
                       !types::is_complex<T>::value,
                       typename std::conditional<std::is_integral<T>::value,
                                                 double, T>::type>::type,
                   types::array<long, std::tuple_size<pS>::value>>
    irfft(types::ndarray<T, pS> const &in_array, types::none_type n, long axis,
          types::str const &norm)
    {
      auto tmp_array = _copy_to_complex(in_array);
      return c2r(tmp_array, -1, axis, norm, false);
    }

    template <class T, class pS>
    types::ndarray<typename std::enable_if<
                       !types::is_complex<T>::value,
                       typename std::conditional<std::is_integral<T>::value,
                                                 double, T>::type>::type,
                   types::array<long, std::tuple_size<pS>::value>>
    irfft(types::ndarray<T, pS> const &in_array, types::none_type n, long axis,
          types::none_type norm)
    {
      auto tmp_array = _copy_to_complex(in_array);
      return c2r(tmp_array, -1, axis, "", false);
    }

    template <class T, class pS>
    types::ndarray<typename std::enable_if<
                       !types::is_complex<T>::value,
                       typename std::conditional<std::is_integral<T>::value,
                                                 double, T>::type>::type,
                   types::array<long, std::tuple_size<pS>::value>>
    irfft(types::ndarray<T, pS> const &in_array, long n, long axis,
          types::none_type norm)
    {
      auto tmp_array = _copy_to_complex(in_array);
      return c2r(tmp_array, n, axis, "", false);
    }

    template <class T, class pS>
    types::ndarray<typename std::enable_if<
                       !types::is_complex<T>::value,
                       typename std::conditional<std::is_integral<T>::value,
                                                 double, T>::type>::type,
                   types::array<long, std::tuple_size<pS>::value>>
    irfft(types::ndarray<T, pS> const &in_array, long n, long axis,
          types::str const &norm)
    {
      auto tmp_array = _copy_to_complex(in_array);
      return c2r(tmp_array, n, axis, norm, false);
    }

    NUMPY_EXPR_TO_NDARRAY0_IMPL(irfft);
  }
}
PYTHONIC_NS_END

#endif

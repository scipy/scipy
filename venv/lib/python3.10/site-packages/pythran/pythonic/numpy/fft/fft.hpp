#ifndef PYTHONIC_NUMPY_FFT_FFT_HPP
#define PYTHONIC_NUMPY_FFT_FFT_HPP

#include "pythonic/include/numpy/fft/fft.hpp"
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
    namespace details
    {
      inline types::str normalize_norm(types::none_type const &)
      {
        return "backward";
      }

      template <class T>
      inline types::str normalize_norm(T const &norm)
      {
        return norm;
      }

      inline long normalize_n(types::none_type const &)
      {
        return -1;
      }

      inline long normalize_n(long const &n)
      {
        return n;
      }
    }

    template <class T, class pS, class N, class Norm>
    types::ndarray<
        typename std::enable_if<types::is_complex<T>::value, T>::type,
        types::array<long, std::tuple_size<pS>::value>>
    fft(types::ndarray<T, pS> const &in_array, N const& n, long axis,
        Norm const &norm)
    {
      return c2c(in_array, details::normalize_n(n), axis, details::normalize_norm(norm), true);
    }

    template <class T, class pS, class N, class Norm>
    types::ndarray<typename std::enable_if<std::is_floating_point<T>::value,
                                           std::complex<T>>::type,
                   types::array<long, std::tuple_size<pS>::value>>
    fft(types::ndarray<T, pS> const &in_array, N const & n, long axis,
        Norm const &norm)
    {
      return r2c(in_array, details::normalize_n(n), axis, details::normalize_norm(norm), true, true);
    }

    template <class T, class pS, class N, class Norm>
    types::ndarray<typename std::enable_if<std::is_integral<T>::value,
                                           std::complex<double>>::type,
                   types::array<long, std::tuple_size<pS>::value>>
    fft(types::ndarray<T, pS> const &in_array, N const& n, long axis,
        Norm const &norm)
    {
      auto tmp_array = _copy_to_double(in_array);
      return fft(tmp_array, n, axis, norm);
    }

    NUMPY_EXPR_TO_NDARRAY0_IMPL(fft);
  }
}
PYTHONIC_NS_END

#endif

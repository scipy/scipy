/** SciPy-only enhancements to pypocketfft

 - Optimized full-sized fft for real input
 - Normalization factor calculation

*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <complex>

#include "pocketfft.h"

// Only instantiate long double transforms if they offer more precision
using longfloat_t = typename std::conditional<
  sizeof(long double) == sizeof(double), double, long double>::type;

template<typename T> struct array_iterator
  {
  T * data, * cur;

  pocketfft::stride_t strides;
  pocketfft::shape_t shape, idx;

  array_iterator(T * data_, pocketfft::stride_t strides_,
                 pocketfft::shape_t shape_):
    data(data_),
    cur(data_),
    strides(strides_),
    shape(shape_),
    idx(shape.size())
    {
    }

  void reset()
    {
    std::fill(idx.begin(), idx.end(), 0);
    cur = data;
    }

  /** Advance by 1, wrapping around into higher dimensions as necessary */
  bool advance()
    {
    for (ptrdiff_t i = shape.size() - 1; i >= 0; --i)
      {
      ++idx[i];
      cur += strides[i];

      if (idx[i] < shape[i])
        return true;

      cur -= idx[i] * strides[i];
      idx[i] = 0;
      }
    return false;
    }

  bool empty()
    {
    return std::any_of(shape.begin(), shape.end(),
                       [](size_t len){ return len == 0; });
    }

  /** Take a slice along a single axis and reset index to (0,0,...,0) */
  void slice(size_t dim, ptrdiff_t begin, ptrdiff_t end, ptrdiff_t step=1)
    {
    data += begin * strides[dim];
    shape[dim] = end - begin;
    strides[dim] *= step;
    reset();
    }
  };

enum class norm_t
  {
  none,  // No normalization
  ortho, // Orthonormal transform
  size   // Normalize by the fft size
  };

template<typename T> T norm_fct(norm_t norm, const pocketfft::shape_t & shape,
                                const pocketfft::shape_t & axes)
  {
  if (norm == norm_t::none) return T(1);

  size_t N = 1;
  for (auto a : axes)
    N *= shape[a];

  auto nl = static_cast<long double>(N);
  if (norm == norm_t::size) return static_cast<T>(1.L / nl);
  if (norm == norm_t::ortho) return static_cast<T>(1.L / std::sqrt(nl));

  throw std::invalid_argument("invalid normalization type");
  }


pybind11::array asfarray(const pybind11::array & a, bool & inplace)
  {
    auto dtype = a.dtype();
    switch (dtype.kind())
      {
      case 'c': // complex
        return a;

      case 'f': // floating
        {
        static const auto f16 = pybind11::dtype("float16");
        if (!dtype.is(f16))
          return a;

        inplace = true;
        return a.cast<pybind11::array_t<float>>();
        }

      default:
        {
        inplace = true;
        return a.cast<pybind11::array_t<double>>();
        }
      }
  }

pybind11::array asfarray(const pybind11::array & a)
  {
  bool ignored;
  return asfarray(a, ignored);
  }


namespace scipy {
// Copied from pypocketfft.cxx
pocketfft::shape_t copy_shape(const pybind11::array &arr)
{
  pocketfft::shape_t res(size_t(arr.ndim()));
  for (size_t i=0; i<res.size(); ++i)
    res[i] = size_t(arr.shape(int(i)));
  return res;
}

// Copied from pypocketfft.cxx
pocketfft::stride_t copy_strides(const pybind11::array &arr)
{
  pocketfft::stride_t res(size_t(arr.ndim()));
  for (size_t i=0; i<res.size(); ++i)
    res[i] = arr.strides(int(i));
  return res;
}
}

template<typename T> pybind11::array xfftn_real(const pybind11::array &in,
  const pocketfft::shape_t & axes, norm_t norm, bool /*inplace*/, bool fwd,
  size_t nthreads)
  {
  auto dims_in = scipy::copy_shape(in);
  auto dims_out = dims_in;

  const auto full_len = dims_in[axes.back()];
  const auto rfft_len = full_len / 2 + 1;

  dims_out[axes.back()] = rfft_len;
  pybind11::array_t<std::complex<T>> res(dims_in);

  // Perform real fft on the last axis
  auto strides_in = scipy::copy_strides(in);
  auto strides_out = scipy::copy_strides(res);
  auto * in_data = reinterpret_cast<const T *>(in.data());
  auto * data = res.mutable_data();
  const auto size = static_cast<size_t>(res.size());

  {
  pybind11::gil_scoped_release release;
  auto fct = norm_fct<T>(norm, dims_in, axes);

  pocketfft::r2c(dims_in, strides_in, strides_out, axes.back(), in_data, data,
                 fct, nthreads);

  // For real input, ifft is the conjugate of the fft
  if (!fwd)
    {
    for (size_t i = 0; i < size; i++)
      data[i] = conj(data[i]);
    }

  // Perform complex (i)fft on remaining axes
  if (axes.size() > 1)
    {
    pocketfft::shape_t new_axes(axes.begin(), axes.end()-1);
    pocketfft::c2c(dims_out, strides_out, strides_out, new_axes, fwd, data, data,
                   T(1), nthreads);
    }

  // Convert strides from bytes to objects
  for (auto & s : strides_out)
    s /= sizeof(std::complex<T>);

  // Use hermitian symmetry to fill remaining output
  array_iterator<std::complex<T>> work(data, strides_out, dims_in);
  work.slice(axes.back(), 1, (full_len + 1) / 2);

  if (!work.empty())
    {
    do {
      auto mirror_offset =
        dims_in[axes.back()] - 2 - 2*work.idx[axes.back()];
      work.cur[mirror_offset * work.strides[axes.back()]] =
        conj(work.cur[0]);
      } while (work.advance());
    }


  // FFT(conj(x)) is the reverse of conj(FFT(x)) so for each complex transform,
  // we reverse the output along the FFT axis
  for (auto itr = axes.begin(); itr != prev(axes.end()); ++itr)
    {
    auto m_axis = *itr;

    array_iterator<std::complex<T>> work(data, strides_out, dims_in);
    work.slice(axes.back(), rfft_len, full_len);
    work.slice(m_axis, 1, (dims_in[m_axis]+1) / 2);

    if (work.empty())
      continue;

    do {
      auto mirror_offset = dims_in[m_axis] - 2 - 2*work.idx[m_axis];
      swap(work.cur[0], work.cur[mirror_offset * work.strides[m_axis]]);
      } while(work.advance());
    }

  } // Aquire GIL

  return res;
  }

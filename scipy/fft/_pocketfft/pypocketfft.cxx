/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*
 *  Python interface.
 *
 *  Copyright (C) 2019 Max-Planck-Society
 *  Copyright (C) 2019 Peter Bell
 *  \author Martin Reinecke
 *  \author Peter Bell
 */

// Note: do not reorder, pybind11.h must come first (it includes Python.h)
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "pocketfft_hdronly.h"

namespace {

using pocketfft::shape_t;
using pocketfft::stride_t;
using std::size_t;
using std::ptrdiff_t;

namespace py = pybind11;

// Only instantiate long double transforms if they offer more precision
using ldbl_t = typename std::conditional<
  sizeof(long double)==sizeof(double), double, long double>::type;

using c64 = std::complex<float>;
using c128 = std::complex<double>;
using clong = std::complex<ldbl_t>;
using f32 = float;
using f64 = double;
using flong = ldbl_t;

shape_t copy_shape(const py::array &arr)
  {
  shape_t res(size_t(arr.ndim()));
  for (size_t i=0; i<res.size(); ++i)
    res[i] = size_t(arr.shape(int(i)));
  return res;
  }

stride_t copy_strides(const py::array &arr)
  {
  stride_t res(size_t(arr.ndim()));
  for (size_t i=0; i<res.size(); ++i)
    res[i] = arr.strides(int(i));
  return res;
  }

shape_t makeaxes(const py::array &in, const py::object &axes)
  {
  if (axes.is_none())
    {
    shape_t res(size_t(in.ndim()));
    for (size_t i=0; i<res.size(); ++i)
      res[i]=i;
    return res;
    }
  auto tmp=axes.cast<std::vector<ptrdiff_t>>();
  auto ndim = in.ndim();
  if ((tmp.size()>size_t(ndim)) || (tmp.size()==0))
    throw std::runtime_error("bad axes argument");
  for (auto& sz: tmp)
    {
    if (sz<0)
      sz += ndim;
    if ((sz>=ndim) || (sz<0))
      throw std::invalid_argument("axes exceeds dimensionality of output");
    }
  return shape_t(tmp.begin(), tmp.end());
  }

#define DISPATCH(arr, T1, T2, T3, func, args) \
  { \
  if (py::isinstance<py::array_t<T1>>(arr)) return func<double> args; \
  if (py::isinstance<py::array_t<T2>>(arr)) return func<float> args;  \
  if (py::isinstance<py::array_t<T3>>(arr)) return func<ldbl_t> args; \
  throw std::runtime_error("unsupported data type"); \
  }

template<typename T> T norm_fct(int inorm, size_t N)
  {
  if (inorm==0) return T(1);
  if (inorm==2) return T(1/ldbl_t(N));
  if (inorm==1) return T(1/sqrt(ldbl_t(N)));
  throw std::invalid_argument("invalid value for inorm (must be 0, 1, or 2)");
  }

template<typename T> T norm_fct(int inorm, const shape_t &shape,
  const shape_t &axes, size_t fct=1, int delta=0)
  {
  if (inorm==0) return T(1);
  size_t N(1);
  for (auto a: axes)
    N *= fct * size_t(int64_t(shape[a])+delta);
  return norm_fct<T>(inorm, N);
  }

template<typename T> py::array_t<T> prepare_output(py::object &out_,
  shape_t &dims)
  {
  if (out_.is_none()) return py::array_t<T>(dims);
  auto tmp = out_.cast<py::array_t<T>>();
  if (!tmp.is(out_)) // a new object was created during casting
    throw std::runtime_error("unexpected data type for output array");
  return tmp;
  }

template<typename T> py::array c2c_internal(const py::array &in,
  const py::object &axes_, bool forward, int inorm, py::object &out_,
  size_t nthreads)
  {
  auto axes = makeaxes(in, axes_);
  auto dims(copy_shape(in));
  auto res = prepare_output<std::complex<T>>(out_, dims);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const std::complex<T> *>(in.data());
  auto d_out=reinterpret_cast<std::complex<T> *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, dims, axes);
  pocketfft::c2c(dims, s_in, s_out, axes, forward, d_in, d_out, fct, nthreads);
  }
  return std::move(res);
  }

template<typename T> py::array c2c_sym_internal(const py::array &in,
  const py::object &axes_, bool forward, int inorm, py::object &out_,
  size_t nthreads)
  {
  auto axes = makeaxes(in, axes_);
  auto dims(copy_shape(in));
  auto res = prepare_output<std::complex<T>>(out_, dims);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<std::complex<T> *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, dims, axes);
  pocketfft::r2c(dims, s_in, s_out, axes, forward, d_in, d_out, fct, nthreads);
  // now fill in second half
  using namespace pocketfft::detail;
  ndarr<std::complex<T>> ares(res.mutable_data(), dims, s_out);
  rev_iter iter(ares, axes);
  while(iter.remaining()>0)
    {
    auto v = ares[iter.ofs()];
    ares[iter.rev_ofs()] = conj(v);
    iter.advance();
    }
  }
  return std::move(res);
  }

py::array c2c(const py::array &a, const py::object &axes_, bool forward,
  int inorm, py::object &out_, size_t nthreads)
  {
  if (a.dtype().kind() == 'c')
    DISPATCH(a, c128, c64, clong, c2c_internal, (a, axes_, forward,
             inorm, out_, nthreads))

  DISPATCH(a, f64, f32, flong, c2c_sym_internal, (a, axes_, forward,
           inorm, out_, nthreads))
  }

template<typename T> py::array r2c_internal(const py::array &in,
  const py::object &axes_, bool forward, int inorm, py::object &out_,
  size_t nthreads)
  {
  auto axes = makeaxes(in, axes_);
  auto dims_in(copy_shape(in)), dims_out(dims_in);
  dims_out[axes.back()] = (dims_out[axes.back()]>>1)+1;
  py::array res = prepare_output<std::complex<T>>(out_, dims_out);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<std::complex<T> *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, dims_in, axes);
  pocketfft::r2c(dims_in, s_in, s_out, axes, forward, d_in, d_out, fct,
    nthreads);
  }
  return res;
  }

py::array r2c(const py::array &in, const py::object &axes_, bool forward,
  int inorm, py::object &out_, size_t nthreads)
  {
  DISPATCH(in, f64, f32, flong, r2c_internal, (in, axes_, forward, inorm, out_,
    nthreads))
  }

template<typename T> py::array r2r_fftpack_internal(const py::array &in,
  const py::object &axes_, bool real2hermitian, bool forward, int inorm,
  py::object &out_, size_t nthreads)
  {
  auto axes = makeaxes(in, axes_);
  auto dims(copy_shape(in));
  py::array res = prepare_output<T>(out_, dims);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<T *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, dims, axes);
  pocketfft::r2r_fftpack(dims, s_in, s_out, axes, real2hermitian, forward,
    d_in, d_out, fct, nthreads);
  }
  return res;
  }

py::array r2r_fftpack(const py::array &in, const py::object &axes_,
  bool real2hermitian, bool forward, int inorm, py::object &out_,
  size_t nthreads)
  {
  DISPATCH(in, f64, f32, flong, r2r_fftpack_internal, (in, axes_,
    real2hermitian, forward, inorm, out_, nthreads))
  }

template<typename T> py::array dct_internal(const py::array &in,
  const py::object &axes_, int type, int inorm, py::object &out_,
  size_t nthreads, bool ortho)
  {
  auto axes = makeaxes(in, axes_);
  auto dims(copy_shape(in));
  py::array res = prepare_output<T>(out_, dims);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<T *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  T fct = (type==1) ? norm_fct<T>(inorm, dims, axes, 2, -1)
                    : norm_fct<T>(inorm, dims, axes, 2);
  pocketfft::dct(dims, s_in, s_out, axes, type, d_in, d_out, fct, ortho,
    nthreads);
  }
  return res;
  }

py::array dct(const py::array &in, int type, const py::object &axes_,
  int inorm, py::object &out_, size_t nthreads, const py::object & ortho_obj)
  {
  bool ortho=inorm==1;
  if (!ortho_obj.is_none())
    ortho=ortho_obj.cast<bool>();

  if ((type<1) || (type>4)) throw std::invalid_argument("invalid DCT type");
  DISPATCH(in, f64, f32, flong, dct_internal, (in, axes_, type, inorm, out_,
    nthreads, ortho))
  }

template<typename T> py::array dst_internal(const py::array &in,
  const py::object &axes_, int type, int inorm, py::object &out_,
  size_t nthreads, bool ortho)
  {
  auto axes = makeaxes(in, axes_);
  auto dims(copy_shape(in));
  py::array res = prepare_output<T>(out_, dims);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<T *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  T fct = (type==1) ? norm_fct<T>(inorm, dims, axes, 2, 1)
                    : norm_fct<T>(inorm, dims, axes, 2);
  pocketfft::dst(dims, s_in, s_out, axes, type, d_in, d_out, fct, ortho,
    nthreads);
  }
  return res;
  }

py::array dst(const py::array &in, int type, const py::object &axes_,
  int inorm, py::object &out_, size_t nthreads, const py::object &ortho_obj)
  {
  bool ortho=inorm==1;
  if (!ortho_obj.is_none())
    ortho=ortho_obj.cast<bool>();

  if ((type<1) || (type>4)) throw std::invalid_argument("invalid DST type");
  DISPATCH(in, f64, f32, flong, dst_internal, (in, axes_, type, inorm,
    out_, nthreads, ortho))
  }

template<typename T> py::array c2r_internal(const py::array &in,
  const py::object &axes_, size_t lastsize, bool forward, int inorm,
  py::object &out_, size_t nthreads)
  {
  auto axes = makeaxes(in, axes_);
  size_t axis = axes.back();
  shape_t dims_in(copy_shape(in)), dims_out=dims_in;
  if (lastsize==0) lastsize=2*dims_in[axis]-1;
  if ((lastsize/2) + 1 != dims_in[axis])
    throw std::invalid_argument("bad lastsize");
  dims_out[axis] = lastsize;
  py::array res = prepare_output<T>(out_, dims_out);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const std::complex<T> *>(in.data());
  auto d_out=reinterpret_cast<T *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, dims_out, axes);
  pocketfft::c2r(dims_out, s_in, s_out, axes, forward, d_in, d_out, fct,
    nthreads);
  }
  return res;
  }

py::array c2r(const py::array &in, const py::object &axes_, size_t lastsize,
  bool forward, int inorm, py::object &out_, size_t nthreads)
  {
  DISPATCH(in, c128, c64, clong, c2r_internal, (in, axes_, lastsize, forward,
    inorm, out_, nthreads))
  }

template<typename T> py::array separable_hartley_internal(const py::array &in,
  const py::object &axes_, int inorm, py::object &out_, size_t nthreads)
  {
  auto dims(copy_shape(in));
  py::array res = prepare_output<T>(out_, dims);
  auto axes = makeaxes(in, axes_);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<T *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, dims, axes);
  pocketfft::r2r_separable_hartley(dims, s_in, s_out, axes, d_in, d_out, fct,
    nthreads);
  }
  return res;
  }

py::array separable_hartley(const py::array &in, const py::object &axes_,
  int inorm, py::object &out_, size_t nthreads)
  {
  DISPATCH(in, f64, f32, flong, separable_hartley_internal, (in, axes_, inorm,
    out_, nthreads))
  }

template<typename T> py::array genuine_hartley_internal(const py::array &in,
  const py::object &axes_, int inorm, py::object &out_, size_t nthreads)
  {
  auto dims(copy_shape(in));
  py::array res = prepare_output<T>(out_, dims);
  auto axes = makeaxes(in, axes_);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<T *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, dims, axes);
  pocketfft::r2r_genuine_hartley(dims, s_in, s_out, axes, d_in, d_out, fct,
    nthreads);
  }
  return res;
  }

py::array genuine_hartley(const py::array &in, const py::object &axes_,
  int inorm, py::object &out_, size_t nthreads)
  {
  DISPATCH(in, f64, f32, flong, genuine_hartley_internal, (in, axes_, inorm,
    out_, nthreads))
  }

// Export good_size in raw C-API to reduce overhead (~4x faster)
PyObject * good_size(PyObject * /*self*/, PyObject * args, PyObject * kwargs)
  {
  Py_ssize_t n_ = -1;
  int real = false;
  static const char * keywords[] = {"target", "real", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "n|p:good_size",
                                   (char **) keywords, &n_, &real))
    return nullptr;

  if (n_<0)
    {
    PyErr_SetString(PyExc_ValueError, "Target length must be positive");
    return nullptr;
    }
  if ((n_-1) > static_cast<Py_ssize_t>(std::numeric_limits<size_t>::max() / 11))
    {
    PyErr_Format(PyExc_ValueError,
                 "Target length is too large to perform an FFT: %zi", n_);
    return nullptr;
    }
  const auto n = static_cast<size_t>(n_);
  using namespace pocketfft::detail;
  return PyLong_FromSize_t(
    real ? util::good_size_real(n) : util::good_size_cmplx(n));
  }

// Export prev_good_size in raw C-API to reduce overhead
PyObject * prev_good_size(PyObject * /*self*/, PyObject * args, PyObject * kwargs)
  {
  Py_ssize_t n_ = -1;
  int real = false;
  static const char * keywords[] = {"target", "real", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "n|p:prev_good_size",
                                   (char **) keywords, &n_, &real))
    return nullptr;

  if (n_<0)
    {
    PyErr_SetString(PyExc_ValueError, "Target length must be positive");
    return nullptr;
    }
  if ((n_-1) > static_cast<Py_ssize_t>(std::numeric_limits<size_t>::max() / 11))
    {
    PyErr_Format(PyExc_ValueError,
                 "Target length is too large to perform an FFT: %zi", n_);
    return nullptr;
    }
  const auto n = static_cast<size_t>(n_);
  using namespace pocketfft::detail;
  return PyLong_FromSize_t(
    real ? util::prev_good_size_real(n) : util::prev_good_size_cmplx(n));
  }

const char *pypocketfft_DS = R"""(Fast Fourier and Hartley transforms.

This module supports
- single, double, and long double precision
- complex and real-valued transforms
- multi-dimensional transforms

For two- and higher-dimensional transforms the code will use SSE2 and AVX
vector instructions for faster execution if these are supported by the CPU and
were enabled during compilation.
)""";

const char *c2c_DS = R"""(Performs a complex FFT.

Parameters
----------
a : numpy.ndarray (any complex or real type)
    The input data. If its type is real, a more efficient real-to-complex
    transform will be used.
axes : list of integers
    The axes along which the FFT is carried out.
    If not set, all axes will be transformed.
forward : bool
    If `True`, a negative sign is used in the exponent, else a positive one.
inorm : int
    Normalization type
      0 : no normalization
      1 : divide by sqrt(N)
      2 : divide by N
    where N is the product of the lengths of the transformed axes.
out : numpy.ndarray (same shape as `a`, complex type with same accuracy as `a`)
    May be identical to `a`, but if it isn't, it must not overlap with `a`.
    If None, a new array is allocated to store the output.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
numpy.ndarray (same shape as `a`, complex type with same accuracy as `a`)
    The transformed data.
)""";

const char *r2c_DS = R"""(Performs an FFT whose input is strictly real.

Parameters
----------
a : numpy.ndarray (any real type)
    The input data
axes : list of integers
    The axes along which the FFT is carried out.
    If not set, all axes will be transformed in ascending order.
forward : bool
    If `True`, a negative sign is used in the exponent, else a positive one.
inorm : int
    Normalization type
      0 : no normalization
      1 : divide by sqrt(N)
      2 : divide by N
    where N is the product of the lengths of the transformed input axes.
out : numpy.ndarray (complex type with same accuracy as `a`)
    For the required shape, see the `Returns` section.
    Must not overlap with `a`.
    If None, a new array is allocated to store the output.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
numpy.ndarray (complex type with same accuracy as `a`)
    The transformed data. The shape is identical to that of the input array,
    except for the axis that was transformed last. If the length of that axis
    was n on input, it is n//2+1 on output.
)""";

const char *c2r_DS = R"""(Performs an FFT whose output is strictly real.

Parameters
----------
a : numpy.ndarray (any complex type)
    The input data
axes : list of integers
    The axes along which the FFT is carried out.
    If not set, all axes will be transformed in ascending order.
lastsize : the output size of the last axis to be transformed.
    If the corresponding input axis has size n, this can be 2*n-2 or 2*n-1.
forward : bool
    If `True`, a negative sign is used in the exponent, else a positive one.
inorm : int
    Normalization type
      0 : no normalization
      1 : divide by sqrt(N)
      2 : divide by N
    where N is the product of the lengths of the transformed output axes.
out : numpy.ndarray (real type with same accuracy as `a`)
    For the required shape, see the `Returns` section.
    Must not overlap with `a`.
    If None, a new array is allocated to store the output.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
numpy.ndarray (real type with same accuracy as `a`)
    The transformed data. The shape is identical to that of the input array,
    except for the axis that was transformed last, which has now `lastsize`
    entries.
)""";

const char *r2r_fftpack_DS = R"""(Performs a real-valued FFT using the FFTPACK storage scheme.

Parameters
----------
a : numpy.ndarray (any real type)
    The input data
axes : list of integers
    The axes along which the FFT is carried out.
    If not set, all axes will be transformed.
real2hermitian : bool
    if True, the input is purely real and the output will have Hermitian
    symmetry and be stored in FFTPACK's halfcomplex ordering, otherwise the
    opposite.
forward : bool
    If `True`, a negative sign is used in the exponent, else a positive one.
inorm : int
    Normalization type
      0 : no normalization
      1 : divide by sqrt(N)
      2 : divide by N
    where N is the length of `axis`.
out : numpy.ndarray (same shape and data type as `a`)
    May be identical to `a`, but if it isn't, it must not overlap with `a`.
    If None, a new array is allocated to store the output.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
numpy.ndarray (same shape and data type as `a`)
    The transformed data. The shape is identical to that of the input array.
)""";

const char *separable_hartley_DS = R"""(Performs a separable Hartley transform.
For every requested axis, a 1D forward Fourier transform is carried out, and
the real and imaginary parts of the result are added before the next axis is
processed.

Parameters
----------
a : numpy.ndarray (any real type)
    The input data
axes : list of integers
    The axes along which the transform is carried out.
    If not set, all axes will be transformed.
inorm : int
    Normalization type
      0 : no normalization
      1 : divide by sqrt(N)
      2 : divide by N
    where N is the product of the lengths of the transformed axes.
out : numpy.ndarray (same shape and data type as `a`)
    May be identical to `a`, but if it isn't, it must not overlap with `a`.
    If None, a new array is allocated to store the output.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
numpy.ndarray (same shape and data type as `a`)
    The transformed data
)""";

const char *genuine_hartley_DS = R"""(Performs a full Hartley transform.
A full Fourier transform is carried out over the requested axes, and the
sum of real and imaginary parts of the result is stored in the output
array. For a single transformed axis, this is identical to `separable_hartley`,
but when transforming multiple axes, the results are different.

Parameters
----------
a : numpy.ndarray (any real type)
    The input data
axes : list of integers
    The axes along which the transform is carried out.
    If not set, all axes will be transformed.
inorm : int
    Normalization type
      0 : no normalization
      1 : divide by sqrt(N)
      2 : divide by N
    where N is the product of the lengths of the transformed axes.
out : numpy.ndarray (same shape and data type as `a`)
    May be identical to `a`, but if it isn't, it must not overlap with `a`.
    If None, a new array is allocated to store the output.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
numpy.ndarray (same shape and data type as `a`)
    The transformed data
)""";

const char *dct_DS = R"""(Performs a discrete cosine transform.

Parameters
----------
a : numpy.ndarray (any real type)
    The input data
type : integer
    the type of DCT. Must be in [1; 4].
axes : list of integers
    The axes along which the transform is carried out.
    If not set, all axes will be transformed.
inorm : int
    Normalization type
      0 : no normalization
      1 : divide by sqrt(N)
      2 : divide by N
    where N is the product of n_i for every transformed axis i.
    n_i is 2*(<axis_length>-1 for type 1 and 2*<axis length>
    for types 2, 3, 4.
    Making the transform orthogonal involves the following additional steps
    for every 1D sub-transform:
      Type 1 : multiply first and last input value by sqrt(2)
               divide first and last output value by sqrt(2)
      Type 2 : divide first output value by sqrt(2)
      Type 3 : multiply first input value by sqrt(2)
      Type 4 : nothing
out : numpy.ndarray (same shape and data type as `a`)
    May be identical to `a`, but if it isn't, it must not overlap with `a`.
    If None, a new array is allocated to store the output.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).
ortho: bool
    Orthogonalize transform (defaults to ``inorm=1``)

Returns
-------
numpy.ndarray (same shape and data type as `a`)
    The transformed data
)""";

const char *dst_DS = R"""(Performs a discrete sine transform.

Parameters
----------
a : numpy.ndarray (any real type)
    The input data
type : integer
    the type of DST. Must be in [1; 4].
axes : list of integers
    The axes along which the transform is carried out.
    If not set, all axes will be transformed.
inorm : int
    Normalization type
      0 : no normalization
      1 : divide by sqrt(N)
      2 : divide by N
    where N is the product of n_i for every transformed axis i.
    n_i is 2*(<axis_length>+1 for type 1 and 2*<axis length>
    for types 2, 3, 4.
    Making the transform orthogonal involves the following additional steps
    for every 1D sub-transform:
      Type 1 : nothing
      Type 2 : divide first output value by sqrt(2)
      Type 3 : multiply first input value by sqrt(2)
      Type 4 : nothing
out : numpy.ndarray (same shape and data type as `a`)
    May be identical to `a`, but if it isn't, it must not overlap with `a`.
    If None, a new array is allocated to store the output.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).
ortho: bool
    Orthogonalize transform (defaults to ``inorm=1``)

Returns
-------
numpy.ndarray (same shape and data type as `a`)
    The transformed data
)""";

const char * good_size_DS = R"""(Returns a good length to pad an FFT to.

Parameters
----------
target : int
    Minimum transform length
real : bool, optional
    True if either input or output of FFT should be fully real.

Returns
-------
out : int
    The smallest fast size >= n

)""";


const char * prev_good_size_DS = R"""(Returns the largest FFT length less than target length.

Parameters
----------
target : int
    Maximum transform length
real : bool, optional
    True if either input or output of FFT should be fully real.

Returns
-------
out : int
    The largest fast length <= n

)""";

} // unnamed namespace

PYBIND11_MODULE(pypocketfft, m, py::mod_gil_not_used())
  {
  using namespace pybind11::literals;

  auto None = py::none();

  m.doc() = pypocketfft_DS;
  m.def("c2c", c2c, c2c_DS, "a"_a, "axes"_a=None, "forward"_a=true,
    "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
  m.def("r2c", r2c, r2c_DS, "a"_a, "axes"_a=None, "forward"_a=true,
    "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
  m.def("c2r", c2r, c2r_DS, "a"_a, "axes"_a=None, "lastsize"_a=0,
    "forward"_a=true, "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
  m.def("r2r_fftpack", r2r_fftpack, r2r_fftpack_DS, "a"_a, "axes"_a,
    "real2hermitian"_a, "forward"_a, "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
  m.def("separable_hartley", separable_hartley, separable_hartley_DS, "a"_a,
    "axes"_a=None, "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
  m.def("genuine_hartley", genuine_hartley, genuine_hartley_DS, "a"_a,
    "axes"_a=None, "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
  m.def("dct", dct, dct_DS, "a"_a, "type"_a, "axes"_a=None, "inorm"_a=0,
    "out"_a=None, "nthreads"_a=1, "ortho"_a=None);
  m.def("dst", dst, dst_DS, "a"_a, "type"_a, "axes"_a=None, "inorm"_a=0,
    "out"_a=None, "nthreads"_a=1, "ortho"_a=None);

  static PyMethodDef good_size_meth[] =
    {{"good_size", (PyCFunction)good_size,
      METH_VARARGS | METH_KEYWORDS, good_size_DS}, {0}};
  PyModule_AddFunctions(m.ptr(), good_size_meth);

  static PyMethodDef prev_good_size_meth[] =
    {{"prev_good_size", (PyCFunction)prev_good_size,
      METH_VARARGS | METH_KEYWORDS, prev_good_size_DS}, {0}};
  PyModule_AddFunctions(m.ptr(), prev_good_size_meth);
  }

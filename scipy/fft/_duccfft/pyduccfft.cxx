/*
 * This file is part of duccfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*
 *  Python interface.
 *
 *  Copyright (C) 2019-2025 Max-Planck-Society
 *  Copyright (C) 2019 Peter Bell
 *  \author Martin Reinecke
 *  \author Peter Bell
 */

#include <optional>
#include <complex>

#include "ducc0/fft/fft.h"
#include "ducc0/fft/fftnd_impl.h"
#include "ducc0/bindings/pybind_utils.h"

namespace {

using shape_t = ducc0::fmav_info::shape_t;
using stride_t = ducc0::fmav_info::stride_t;
using std::size_t;
using std::ptrdiff_t;

namespace py = ducc0::py;
using ducc0::CNpArr;
using ducc0::NpArr;
using ducc0::OptNpArr;

// Only instantiate long double transforms if they offer more precision
using ldbl_t = typename std::conditional<
  sizeof(long double)==sizeof(double), double, long double>::type;

using c64 = std::complex<float>;
using c128 = std::complex<double>;
using clong = std::complex<ldbl_t>;
using f32 = float;
using f64 = double;
using flong = ldbl_t;

using OptAxes = std::optional<std::vector<ptrdiff_t>>;

DUCC0_NOINLINE static size_t prev_good_size_cmplx(size_t n)
  {
    if (n<=12) return n;

    size_t bestfound = 1;
    for (size_t f11 = 1;f11 <= n; f11 *= 11)
      for (size_t f117 = f11; f117 <= n; f117 *= 7)
        for (size_t f1175 = f117; f1175 <= n; f1175 *= 5)
        {
          size_t x = f1175;
          while (x*2 <= n) x *= 2;
          if (x > bestfound) bestfound = x;
          while (true) 
          {
            if (x * 3 <= n) x *= 3;
            else if (x % 2 == 0) x /= 2;
            else break;
              
            if (x > bestfound) bestfound = x;
          }
        }
    return bestfound;
  }

/* returns the largest composite of 2, 3, 5 which is <= n */
DUCC0_NOINLINE static size_t prev_good_size_real(size_t n)
  {
    if (n<=6) return n;

    size_t bestfound = 1;
    for (size_t f5 = 1; f5 <= n; f5 *= 5)
    {
      size_t x = f5;
      while (x*2 <= n) x *= 2;
      if (x > bestfound) bestfound = x;
      while (true) 
      {
        if (x * 3 <= n) x *= 3;
        else if (x % 2 == 0) x /= 2;
        else break;
      
        if (x > bestfound) bestfound = x;
      }
    }
    return bestfound;
  }

shape_t makeaxes(const CNpArr &in, const OptAxes &axes)
  {
  if (!axes)
    {
    shape_t res(size_t(in.ndim()));
    for (size_t i=0; i<res.size(); ++i)
      res[i]=i;
    return res;
    }
  auto tmp=axes.value();
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
  if (ducc0::isPyarr<T1>(arr)) return func<double> args; \
  if (ducc0::isPyarr<T2>(arr)) return func<float> args;  \
  if (ducc0::isPyarr<T3>(arr)) return func<ldbl_t> args; \
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

template<typename T> NpArr c2c_internal(const CNpArr &in_,
  const OptAxes &axes_, bool forward, int inorm, const OptNpArr &out_,
  size_t nthreads)
  {
  auto axes = makeaxes(in_, axes_);
  auto in = ducc0::to_cfmav<std::complex<T>>(in_);
  auto [res_, res] = ducc0::get_OptNpArr_and_vfmav<std::complex<T>>(out_, in.shape(), "out");
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, in.shape(), axes);
  ducc0::c2c(in, res, axes, forward, fct, nthreads);
  }
  return std::move(res_);
  }

template<typename T> NpArr c2c_sym_internal(const CNpArr &in_,
  const OptAxes &axes_, bool forward, int inorm, const OptNpArr &out_,
  size_t nthreads)
  {
  auto axes = makeaxes(in_, axes_);
  auto in = ducc0::to_cfmav<T>(in_);
  auto [res_, res] = ducc0::get_OptNpArr_and_vfmav<std::complex<T>>(out_, in.shape(), "out");
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, in.shape(), axes);
  // select proper sub-array for FFT
  auto shp_half = res.shape();
  shp_half[axes.back()] = shp_half[axes.back()]/2+1;
  ducc0::vfmav<std::complex<T>> aout_half(res, shp_half, res.stride());
  ducc0::r2c(in, aout_half, axes, forward, fct, nthreads);
  // now fill in second half
  using namespace ducc0::detail_fft;
  hermiteHelper(0, 0, 0, 0, res, res, axes, [](const complex<T> &c, complex<T> &, complex<T> &c1)
    {
    c1 = conj(c);
    }, nthreads);
  }
  return std::move(res_);
  }

NpArr c2c(const CNpArr &a, const OptAxes &axes_, bool forward,
  int inorm, const OptNpArr &out_, size_t nthreads)
  {
  if (a.dtype().kind() == 'c')
    DISPATCH(a, c128, c64, clong, c2c_internal, (a, axes_, forward,
             inorm, out_, nthreads))

  DISPATCH(a, f64, f32, flong, c2c_sym_internal, (a, axes_, forward,
           inorm, out_, nthreads))
  }

template<typename T> NpArr r2c_internal(const CNpArr &in_,
  const OptAxes &axes_, bool forward, int inorm, const OptNpArr &out_,
  size_t nthreads)
  {
  auto axes = makeaxes(in_, axes_);
  auto in = ducc0::to_cfmav<T>(in_);
  auto dims_out(in.shape());
  dims_out[axes.back()] = (dims_out[axes.back()]>>1)+1;
  auto [res_, res] = ducc0::get_OptNpArr_and_vfmav<std::complex<T>>(out_, dims_out, "out");
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, in.shape(), axes);
  ducc0::r2c(in, res, axes, forward, fct, nthreads);
  }
  return res_;
  }

NpArr r2c(const CNpArr &in, const OptAxes &axes_, bool forward,
  int inorm, const OptNpArr &out_, size_t nthreads)
  {
  DISPATCH(in, f64, f32, flong, r2c_internal, (in, axes_, forward, inorm, out_,
    nthreads))
  }

template<typename T> NpArr r2r_fftpack_internal(const CNpArr &in_,
  const OptAxes &axes_, bool real2hermitian, bool forward, int inorm,
  const OptNpArr &out_, size_t nthreads)
  {
  auto axes = makeaxes(in_, axes_);
  auto in = ducc0::to_cfmav<T>(in_);
  auto [res_, res] = ducc0::get_OptNpArr_and_vfmav<T>(out_, in.shape(), "out");
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, in.shape(), axes);
  ducc0::r2r_fftpack(in, res, axes, real2hermitian, forward, fct, nthreads);
  }
  return res_;
  }

NpArr r2r_fftpack(const CNpArr &in, const OptAxes &axes_,
  bool real2hermitian, bool forward, int inorm, const OptNpArr &out_,
  size_t nthreads)
  {
  DISPATCH(in, f64, f32, flong, r2r_fftpack_internal, (in, axes_,
    real2hermitian, forward, inorm, out_, nthreads))
  }

template<typename T> NpArr dct_internal(const CNpArr &in_,
  const OptAxes &axes_, int type, int inorm, const OptNpArr &out_,
  size_t nthreads, bool ortho)
  {
  auto axes = makeaxes(in_, axes_);
  auto in = ducc0::to_cfmav<T>(in_);
  auto [res_, res] = ducc0::get_OptNpArr_and_vfmav<T>(out_, in.shape(), "out");
  {
  py::gil_scoped_release release;
  T fct = (type==1) ? norm_fct<T>(inorm, in.shape(), axes, 2, -1)
                    : norm_fct<T>(inorm, in.shape(), axes, 2);
  ducc0::dct(in, res, axes, type, fct, ortho, nthreads);
  }
  return res_;
  }

NpArr dct(const CNpArr &in, int type, const OptAxes &axes_,
  int inorm, const OptNpArr &out_, size_t nthreads, const py::object &ortho_obj)
  {
  bool ortho=inorm==1;
  if (!ortho_obj.is_none())
    ortho=ortho_obj.cast<bool>();

  if ((type<1) || (type>4)) throw std::invalid_argument("invalid DCT type");
  DISPATCH(in, f64, f32, flong, dct_internal, (in, axes_, type, inorm, out_,
    nthreads, ortho))
  }

template<typename T> NpArr dst_internal(const CNpArr &in_,
  const OptAxes &axes_, int type, int inorm, const OptNpArr &out_,
  size_t nthreads, bool ortho)
  {
  auto axes = makeaxes(in_, axes_);
  auto in = ducc0::to_cfmav<T>(in_);
  auto [res_, res] = ducc0::get_OptNpArr_and_vfmav<T>(out_, in.shape(), "out");
  {
  py::gil_scoped_release release;
  T fct = (type==1) ? norm_fct<T>(inorm, in.shape(), axes, 2, 1)
                    : norm_fct<T>(inorm, in.shape(), axes, 2);
  ducc0::dst(in, res, axes, type, fct, ortho, nthreads);
  }
  return res_;
  }

NpArr dst(const CNpArr &in, int type, const OptAxes &axes_,
  int inorm, const OptNpArr &out_, size_t nthreads, const py::object &ortho_obj)
  {
  bool ortho=inorm==1;
  if (!ortho_obj.is_none())
    ortho=ortho_obj.cast<bool>();

  if ((type<1) || (type>4)) throw std::invalid_argument("invalid DST type");
  DISPATCH(in, f64, f32, flong, dst_internal, (in, axes_, type, inorm,
    out_, nthreads, ortho))
  }

template<typename T> NpArr c2r_internal(const CNpArr &in_,
  const OptAxes &axes_, size_t lastsize, bool forward, int inorm,
  const OptNpArr &out_, size_t nthreads)
  {
  auto axes = makeaxes(in_, axes_);
  size_t axis = axes.back();
  auto in = ducc0::to_cfmav<std::complex<T>>(in_);
  shape_t dims_out=in.shape();
  if (lastsize==0) lastsize=2*in.shape(axis)-1;
  if ((lastsize/2) + 1 != in.shape(axis))
    throw std::invalid_argument("bad lastsize");
  dims_out[axis] = lastsize;
  auto [res_, res] = ducc0::get_OptNpArr_and_vfmav<T>(out_, dims_out, "out");
  {
  py::gil_scoped_release release;
  T fct = norm_fct<T>(inorm, dims_out, axes);
  ducc0::c2r(in, res, axes, forward, fct, nthreads);
  }
  return res_;
  }

NpArr c2r(const CNpArr &in, const OptAxes &axes_, size_t lastsize,
  bool forward, int inorm, const OptNpArr &out_, size_t nthreads)
  {
  DISPATCH(in, c128, c64, clong, c2r_internal, (in, axes_, lastsize, forward,
    inorm, out_, nthreads))
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
  return PyLong_FromSize_t(
    real ? ducc0::good_size_real(n) : ducc0::good_size_complex(n));
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
  return PyLong_FromSize_t(
    real ? prev_good_size_real(n) : prev_good_size_cmplx(n));
  }

const char *pyduccfft_DS = R"""(Fast Fourier and trogonomeric transforms.

This module supports
- single, double, and long double precision
- complex and real-valued transforms
- multi-dimensional transforms

The code will use SSE2 and AVX vector instructions for faster execution if
these are supported by the CPU and were enabled during compilation.
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

PYBIND11_MODULE(pyduccfft, m, py::mod_gil_not_used())
  {
  using namespace pybind11::literals;

  auto None = py::none();

  m.doc() = pyduccfft_DS;
  m.def("c2c", c2c, c2c_DS, "a"_a, "axes"_a=None, "forward"_a=true,
    "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
  m.def("r2c", r2c, r2c_DS, "a"_a, "axes"_a=None, "forward"_a=true,
    "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
  m.def("c2r", c2r, c2r_DS, "a"_a, "axes"_a=None, "lastsize"_a=0,
    "forward"_a=true, "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
  m.def("r2r_fftpack", r2r_fftpack, r2r_fftpack_DS, "a"_a, "axes"_a,
    "real2hermitian"_a, "forward"_a, "inorm"_a=0, "out"_a=None, "nthreads"_a=1);
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

/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*
 *  Python interface.
 *
 *  Copyright (C) 2019 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "pocketfft.h"
#include "scipy_features.h"

//
// Python interface
//

namespace {

using namespace std;
using namespace pocketfft;

namespace py = pybind11;

auto c64 = py::dtype("complex64");
auto c128 = py::dtype("complex128");
auto c256 = py::dtype("longcomplex");
auto f32 = py::dtype("float32");
auto f64 = py::dtype("float64");
auto f128 = py::dtype("longfloat");

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

shape_t makeaxes(const py::array &in, py::object axes)
  {
  if (axes.is(py::none()))
    {
    shape_t res(size_t(in.ndim()));
    for (size_t i=0; i<res.size(); ++i)
      res[i]=i;
    return res;
    }
  auto tmp=axes.cast<std::vector<ptrdiff_t>>();
  auto ndim = in.ndim();
  if ((tmp.size()>size_t(ndim)) || (tmp.size()==0))
    throw runtime_error("bad axes argument");
  for (auto& sz: tmp)
    {
    if (sz<0)
      sz += ndim;
    if ((sz>=ndim) || (sz<0))
      throw invalid_argument("axes exceeds dimensionality of output");
    }
  return shape_t(tmp.begin(), tmp.end());
  }

#define DISPATCH(arr, T1, T2, T3, func, args) \
  auto dtype = arr.dtype(); \
  if (dtype.is(T1)) return func<double> args; \
  if (dtype.is(T2)) return func<float> args; \
  if (dtype.is(T3)) return func<longfloat_t> args; \
  throw runtime_error("unsupported data type");

template<typename T> py::array xfftn_internal(const py::array &in,
  const shape_t &axes, norm_t norm, bool inplace, bool fwd, size_t nthreads)
  {
  auto dims(copy_shape(in));
  py::array res = inplace ? in : py::array_t<complex<T>>(dims);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const complex<T> *>(in.data());
  auto d_out=reinterpret_cast<complex<T> *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  auto fct = norm_fct<T>(norm, dims, axes);
  c2c(dims, s_in, s_out, axes, fwd, d_in, d_out, fct, nthreads);
  }
  return res;
  }

py::array xfftn(const py::array &in, py::object axes, norm_t norm, bool inplace,
                bool fwd, size_t nthreads)
  {
    py::array a = asfarray(in, inplace);
    auto dtype = a.dtype();
#define X_(NP_TYPE, FUNC)                                               \
    if (dtype.is(NP_TYPE)) \
      return FUNC(a, makeaxes(a, axes), norm, inplace, fwd, nthreads)
    X_(c64,  xfftn_internal<float>);
    X_(c128, xfftn_internal<double>);
    X_(c256, xfftn_internal<longfloat_t>);
    X_(f32,  xfftn_real<float>);
    X_(f64,  xfftn_real<double>);
    X_(f128, xfftn_real<longfloat_t>);
#undef X_
    throw runtime_error("unsupported data type");
  }

py::array fftn(const py::array &a, py::object axes, norm_t norm, bool inplace,
  size_t nthreads)
  { return xfftn(a, axes, norm, inplace, true, nthreads); }

py::array ifftn(const py::array &a, py::object axes, norm_t norm, bool inplace,
  size_t nthreads)
  { return xfftn(a, axes, norm, inplace, false, nthreads); }

template<typename T> py::array rfftn_internal(const py::array &in,
  py::object axes_, norm_t norm, size_t nthreads)
  {
  auto axes = makeaxes(in, axes_);
  auto dims_in(copy_shape(in)), dims_out(dims_in);
  dims_out[axes.back()] = (dims_out[axes.back()]>>1)+1;
  py::array res = py::array_t<complex<T>>(dims_out);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<complex<T> *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  auto fct = norm_fct<T>(norm, dims_in, axes);
  r2c(dims_in, s_in, s_out, axes, d_in, d_out, fct, nthreads);
  }
  return res;
  }

py::array rfftn(const py::array &in, py::object axes_, norm_t norm,
  size_t nthreads)
  {
  py::array a = asfarray(in);
  DISPATCH(a, f64, f32, f128, rfftn_internal, (a, axes_, norm, nthreads))
  }

template<typename T> py::array xrfft_scipy(const py::array &in,
  size_t axis, norm_t norm, bool inplace, bool fwd, size_t nthreads)
  {
  auto dims(copy_shape(in));
  py::array res = inplace ? in : py::array_t<T>(dims);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<T *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  auto fct = norm_fct<T>(norm, dims, {axis});
  r2r_fftpack(dims, s_in, s_out, axis, fwd, d_in, d_out, fct, nthreads);
  }
  return res;
  }

py::array rfft_scipy(const py::array &in, size_t axis, norm_t norm, bool inplace,
  size_t nthreads)
  {
  py::array a = asfarray(in, inplace);
  DISPATCH(a, f64, f32, f128, xrfft_scipy, (a, axis, norm, inplace, true,
    nthreads))
  }

py::array irfft_scipy(const py::array &in, size_t axis, norm_t norm,
  bool inplace, size_t nthreads)
  {
  py::array a = asfarray(in, inplace);
  DISPATCH(a, f64, f32, f128, xrfft_scipy, (a, axis, norm, inplace, false,
    nthreads))
  }
template<typename T> py::array irfftn_internal(const py::array &in,
  py::object axes_, size_t lastsize, norm_t norm, size_t nthreads)
  {
  auto axes = makeaxes(in, axes_);
  size_t axis = axes.back();
  shape_t dims_in(copy_shape(in)), dims_out=dims_in;
  if (lastsize==0) lastsize=2*dims_in[axis]-1;
  if ((lastsize/2) + 1 != dims_in[axis])
    throw runtime_error("bad lastsize");
  dims_out[axis] = lastsize;
  py::array res = py::array_t<T>(dims_out);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const complex<T> *>(in.data());
  auto d_out=reinterpret_cast<T *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  auto fct = norm_fct<T>(norm, dims_out, axes);
  c2r(dims_out, s_in, s_out, axes, d_in, d_out, fct, nthreads);
  }
  return res;
  }

py::array irfftn(const py::array &in, py::object axes_, size_t lastsize,
  norm_t norm, size_t nthreads)
  {
  py::array a = asfarray(in);
  DISPATCH(a, c128, c64, c256, irfftn_internal, (a, axes_, lastsize, norm,
    nthreads))
  }

template<typename T> py::array hartley_internal(const py::array &in,
  py::object axes_, norm_t norm, bool inplace, size_t nthreads)
  {
  auto dims(copy_shape(in));
  py::array res = inplace ? in : py::array_t<T>(dims);
  auto axes = makeaxes(in, axes_);
  auto s_in=copy_strides(in);
  auto s_out=copy_strides(res);
  auto d_in=reinterpret_cast<const T *>(in.data());
  auto d_out=reinterpret_cast<T *>(res.mutable_data());
  {
  py::gil_scoped_release release;
  auto fct = norm_fct<T>(norm, dims, axes);
  r2r_hartley(dims, s_in, s_out, axes, d_in, d_out, fct, nthreads);
  }
  return res;
  }

py::array hartley(const py::array &in, py::object axes_, norm_t norm,
  bool inplace, size_t nthreads)
  {
  DISPATCH(in, f64, f32, f128, hartley_internal, (in, axes_, norm, inplace,
    nthreads))
  }

template<typename T>py::array complex2hartley(const py::array &in,
  const py::array &tmp, py::object axes_, bool inplace)
  {
  using namespace pocketfft::detail;
  size_t ndim = size_t(in.ndim());
  auto dims_out(copy_shape(in));
  py::array out = inplace ? in : py::array_t<T>(dims_out);
  ndarr<cmplx<T>> atmp(tmp.data(), copy_shape(tmp), copy_strides(tmp));
  ndarr<T> aout(out.mutable_data(), copy_shape(out), copy_strides(out));
  auto axes = makeaxes(in, axes_);
  {
  py::gil_scoped_release release;
  size_t axis = axes.back();
  multi_iter<1,cmplx<T>,T> it(atmp, aout, axis);
  vector<bool> swp(ndim,false);
  for (auto i: axes)
    if (i!=axis)
      swp[i] = true;
  while(it.remaining()>0)
    {
    ptrdiff_t rofs = 0;
    for (size_t i=0; i<it.pos.size(); ++i)
      {
      if (i==axis) continue;
      if (!swp[i])
        rofs += ptrdiff_t(it.pos[i])*it.oarr.stride(i);
      else
        {
        auto x = ptrdiff_t((it.pos[i]==0) ? 0 : it.iarr.shape(i)-it.pos[i]);
        rofs += x*it.oarr.stride(i);
        }
      }
    it.advance(1);
    for (size_t i=0; i<it.length_in(); ++i)
      {
      auto re = it.in(i).r;
      auto im = it.in(i).i;
      auto rev_i = ptrdiff_t((i==0) ? 0 : it.length_out()-i);
      it.out(i) = re+im;
      aout[rofs + rev_i*it.stride_out()] = re-im;
      }
    }
  }
  return out;
  }

py::array mycomplex2hartley(const py::array &in,
  const py::array &tmp, py::object axes_, bool inplace)
  {
  DISPATCH(in, f64, f32, f128, complex2hartley, (in, tmp, axes_, inplace))
  }

py::array hartley2(const py::array &in, py::object axes_, norm_t norm,
  bool inplace, size_t nthreads)
  {
  return mycomplex2hartley(in, rfftn(in, axes_, norm, nthreads), axes_,
    inplace);
  }

const char *pypocketfft_DS = R"""(Fast Fourier and Hartley transforms.

This module supports
- single and double precision
- complex and real-valued transforms
- multi-dimensional transforms

For two- and higher-dimensional transforms the code will use SSE2 and AVX
vector instructions for faster execution if these are supported by the CPU and
were enabled during compilation.
)""";

const char *fftn_DS = R"""(
Performs a forward complex FFT.

Parameters
----------
a : numpy.ndarray (np.complex64 or np.complex128)
    The input data
axes : list of integers
    The axes along which the FFT is carried out.
    If not set, all axes will be transformed.
fct : float
    Normalization factor
inplace : bool
    if False, returns the result in a new array and leaves the input unchanged.
    if True, stores the result in the input array and returns a handle to it.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
np.ndarray (same shape and data type as a)
    The transformed data.
)""";

const char *ifftn_DS = R"""(Performs a backward complex FFT.

Parameters
----------
a : numpy.ndarray (np.complex64 or np.complex128)
    The input data
axes : list of integers
    The axes along which the FFT is carried out.
    If not set, all axes will be transformed.
fct : float
    Normalization factor
inplace : bool
    if False, returns the result in a new array and leaves the input unchanged.
    if True, stores the result in the input array and returns a handle to it.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
np.ndarray (same shape and data type as a)
    The transformed data
)""";

const char *rfftn_DS = R"""(Performs a forward real-valued FFT.

Parameters
----------
a : numpy.ndarray (np.float32 or np.float64)
    The input data
axes : list of integers
    The axes along which the FFT is carried out.
    If not set, all axes will be transformed in ascending order.
fct : float
    Normalization factor
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
np.ndarray (np.complex64 or np.complex128)
    The transformed data. The shape is identical to that of the input array,
    except for the axis that was transformed last. If the length of that axis
    was n on input, it is n//2+1 on output.
)""";

const char *rfft_scipy_DS = R"""(Performs a forward real-valued FFT.

Parameters
----------
a : numpy.ndarray (np.float32 or np.float64)
    The input data
axis : int
    The axis along which the FFT is carried out.
fct : float
    Normalization factor
inplace : bool
    if False, returns the result in a new array and leaves the input unchanged.
    if True, stores the result in the input array and returns a handle to it.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
np.ndarray (np.float32 or np.float64)
    The transformed data. The shape is identical to that of the input array.
    Along the transformed axis, values are arranged in
    FFTPACK half-complex order, i.e. `a[0].re, a[1].re, a[1].im, a[2].re ...`.
)""";

const char *irfftn_DS = R"""(Performs a backward real-valued FFT.

Parameters
----------
a : numpy.ndarray (np.complex64 or np.complex128)
    The input data
axes : list of integers
    The axes along which the FFT is carried out.
    If not set, all axes will be transformed in ascending order.
lastsize : the output size of the last axis to be transformed.
    If the corresponding input axis has size n, this can be 2*n-2 or 2*n-1.
fct : float
    Normalization factor
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
np.ndarray (np.float32 or np.float64)
    The transformed data. The shape is identical to that of the input array,
    except for the axis that was transformed last, which has now `lastsize`
    entries.
)""";

const char *irfft_scipy_DS = R"""(Performs a backward real-valued FFT.

Parameters
----------
a : numpy.ndarray (np.float32 or np.float64)
    The input data. Along the transformed axis, values are expected in
    FFTPACK half-complex order, i.e. `a[0].re, a[1].re, a[1].im, a[2].re ...`.
axis : int
    The axis along which the FFT is carried out.
fct : float
    Normalization factor
inplace : bool
    if False, returns the result in a new array and leaves the input unchanged.
    if True, stores the result in the input array and returns a handle to it.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
np.ndarray (np.float32 or np.float64)
    The transformed data. The shape is identical to that of the input array.
)""";

const char *hartley_DS = R"""(Performs a Hartley transform.
For every requested axis, a 1D forward Fourier transform is carried out,
and the sum of real and imaginary parts of the result is stored in the output
array.

Parameters
----------
a : numpy.ndarray (np.float32 or np.float64)
    The input data
axes : list of integers
    The axes along which the transform is carried out.
    If not set, all axes will be transformed.
fct : float
    Normalization factor
inplace : bool
    if False, returns the result in a new array and leaves the input unchanged.
    if True, stores the result in the input array and returns a handle to it.
nthreads : int
    Number of threads to use. If 0, use the system default (typically governed
    by the `OMP_NUM_THREADS` environment variable).

Returns
-------
np.ndarray (same shape and data type as a)
    The transformed data
)""";

} // unnamed namespace

PYBIND11_MODULE(pypocketfft, m)
  {
  using namespace pybind11::literals;

  m.doc() = pypocketfft_DS;
  py::enum_<norm_t>(m, "norm_t")
    .value("none",  norm_t::none)
    .value("ortho", norm_t::ortho)
    .value("size",  norm_t::size);
  m.def("fftn",&fftn, fftn_DS, "a"_a, "axes"_a=py::none(), "norm"_a=norm_t::none,
        "inplace"_a=false, "nthreads"_a=1);
  m.def("ifftn",&ifftn, ifftn_DS, "a"_a, "axes"_a=py::none(),
        "norm"_a=norm_t::none, "inplace"_a=false, "nthreads"_a=1);
  m.def("rfftn",&rfftn, rfftn_DS, "a"_a, "axes"_a=py::none(),
        "norm"_a=norm_t::none, "nthreads"_a=1);
  m.def("rfft_scipy",&rfft_scipy, rfft_scipy_DS, "a"_a, "axis"_a,
        "norm"_a=norm_t::none, "inplace"_a=false, "nthreads"_a=1);
  m.def("irfftn",&irfftn, irfftn_DS, "a"_a, "axes"_a=py::none(), "lastsize"_a=0,
    "norm"_a=norm_t::none, "nthreads"_a=1);
  m.def("irfft_scipy",&irfft_scipy, irfft_scipy_DS, "a"_a, "axis"_a,
        "norm"_a=norm_t::none, "inplace"_a=false, "nthreads"_a=1);
  m.def("hartley",&hartley, hartley_DS, "a"_a, "axes"_a=py::none(),
        "norm"_a=norm_t::none, "inplace"_a=false, "nthreads"_a=1);
  m.def("hartley2",&hartley2, "a"_a, "axes"_a=py::none(), "norm"_a=norm_t::none,
    "inplace"_a=false, "nthreads"_a=1);
  m.def("complex2hartley",&mycomplex2hartley, "in"_a, "tmp"_a, "axes"_a,
    "inplace"_a=false);
  }

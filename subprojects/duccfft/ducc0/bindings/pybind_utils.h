/* SPDX-License-Identifier: BSD-3-Clause OR GPL-2.0-or-later */

/*
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.
* Neither the name of the copyright holder nor the names of its contributors may
  be used to endorse or promote products derived from this software without
  specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* Copyright (C) 2020-2025 Max-Planck-Society
   Author: Martin Reinecke */

#ifndef DUCC0_PYBIND_UTILS_H
#define DUCC0_PYBIND_UTILS_H

#include <cstddef>
#include <string>
#include <array>
#include <vector>
#include <optional>
#include <variant>
#include <tuple>
#ifdef DUCC0_USE_NANOBIND
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/function.h>
#else
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#endif

#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/mav.h"
#include "ducc0/infra/misc_utils.h"

namespace ducc0 {

#ifdef DUCC0_USE_NANOBIND
namespace py = nanobind;
#else
namespace py = pybind11;
#endif

namespace detail_pybind {

using namespace std;

using shape_t=fmav_info::shape_t;
using stride_t=fmav_info::stride_t;

static const auto None = py::none();

#ifdef DUCC0_USE_NANOBIND
using NpArr = py::ndarray<py::numpy, py::device::cpu>;
using CNpArr = py::ndarray<py::numpy, py::ro, py::device::cpu>;
template<typename T> using NpArrT = py::ndarray<py::numpy, py::device::cpu, T>;
template<typename T> using CNpArrT = py::ndarray<py::numpy, py::ro, py::device::cpu, T>;
#else
using NpArr = py::array;
using CNpArr = py::array;
template<typename T> using NpArrT = py::array_t<T>;
template<typename T> using CNpArrT = py::array_t<T>;
#endif

using OptNpArr = optional<NpArr>;
using OptCNpArr = optional<CNpArr>;

static inline string makeSpec(const string &name)
  { return (name=="") ? "" : name+": "; }

template<typename T> bool isPyarr(const CNpArr &obj)
#ifdef DUCC0_USE_NANOBIND
  { return obj.dtype()==py::dtype<T>(); }
#else
  { return py::isinstance<py::array_t<T>>(obj); }
#endif

static inline shape_t copy_shape(const CNpArr &arr, const string &/*spec*/="")
  {
  shape_t res(size_t(arr.ndim()));
  for (size_t i=0; i<res.size(); ++i)
    res[i] = size_t(arr.shape(int(i)));
  return res;
  }

template<typename T, bool rw> stride_t copy_strides(const CNpArr &arr,
  const string &spec="")
  {
  stride_t res(size_t(arr.ndim()));
  bool zerosized = false;
  if constexpr(rw)
    for (size_t i=0; i<res.size(); ++i)
      if (arr.shape(i)==0)
        zerosized=true;
  for (size_t i=0; i<res.size(); ++i)
    {
#ifdef DUCC0_USE_NANOBIND
    auto tmp = arr.stride(int(i));
    res[i] = tmp;
#else
    auto tmp = arr.strides(int(i));
    constexpr auto st = ptrdiff_t(sizeof(T));
    MR_assert((tmp/st)*st==tmp, spec, "bad stride");
    res[i] = tmp/st;
#endif
    if constexpr(rw)
      if (!zerosized)  // if the array has no elements, we needn't worry
        MR_assert((arr.shape(int(i))==1) || (tmp!=0),
          spec, "detected zero stride in writable array");
    }
  return res;
  }

template<typename T> cfmav<T> to_cfmav(const CNpArr &obj, const string &name="")
  {
  const auto spec = makeSpec(name);
  MR_assert(isPyarr<const T>(obj), "data type mismatch");
  return cfmav<T>(reinterpret_cast<const T *>(obj.data()),
    copy_shape(obj, spec), copy_strides<T,false>(obj, spec));
  }
template<typename T> cfmav<T> to_cfmav(const CNpArrT<T> &obj, const string &name="")
  { return to_cfmav<T>(CNpArr(obj), name); }

template<typename T, size_t ndim> cmav<T,ndim> to_cmav(const CNpArr &obj,
  const string &name="")
  { return cmav<T,ndim>(to_cfmav<T>(obj, name)); }

template<typename T, size_t ndim> cmav<T,ndim> to_cmav(const CNpArrT<T> &obj,
  const string &name="")
  { return cmav<T,ndim>(to_cfmav<T>(obj, name)); }

static inline auto extend_axes(fmav_info &info, size_t ndim, const string &name="")
  {
  const auto spec = makeSpec(name);
  MR_assert(info.ndim()<=ndim, spec, "array has too many dimensions");
  shape_t newshape(ndim, 1);
  stride_t newstride(ndim, 0);
  size_t add=ndim-info.ndim();
  for (size_t i=0; i<info.ndim(); ++i)
    { newshape[i+add]=info.shape(i); newstride[i+add]=info.stride(i); }
  return make_tuple(newshape, newstride);
  }

template<typename T> cfmav<T> to_cfmav_with_optional_leading_dimensions(const CNpArr &obj, size_t ndim,
  const string &name="")
  {
  auto tmp = to_cfmav<T>(obj, name);
  auto [newshape, newstride] = extend_axes(tmp, ndim, name);
  return cfmav<T>(tmp.data(), newshape, newstride);
  }
template<typename T, size_t ndim> cmav<T,ndim> to_cmav_with_optional_leading_dimensions(const CNpArr &obj,
  const string &name="")
  { return cmav<T,ndim>(to_cfmav_with_optional_leading_dimensions<T>(obj, ndim, name)); }

#ifdef DUCC0_USE_NANOBIND  // with nanobind, we need extra functions working on nonconst Python arrays
template<typename T> bool isPyarr(const NpArr &obj)
  { return isPyarr<T>(CNpArr(obj)); }
template<typename T> cfmav<T> to_cfmav(const NpArr &obj, const string &name="")
  { return to_cfmav<T>(CNpArr(obj), name); }
template<typename T, size_t ndim> cmav<T,ndim> to_cmav(const NpArr &obj, const string &name="")
  { return to_cmav<T,ndim>(CNpArr(obj), name); }
template<typename T> cfmav<T> to_cfmav_with_optional_leading_dimensions(const NpArr &obj, size_t ndim,
  const string &name="")
  { return to_cfmav_with_optional_leading_dimensions<T>(CNpArr(obj), ndim, name); }
template<typename T, size_t ndim> cmav<T,ndim> to_cmav_with_optional_leading_dimensions(const NpArr &obj,
  const string &name="")
  { return to_cmav_with_optional_leading_dimensions<T, ndim>(CNpArr(obj), name); }
#endif

template<typename T> vfmav<T> to_vfmav(const NpArr &obj, const string &name="")
  {
  const auto spec = makeSpec(name);
  MR_assert(isPyarr<T>(obj), "data type mismatch");
#ifdef DUCC0_USE_NANOBIND
  return vfmav<T>(reinterpret_cast<T *>(obj.data()),
    copy_shape(CNpArr(obj), spec), copy_strides<T,true>(CNpArr(obj), spec));
#else
  auto arr = NpArrT<T>(obj);
  return vfmav<T>(reinterpret_cast<T *>(arr.mutable_data()),
    copy_shape(CNpArr(obj), spec), copy_strides<T,true>(CNpArr(obj), spec));
#endif
  }
  template<typename T> vfmav<T> to_vfmav(const NpArrT<T> &obj, const string &name="")
  { return to_vfmav<T>(NpArr(obj), name); }


template<typename T, size_t ndim> vmav<T,ndim> to_vmav(const NpArr &obj,
  const string &name="")
  { return vmav<T,ndim>(to_vfmav<T>(obj, name)); }
template<typename T, size_t ndim> vmav<T,ndim> to_vmav(const NpArrT<T> &obj,
  const string &name="")
  { return to_vmav<T,ndim>(NpArr(obj), name); }

template<typename T> vfmav<T> to_vfmav_with_optional_leading_dimensions(const NpArr &obj, size_t ndim,
  const string &name="")
  {
  auto tmp = to_vfmav<T>(obj, name);
  auto [newshape, newstride] = extend_axes(tmp, ndim, name);
  return vfmav<T>(tmp.data(), newshape, newstride);
  }
template<typename T, size_t ndim> vmav<T,ndim> to_vmav_with_optional_leading_dimensions(const NpArr &obj,
  const string &name="")
  { return vmav<T,ndim>(to_vfmav_with_optional_leading_dimensions<T>(obj, ndim, name)); }

template<typename T> void zero_Pyarr(const NpArr &arr, size_t nthreads=1)
  { mav_apply([](T &v){ v=T(0); }, nthreads, to_vfmav<T>(arr)); }

template<typename T> NpArr make_Pyarr(const shape_t &dims, bool zero=false, size_t nthreads=1)
  {
#ifdef DUCC0_USE_NANOBIND
  auto *res = new vfmav<T>(dims, PAGE_IN(nthreads));
  py::capsule owner(res, [](void *p) noexcept {
      delete reinterpret_cast<vfmav<T> *>(p);
    });
  NpArr res_(NpArrT<T>(res->data(), dims.size(), dims.data(), owner));
#else
  auto res_=NpArr(NpArrT<T>(dims));
  page_in_memory(reinterpret_cast<T *>(res_.mutable_data()), res_.size(), nthreads);
#endif
  if (zero) zero_Pyarr<T>(res_, nthreads);
  return res_;
  }
template<typename T, size_t ndim> NpArr make_Pyarr
  (const array<size_t,ndim> &dims, bool zero=false)
  { return make_Pyarr<T>(shape_t(dims.begin(), dims.end()), zero); }
template<typename T, size_t ndim> auto make_Pyarr_and_vmav
  (const shape_t &dims, bool zero=false, size_t nthreads=1)
  {
  auto res_py = make_Pyarr<T>(dims, zero, nthreads);
  auto res_mav = to_vmav<T,ndim>(res_py);
  return std::make_tuple(res_py, res_mav);
  }
template<typename T> auto make_Pyarr_and_vfmav
  (const shape_t &dims, bool zero=false, size_t nthreads=1)
  {
  auto res_py = make_Pyarr<T>(dims, zero, nthreads);
  auto res_vfmav = to_vfmav<T>(res_py);
  return std::make_tuple(res_py, res_vfmav);
  }

template<typename T> NpArr make_noncritical_Pyarr(const shape_t &shape, bool zero=false, size_t nthreads=1)
  {
  auto ndim = shape.size();
  if (ndim==1) return make_Pyarr<T>(shape);
  auto shape2 = noncritical_shape(shape, sizeof(T));
  NpArr res;
#ifdef DUCC0_USE_NANOBIND
  auto *tmp = new vfmav<T>(shape2, PAGE_IN(nthreads));
  py::capsule owner(tmp, [](void *p) noexcept {
      delete reinterpret_cast<vfmav<T> *>(p);
    });
  // nanobind strides are always int64_t, but on some platforms ptrdiff_t is a
  // different type, so we must be careful
  if constexpr(is_same<ptrdiff_t, int64_t>::value)
    {
    NpArrT<T> res_(tmp->data(), shape.size(), shape.data(), owner, tmp->stride().data());
    res = NpArr(res_);
    }
  else
    {
    std::vector<int64_t> stmp;
    for (auto x: tmp->stride()) stmp.push_back(int64_t(x));
    NpArrT<T> res_(tmp->data(), shape.size(), shape.data(), owner, stmp.data());
    res = NpArr(res_);
    }
#else
  NpArrT<T> tmp(shape2);
  page_in_memory(reinterpret_cast<T *>(tmp.mutable_data()), tmp.size(), nthreads);
  py::list slices;
  for (size_t i=0; i<ndim; ++i)
    slices.append(py::slice(0, shape[i], 1));
  NpArrT<T> res_(tmp[py::tuple(slices)]);
  res = NpArr(res_);
#endif
  if (zero) zero_Pyarr<T>(res, nthreads);
  return res;
  }

template<typename T> NpArr get_OptNpArr(const OptNpArr &arr_,
  const shape_t &dims, const string &name="", size_t nthreads=1)
  {
  if (!arr_) return make_Pyarr<T>(dims, false, nthreads);
  const auto spec = makeSpec(name);
  auto val = arr_.value();
  MR_assert(isPyarr<T>(val), spec, "incorrect data type");
  MR_assert(dims.size()==size_t(val.ndim()), spec, "dimension mismatch");
  for (size_t i=0; i<dims.size(); ++i)
    MR_assert(dims[i]==size_t(val.shape(int(i))), spec, "dimension mismatch");
  return val;
  }
template<typename T> auto get_OptNpArr_and_vfmav(const OptNpArr &arr_,
  const shape_t &dims, const string &name="", size_t nthreads=1)
  {
  if (!arr_) return make_Pyarr_and_vfmav<T>(dims, false, nthreads);
  const auto spec = makeSpec(name);
  auto val = arr_.value();
  MR_assert(isPyarr<T>(val), spec, "incorrect data type");
  MR_assert(dims.size()==size_t(val.ndim()), spec, "dimension mismatch");
  for (size_t i=0; i<dims.size(); ++i)
    MR_assert(dims[i]==size_t(val.shape(int(i))), spec, "dimension mismatch");
  auto res_vfmav = to_vfmav<T>(val);
  return std::make_tuple(val, res_vfmav);
  }
 template<typename T, size_t ndim> auto get_OptNpArr_and_vmav(const OptNpArr &arr_,
  const shape_t &dims, const string &name="", size_t nthreads=1)
  {
  if (!arr_) return make_Pyarr_and_vmav<T, ndim>(dims, false, nthreads);
  const auto spec = makeSpec(name);
  auto val = arr_.value();
  MR_assert(isPyarr<T>(val), spec, "incorrect data type");
  MR_assert(dims.size()==size_t(val.ndim()), spec, "dimension mismatch");
  MR_assert(dims.size()==ndim, spec, "dimension mismatch");
  for (size_t i=0; i<dims.size(); ++i)
    MR_assert(dims[i]==size_t(val.shape(int(i))), spec, "dimension mismatch");
  auto res_vmav = to_vmav<T,ndim>(val);
  return std::make_tuple(val, res_vmav);
  }
 template<typename T, size_t ndim> auto get_optional_cmav(const OptCNpArr &arr_,
  const typename mav_info<ndim>::shape_t &dims, const T &defaultval, const string &name="")
  {
  if (!arr_) return cmav<T,ndim>::build_uniform(dims, defaultval);
  const auto spec = makeSpec(name);
  auto val = arr_.value();
  MR_assert(isPyarr<T>(val), spec, "incorrect data type");
  MR_assert(dims.size()==size_t(val.ndim()), spec, "dimension mismatch");
  MR_assert(dims.size()==ndim, spec, "dimension mismatch");
  for (size_t i=0; i<dims.size(); ++i)
    MR_assert(dims[i]==size_t(val.shape(int(i))), spec, "dimension mismatch");
  return to_cmav<T,ndim>(val);
  }

template<typename T> NpArr get_OptNpArr_minshape
  (const OptNpArr &arr_, const shape_t &dims, const string &name="", size_t nthreads=1)
  {
  if (!arr_) return make_Pyarr<T>(dims, false, nthreads);
  const auto spec = makeSpec(name);
  auto val = arr_.value();
  MR_assert(isPyarr<T>(val), spec, "incorrect data type");
  MR_assert(dims.size()==size_t(val.ndim()), spec, "dimension mismatch");
  for (size_t i=0; i<dims.size(); ++i)
    MR_assert(dims[i]<=size_t(val.shape(int(i))), spec, "array shape too small");
  return val;
  }

template<typename T> CNpArr get_OptCNpArr(
  const OptCNpArr &arr_, const shape_t &dims, const string &name="")
  {
  if (!arr_) return CNpArr(make_Pyarr<T>(shape_t(dims.size(), 0)));
  const auto spec = makeSpec(name);
  auto val = arr_.value();
  MR_assert(isPyarr<T>(val), spec, "incorrect data type");
  MR_assert(dims.size()==size_t(val.ndim()), spec, "dimension mismatch");
  for (size_t i=0; i<dims.size(); ++i)
    MR_assert(dims[i]==size_t(val.shape(int(i))), spec, "dimension mismatch");
  return val;
  }

#ifdef DUCC0_USE_NANOBIND
inline py::object normalizeDtype(const py::object &dtype)
  {
  static py::object converter = py::module_::import_("numpy").attr("dtype");
  return converter(dtype);
  }
template<typename T> inline py::object Dtype();
template<> inline py::object Dtype<float>()
  { static auto res = normalizeDtype(py::cast("f4")); return res; }
template<> inline py::object Dtype<double>()
  { static auto res = normalizeDtype(py::cast("f8")); return res; }
template<> inline py::object Dtype<complex<float>>()
  { static auto res = normalizeDtype(py::cast("c8")); return res; }
template<> inline py::object Dtype<complex<double>>()
  { static auto res = normalizeDtype(py::cast("c16")); return res; }
template<typename T> bool isDtype(const py::object &dtype)
  { return Dtype<T>().equal(dtype); }
#else
inline py::dtype normalizeDtype(const py::object &dtype)
  {
  static py::object converter = py::module_::import("numpy").attr("dtype");
  return converter(dtype);
  }
template<typename T> py::dtype Dtype()
  { return py::dtype::of<T>(); }
template<typename T> bool isDtype(const py::dtype &dtype)
  { return Dtype<T>().equal(dtype); }
#endif
}

using detail_pybind::NpArr;
using detail_pybind::NpArrT;
using detail_pybind::OptNpArr;
using detail_pybind::CNpArr;
using detail_pybind::CNpArrT;
using detail_pybind::OptCNpArr;
using detail_pybind::None;
using detail_pybind::isPyarr;
using detail_pybind::make_Pyarr;
using detail_pybind::make_Pyarr_and_vmav;
using detail_pybind::make_noncritical_Pyarr;
using detail_pybind::get_OptNpArr;
using detail_pybind::get_OptNpArr_and_vfmav;
using detail_pybind::get_OptNpArr_and_vmav;
using detail_pybind::get_OptNpArr_minshape;
using detail_pybind::get_OptCNpArr;
using detail_pybind::to_cfmav;
using detail_pybind::to_vfmav;
using detail_pybind::to_cmav;
using detail_pybind::to_cmav_with_optional_leading_dimensions;
using detail_pybind::to_cfmav_with_optional_leading_dimensions;
using detail_pybind::to_vmav;
using detail_pybind::to_vmav_with_optional_leading_dimensions;
using detail_pybind::to_vfmav_with_optional_leading_dimensions;
using detail_pybind::normalizeDtype;
using detail_pybind::isDtype;
using detail_pybind::get_optional_cmav;

}

#endif

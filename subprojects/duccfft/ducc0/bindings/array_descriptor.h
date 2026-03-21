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

/* Copyright (C) 2022-2023 Max-Planck-Society
   Author: Martin Reinecke */

#ifndef DUCC0_ARRAY_DESCRIPTOR_H
#define DUCC0_ARRAY_DESCRIPTOR_H

#include <array>
#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/mav.h"
#include "ducc0/bindings/typecode.h"

namespace ducc0 {

namespace detail_array_descriptor {

using namespace std;

struct ArrayDescriptor
  {
  public:
    static constexpr size_t maxdim=10;

    array<uint64_t, maxdim> shape;
    array<int64_t, maxdim> stride;

    void *data;
    uint8_t ndim;
    uint8_t typecode;
//    uint8_t readonly=0;

    ArrayDescriptor() : data(nullptr), ndim(0), typecode(0)/*, readonly(0)*/ {}
    ArrayDescriptor(void *data_, const vector<size_t> &shape_, uint8_t typecode_)
      : data(data_), ndim(shape_.size()), typecode(typecode_)//, readonly(0)
      {
      size_t str = 1;
      for (int i=ndim-1; i>=0; --i)
        {
        shape[i] = shape_[i];
        stride[i] = str;
        str *= shape[i];
        }
      }
    template<typename T> ArrayDescriptor(const cfmav<T> &in)
      : data(const_cast<T *>(in.data())), ndim(in.ndim()),
        typecode(Typecode<T>::value)//, readonly(1)
      {
      MR_assert(ndim<=maxdim, "dimensionality too high");
      for (size_t i=0; i<ndim; ++i)
        {
        shape[i] = in.shape(ndim-1-i);
        stride[i] = in.stride(ndim-1-i);
        }
      }
    template<typename T, size_t ndim2> ArrayDescriptor(const cmav<T,ndim2> &in)
      : data(const_cast<T *>(in.data())), ndim(ndim2),
        typecode(Typecode<T>::value)//, readonly(1)
      {
      MR_assert(ndim<=maxdim, "dimensionality too high");
      for (size_t i=0; i<ndim; ++i)
        {
        shape[i] = in.shape(ndim-1-i);
        stride[i] = in.stride(ndim-1-i);
        }
      }
    template<typename T> ArrayDescriptor(const vfmav<T> &in)
      : data(in.data()), ndim(in.ndim()),
        typecode(Typecode<T>::value)//, readonly(0)
      {
      MR_assert(ndim<=maxdim, "dimensionality too high");
      for (size_t i=0; i<ndim; ++i)
        {
        shape[i] = in.shape(ndim-1-i);
        stride[i] = in.stride(ndim-1-i);
        }
      }
    template<typename T, size_t ndim2> ArrayDescriptor(const vmav<T,ndim2> &in)
      : data(in.data()), ndim(ndim2),
        typecode(Typecode<T>::value)//, readonly(0)
      {
      MR_assert(ndim<=maxdim, "dimensionality too high");
      for (size_t i=0; i<ndim; ++i)
        {
        shape[i] = in.shape(ndim-1-i);
        stride[i] = in.stride(ndim-1-i);
        }
      }
  private:
    template<bool swapdims, typename T1, typename T2> void copy_data
      (T1 &shp, T2 &str) const
      {
      if constexpr (swapdims)
        for (size_t i=0; i<ndim; ++i)
          {
          shp[i] = shape[ndim-1-i];
          str[i] = stride[ndim-1-i];
          }
      else
        for (size_t i=0; i<ndim; ++i)
          {
          shp[i] = shape[i];
          str[i] = stride[i];
          }
      }

    template<bool swapdims, typename T, size_t ndim2> auto prep1() const
      {
      static_assert(ndim2<=maxdim, "dimensionality too high");
      MR_assert(ndim2==ndim, "dimensionality mismatch");
      MR_assert(Typecode<T>::value==typecode, "data type mismatch");
      typename mav_info<ndim2>::shape_t shp;
      typename mav_info<ndim2>::stride_t str;
      copy_data<swapdims>(shp, str);
      return make_tuple(shp, str);
      }
    template<bool swapdims, typename T> auto prep2() const
      {
      MR_assert(Typecode<T>::value==typecode, "data type mismatch");
      typename fmav_info::shape_t shp(ndim);
      typename fmav_info::stride_t str(ndim);
      copy_data<swapdims>(shp, str);
      return make_tuple(shp, str);
      }

  public:
    template<bool swapdims, typename T, size_t ndim> cmav<T,ndim> to_cmav() const
      {
      auto [shp, str] = prep1<swapdims, T, ndim>();
      return cmav<T, ndim>(reinterpret_cast<const T *>(data), shp, str);
      }
    template<bool swapdims, typename T, typename T2, size_t ndim>
      cmav<T2,ndim> to_cmav_with_typecast() const
      {
      static_assert(sizeof(T)==sizeof(T2), "type size mismatch");
      auto [shp, str] = prep1<swapdims, T, ndim>();
      return cmav<T2, ndim>(reinterpret_cast<const T2 *>(data), shp, str);
      }
    template<bool swapdims, typename T, size_t ndim> vmav<T,ndim> to_vmav() const
      {
//      MR_assert(!readonly, "object is read-only");
      auto [shp, str] = prep1<swapdims, T, ndim>();
      return vmav<T, ndim>(reinterpret_cast<T *>(data), shp, str);
      }
    template<bool swapdims, typename T> cfmav<T> to_cfmav() const
      {
      auto [shp, str] = prep2<swapdims, T>();
      return cfmav<T>(reinterpret_cast<const T *>(data), shp, str);
      }
    template<bool swapdims, typename T> vfmav<T> to_vfmav() const
      {
      auto [shp, str] = prep2<swapdims, T>();
      return vfmav<T>(reinterpret_cast<T *>(data), shp, str);
      }

    template<bool swap_content, typename Tin, typename Tout> vector<Tout> to_vector
      () const
      {
      MR_assert(Typecode<Tin>::value==typecode, "data type mismatch");
      MR_assert(ndim==1, "need 1D array for conversion to vector");
      vector<Tout> res;
      res.reserve(shape[0]);
      auto data_ = reinterpret_cast<const Tin *>(data);
      for (size_t i=0; i<shape[0]; ++i)
        res.push_back(swap_content ? data_[(shape[0]-1-i)*stride[0]]
                                   : data_[i*stride[0]]);
      return res;
      }
    template<bool swap_content, typename Tin, typename Tout> vector<Tout> to_vector_subtract_1
      () const
      {
      static_assert(is_integral<Tin>::value, "need an integral type for this");
      MR_assert(Typecode<Tin>::value==typecode, "data type mismatch");
      MR_assert(ndim==1, "need 1D array for conversion to vector");
      vector<Tout> res;
      res.reserve(shape[0]);
      auto data_ = reinterpret_cast<const Tin *>(data);
      for (size_t i=0; i<shape[0]; ++i)
        {
        Tin val = (swap_content ? data_[(shape[0]-1-i)*stride[0]]
                                : data_[i*stride[0]]) - Tin(1);
        res.push_back(Tout(val));
        }
      return res;
      }
    template<bool swap_content, typename Tout> vector<Tout> to_vector_subtract_1
      () const
      {
      if (isTypecode<int32_t>(typecode))
        return to_vector_subtract_1<swap_content, int32_t, Tout>();
      if (isTypecode<int64_t>(typecode))
        return to_vector_subtract_1<swap_content, int64_t, Tout>();
      if (isTypecode<uint32_t>(typecode))
        return to_vector_subtract_1<swap_content, uint32_t, Tout>();
      if (isTypecode<uint64_t>(typecode))
        return to_vector_subtract_1<swap_content, uint64_t, Tout>();
      MR_fail("no suitable type found");
      }
  };

template<typename T, size_t ndim> cmav<T,ndim> subtract_1(const cmav<T,ndim> &inp)
  {
  vmav<T,ndim> res(inp.shape(), UNINITIALIZED);
  mav_apply([](T &v1, const T &v2){v1=v2-T(1);}, 1, res, inp);
  return res;
  }
}

using detail_array_descriptor::ArrayDescriptor;
using detail_array_descriptor::subtract_1;

}

#endif

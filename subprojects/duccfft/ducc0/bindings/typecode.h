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


#ifndef DUCC0_TYPECODE_H
#define DUCC0_TYPECODE_H

#include <cstddef>
#include <type_traits>
#include <complex>

namespace ducc0 {

namespace detail_typecode {

using namespace std;

// bits [0-3]: the type size in bytes - 1
// bits [4-5]: 0 = float, 1 = signed int, 2 = unsigned int (, 3 = bool)
// bits [6-7]: number of primitive variables in one element - 1
constexpr size_t sizeshift = 0x0;
constexpr size_t typeshift = 0x4;
constexpr size_t numshift = 0x6;
constexpr size_t floatcode = 0x0;
constexpr size_t sintcode = 0x1;
constexpr size_t uintcode = 0x2;

template<typename T> class Typecode
  {
  private:
    static constexpr size_t compute()
      {
      static_assert(!is_same<T,bool>::value, "no bools allowed");
      static_assert(is_integral<T>::value||is_floating_point<T>::value,
        "need integral or floating point type");
      static_assert(sizeof(T)<=16, "type size must be at most 16 bytes");
      if constexpr(is_floating_point<T>::value)
        static_assert(is_signed<T>::value,
          "no support for unsigned floating point types");
      return ((sizeof(T)-1)<<sizeshift) +
             ((floatcode*is_floating_point<T>::value
              +sintcode*(is_integral<T>::value&&is_signed<T>::value)
              +uintcode*(is_integral<T>::value&&(!is_signed<T>::value)))<<typeshift);
      }

  public:
    static constexpr size_t value = compute();
  };

template<typename T> class Typecode<complex<T>>
  {
  private:
    static constexpr size_t compute()
      {
      static_assert(is_floating_point<T>::value, "need a floating point type");
      return Typecode<T>::value + (size_t(1)<<numshift);
      }

  public:
    static constexpr size_t value = compute();
  };

inline size_t typeSize(size_t typecode)
  { return ((typecode&15)+1)*((typecode>>6)+1); }

template<typename T> bool isTypecode(size_t typecode)
  { return typecode==Typecode<T>::value; }

}

using detail_typecode::Typecode;
using detail_typecode::typeSize;
using detail_typecode::isTypecode;

}

#endif

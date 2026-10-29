/** \file ducc0/infra/aligned_array.h
 *
 * \copyright Copyright (C) 2019-2026 Max-Planck-Society
 * \author Martin Reinecke
 */

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

#ifndef DUCC0_ALIGNED_ARRAY_H
#define DUCC0_ALIGNED_ARRAY_H

#include <algorithm>
#include <cstddef>
#include <new>

namespace DUCC0_NAMESPACE {

namespace detail_aligned_array {

using namespace std;

/// Bare bones array class.
/** Mostly useful for uninitialized temporary buffers.
 *  \note Since this class operates on raw memory, it should only be used with
 *        POD types, and even then only with caution! */
template<typename T, size_t alignment=alignof(T)> class array_base
  {
  private:
    static_assert((alignment&(alignment-1))==0, "alignment is not a power of 2");
    static_assert(alignment>=alignof(T), "requested alignment is too small");

    T *p;
    size_t sz;

    static T *ralloc(size_t num)
      {
      if (num==0) return nullptr;
      if constexpr(alignment<=alignof(max_align_t))
        return reinterpret_cast<T *>(::operator new(num*sizeof(T)));
      else
        return reinterpret_cast<T *>(::operator new(num*sizeof(T), std::align_val_t(alignment)));
      }
    static void dealloc(T *ptr)
      {
      if constexpr(alignment<=alignof(max_align_t))
        ::operator delete(ptr);
      else
        ::operator delete(ptr, std::align_val_t(alignment));
      }

  public:
    /// Creates a zero-sized array with no associated memory.
    array_base() : p(nullptr), sz(0) {}
    /// Creates an array with \a n entries.
    /** \note Memory is not initialized! */
    explicit array_base(size_t n) : p(ralloc(n)), sz(n) {}
    array_base(const array_base &) = delete;
    array_base(array_base &&other) noexcept
      : p(other.p), sz(other.sz)
      { other.p=nullptr; other.sz=0; }
    ~array_base() { dealloc(p); }

    array_base &operator=(const array_base &) = delete;
    array_base &operator=(array_base &&other)
      {
      std::swap(p, other.p);
      std::swap(sz, other.sz);
      return *this;
      }

    /// If \a n is different from the current size, resizes the array to hold
    /// \a n elements.
    /** \note No data content is copied, the new array is uninitialized! */
    void realloc(size_t n)
      {
      if (n==sz) return;
      dealloc(p);
      p = ralloc(n);
      sz = n;
      }

    /// Returns a writeable reference to the element at index \a idx.
    T &operator[](size_t idx) { return p[idx]; }
    /// Returns a read-only reference to the element at index \a idx.
    const T &operator[](size_t idx) const { return p[idx]; }

    /// Returns a writeable pointer to the array data.
    T *data() { return p; }
    /// Returns a read-only pointer to the array data.
    const T *data() const { return p; }

    /// Returns the size of the array.
    size_t size() const { return sz; }
  };

template<typename T> using quick_array = array_base<T>;
template<typename T> using aligned_array = array_base<T,64>;

}

using detail_aligned_array::aligned_array;
using detail_aligned_array::quick_array;

}

#endif


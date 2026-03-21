/*
 *  This file is part of the MR utility library.
 *
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

/* Copyright (C) 2020-2021 Max-Planck-Society
   Author: Martin Reinecke */

#include <cstdint>

#include "ducc0/infra/types.h"

using namespace std;

namespace ducc0 {

namespace detail_types {

using namespace std;

class Sizemap: public TypeMapper<size_t>
  {
  protected:
    template<typename T, typename... Ts> void addTypes()
      {
      add<T>(sizeof(T));
      if constexpr (sizeof...(Ts)>0) addTypes<Ts...>();
      }

  public:
    Sizemap()
      {
      addTypes<double, float, int, long, size_t, ptrdiff_t,
               int32_t, int64_t, uint32_t, uint64_t>();
      }
   };

Sizemap sizemap;

size_t typesize(const type_index &idx)
  { return sizemap[idx]; }

}}

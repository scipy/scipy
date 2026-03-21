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

#ifndef DUCC0_TYPES_H
#define DUCC0_TYPES_H

#include <typeindex>
#include <cstddef>
#include <unordered_map>
#include "ducc0/infra/error_handling.h"

namespace ducc0 {

namespace detail_types {

using namespace std;

template<typename T> constexpr inline auto tidx()
  { return type_index(typeid(T)); }

template<typename DT> class TypeMapper
  {
  protected:
    unordered_map<type_index, DT> mapping;

  public:
    template<typename T> void add (const DT &dt)
      { mapping[tidx<T>()] = dt; }

    DT operator[](const type_index &idx) const
      {
      auto res = mapping.find(idx);
      MR_assert(res!=mapping.end(), "type not found");
      return res->second;
      }
   };

size_t typesize(const type_index &idx);

}

using detail_types::tidx;
using detail_types::typesize;
using detail_types::TypeMapper;

}

#endif

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

/*
 *  Utilities for conversion between integer coordinates, Morton, and Peano
 *  indices
 *
 *  Copyright (C) 2015-2023 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef DUCC0_SPACE_FILLING_H
#define DUCC0_SPACE_FILLING_H

#include <cstdint>
#include <array>

// PDEP and PEXT should not be used on Zen1 and Zen2 AMD CPUs, since their
// implementation is really slow. On Intel they provide some performance
// advantage, but I think it's better overall to just switch this off.
// Update: AMD Zen3 and later seem to have fast PDEP/PEXT again. Maybe enable
// this again in the near future?
//#define DUCC0_USE_PDEP_PEXT

#if (!(defined(__BMI2__) && defined(__x86_64__)))
#undef DUCC0_USE_PDEP_PEXT
#endif

#ifdef DUCC0_USE_PDEP_PEXT
#include <x86intrin.h>
#endif

namespace ducc0 {

// Naming conventions
// - 2D/3D: function works on 2D/3D integer coordinates
// - 32/64: the largest input or output data type of this function is
//   a 32bit resp. 64bit unsigned integer. Other parmeters may be the same size
//   or smaller, but the total number of relevant bits can never be larger than
//   32 resp. 64

// conversions between integer coordinates in separate variables ("coord")
// and those stored in different parts of a single integer ("block")
inline uint32_t coord2block2D_32(std::array<uint32_t,2> xy)
  { return (xy[0]&0xffff) | (xy[1]<<16); }
inline std::array<uint32_t,2> block2coord2D_32(uint32_t v)
  { return {v&0xffff,v>>16}; }
inline uint32_t coord2block3D_32(std::array<uint32_t,3> xyz)
  { return (xyz[0]&0x3ff) | ((xyz[1]&0x3ff)<<10) | ((xyz[2]&0x3ff)<<20); }
inline std::array<uint32_t,3> block2coord3D_32(uint32_t v)
  { return {v&0x3ff, (v>>10)&0x3ff, (v>>20)&0x3ff}; }

inline uint64_t coord2block2D_64(std::array<uint64_t,2> xy)
  { return (xy[0]&0xffffffff) | (xy[1]<<32); }
inline std::array<uint64_t,2> block2coord2D_64(uint64_t v)
  { return {v&0xffffffff, v>>32}; }
inline uint64_t coord2block3D_64(std::array<uint64_t,3> xyz)
  { return (xyz[0]&0x1fffff) | ((xyz[1]&0x1fffff)<<21) | ((xyz[2]&0x1fffff)<<42); }
inline std::array<uint64_t,3> block2coord3D_64(uint64_t v)
  { return {v&0x1fffff, (v>>21)&0x1fffff, (v>>42)&0x1fffff}; }

#ifndef DUCC0_USE_PDEP_PEXT

// FIXME: these two are only needed for a very specific case in the Healpix
// neighbors() function. Perhaps we can solve this differently?
uint32_t spread_bits_2D_32 (uint32_t v);
uint64_t spread_bits_2D_64 (uint64_t v);

uint32_t block2morton2D_32 (uint32_t v);
uint32_t coord2morton2D_32 (std::array<uint32_t,2> xy);
uint32_t morton2block2D_32 (uint32_t v);
std::array<uint32_t,2> morton2coord2D_32 (uint32_t v);
uint64_t block2morton2D_64 (uint64_t v);
uint64_t coord2morton2D_64 (std::array<uint64_t,2> xy);
uint64_t morton2block2D_64 (uint64_t v);
std::array<uint64_t,2> morton2coord2D_64 (uint64_t v);

uint32_t block2morton3D_32 (uint32_t v);
uint32_t coord2morton3D_32 (std::array<uint32_t,3> xyz);
uint32_t morton2block3D_32 (uint32_t v);
std::array<uint32_t,3> morton2coord3D_32 (uint32_t v);
uint64_t block2morton3D_64 (uint64_t v);
uint64_t coord2morton3D_64 (std::array<uint64_t,3> xyz);
uint64_t morton2block3D_64 (uint64_t v);
std::array<uint64_t,3> morton2coord3D_64 (uint64_t v);

#else

inline uint32_t spread_bits_2D_32 (uint32_t v)
  { return _pdep_u32(v,0x55555555u); }
inline uint64_t spread_bits_2D_64 (uint64_t v)
  { return _pdep_u64(v,0x5555555555555555u); }

inline uint32_t block2morton2D_32 (uint32_t v)
  { return _pdep_u32(v,0x55555555u)|_pdep_u32(v>>16,0xaaaaaaaau); }
inline uint32_t coord2morton2D_32 (std::array<uint32_t,2> xy)
  { return _pdep_u32(xy[0],0x55555555u)|_pdep_u32(xy[1],0xaaaaaaaau); }
inline uint32_t morton2block2D_32 (uint32_t v)
  { return _pext_u32(v,0x55555555u)|(_pext_u32(v,0xaaaaaaaau)<<16); }
inline std::array<uint32_t,2> morton2coord2D_32 (uint32_t v)
  { return {_pext_u32(v,0x55555555u), _pext_u32(v,0xaaaaaaaau)}; }
inline uint64_t block2morton2D_64 (uint64_t v)
  {
  return _pdep_u64(v,0x5555555555555555u)
        |_pdep_u64(v>>32,0xaaaaaaaaaaaaaaaau);
  }
inline uint64_t coord2morton2D_64 (std::array<uint64_t,2> xy)
  { return _pdep_u64(xy[0],0x5555555555555555u)|
           _pdep_u64(xy[1],0xaaaaaaaaaaaaaaaau); }
inline uint64_t morton2block2D_64 (uint64_t v)
  {
  return _pext_u64(v,0x5555555555555555u)
       |(_pext_u64(v,0xaaaaaaaaaaaaaaaau)<<32);
  }
inline std::array<uint64_t,2> morton2coord2D_64 (uint64_t v)
  {
  return {_pext_u64(v,0x5555555555555555u),
          _pext_u64(v,0xaaaaaaaaaaaaaaaau)};
  }

inline uint32_t block2morton3D_32 (uint32_t v)
  {
  return _pdep_u32(v    ,0x09249249u)
        |_pdep_u32(v>>10,0x12492492u)
        |_pdep_u32(v>>20,0x24924924u);
  }
inline uint32_t coord2morton3D_32 (std::array<uint32_t,3> xyz)
  {
  return _pdep_u32(xyz[0],0x09249249u)
        |_pdep_u32(xyz[1],0x12492492u)
        |_pdep_u32(xyz[2],0x24924924u);
  }
inline uint32_t morton2block3D_32 (uint32_t v)
  {
  return _pext_u32(v,0x9249249u)
       |(_pext_u32(v,0x12492492u)<<10)
       |(_pext_u32(v,0x24924924u)<<20);
  }
inline std::array<uint32_t,3> morton2coord3D_32 (uint32_t v)
  {
  return {_pext_u32(v,0x09249249u),
          _pext_u32(v,0x12492492u),
          _pext_u32(v,0x24924924u)};
  }
inline uint64_t block2morton3D_64 (uint64_t v)
  {
  return _pdep_u64(v    ,0x1249249249249249u)
        |_pdep_u64(v>>21,0x2492492492492492u)
        |_pdep_u64(v>>42,0x4924924924924924u);
  }
inline uint64_t coord2morton3D_64 (std::array<uint64_t,3> xyz)
  {
  return _pdep_u64(xyz[0],0x1249249249249249u)
        |_pdep_u64(xyz[1],0x2492492492492492u)
        |_pdep_u64(xyz[2],0x4924924924924924u);
  }
inline uint64_t morton2block3D_64 (uint64_t v)
  {
  return _pext_u64(v,0x1249249249249249u)
       |(_pext_u64(v,0x2492492492492492u)<<21)
       |(_pext_u64(v,0x4924924924924924u)<<42);
  }
inline std::array<uint64_t,3> morton2coord3D_64 (uint64_t v)
  {
  return {_pext_u64(v,0x1249249249249249u),
          _pext_u64(v,0x2492492492492492u),
          _pext_u64(v,0x4924924924924924u)};
  }
#endif

uint32_t morton2peano2D_32(uint32_t v, unsigned bits);
uint32_t peano2morton2D_32(uint32_t v, unsigned bits);

uint64_t morton2peano2D_64(uint64_t v, unsigned bits);
uint64_t peano2morton2D_64(uint64_t v, unsigned bits);

uint32_t morton2peano3D_32(uint32_t v, unsigned bits);
uint32_t peano2morton3D_32(uint32_t v, unsigned bits);

uint64_t morton2peano3D_64(uint64_t v, unsigned bits);
uint64_t peano2morton3D_64(uint64_t v, unsigned bits);

}

#endif

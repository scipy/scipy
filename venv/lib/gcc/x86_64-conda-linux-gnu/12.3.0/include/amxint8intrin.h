/* Copyright (C) 2020-2022 Free Software Foundation, Inc.

   This file is part of GCC.

   GCC is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   GCC is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   Under Section 7 of GPL version 3, you are granted additional
   permissions described in the GCC Runtime Library Exception, version
   3.1, as published by the Free Software Foundation.

   You should have received a copy of the GNU General Public License and
   a copy of the GCC Runtime Library Exception along with this program;
   see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
   <http://www.gnu.org/licenses/>.  */

#if !defined _IMMINTRIN_H_INCLUDED
#error "Never use <amxint8intrin.h> directly; include <immintrin.h> instead."
#endif

#ifndef _AMXINT8INTRIN_H_INCLUDED
#define _AMXINT8INTRIN_H_INCLUDED

#if !defined(__AMX_INT8__)
#pragma GCC push_options
#pragma GCC target("amx-int8")
#define __DISABLE_AMX_INT8__
#endif /* __AMX_INT8__ */

#if defined(__x86_64__)
#define _tile_int8_dp_internal(name,dst,src1,src2)					\
  __asm__ volatile							\
  ("{"#name"\t%%tmm"#src2", %%tmm"#src1", %%tmm"#dst"|"#name"\t%%tmm"#dst", %%tmm"#src1", %%tmm"#src2"}" ::)

#define _tile_dpbssd(dst,src1,src2)					\
  _tile_int8_dp_internal (tdpbssd, dst, src1, src2)

#define _tile_dpbsud(dst,src1,src2)					\
  _tile_int8_dp_internal (tdpbsud, dst, src1, src2)

#define _tile_dpbusd(dst,src1,src2)					\
  _tile_int8_dp_internal (tdpbusd, dst, src1, src2)

#define _tile_dpbuud(dst,src1,src2)					\
  _tile_int8_dp_internal (tdpbuud, dst, src1, src2)

#endif

#ifdef __DISABLE_AMX_INT8__
#undef __DISABLE_AMX_INT8__
#pragma GCC pop_options
#endif /* __DISABLE_AMX_INT8__ */

#endif /* _AMXINT8INTRIN_H_INCLUDED */

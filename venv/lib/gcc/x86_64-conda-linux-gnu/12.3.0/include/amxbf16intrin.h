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
#error "Never use <amxbf16intrin.h> directly; include <immintrin.h> instead."
#endif

#ifndef _AMXBF16INTRIN_H_INCLUDED
#define _AMXBF16INTRIN_H_INCLUDED

#if !defined(__AMX_BF16__)
#pragma GCC push_options
#pragma GCC target("amx-bf16")
#define __DISABLE_AMX_BF16__
#endif /* __AMX_BF16__ */

#if defined(__x86_64__)
#define _tile_dpbf16ps_internal(dst,src1,src2)					\
  __asm__ volatile\
  ("{tdpbf16ps\t%%tmm"#src2", %%tmm"#src1", %%tmm"#dst"|tdpbf16ps\t%%tmm"#dst", %%tmm"#src1", %%tmm"#src2"}" ::)

#define _tile_dpbf16ps(dst,src1,src2)					\
  _tile_dpbf16ps_internal (dst, src1, src2)

#endif

#ifdef __DISABLE_AMX_BF16__
#undef __DISABLE_AMX_BF16__
#pragma GCC pop_options
#endif /* __DISABLE_AMX_BF16__ */

#endif /* _AMXBF16INTRIN_H_INCLUDED */

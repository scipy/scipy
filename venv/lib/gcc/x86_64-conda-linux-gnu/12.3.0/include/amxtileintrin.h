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
#error "Never use <amxtileintrin.h> directly; include <immintrin.h> instead."
#endif

#ifndef _AMXTILEINTRIN_H_INCLUDED
#define _AMXTILEINTRIN_H_INCLUDED

#if !defined(__AMX_TILE__)
#pragma GCC push_options
#pragma GCC target("amx-tile")
#define __DISABLE_AMX_TILE__
#endif /* __AMX_TILE__ */

#if defined(__x86_64__)
extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_tile_loadconfig (const void *__config)
{
  __asm__ volatile ("ldtilecfg\t%X0" :: "m" (*((const void **)__config)));
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_tile_storeconfig (void *__config)
{
  __asm__ volatile ("sttilecfg\t%X0" : "=m" (*((void **)__config)));
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_tile_release (void)
{
  __asm__ volatile ("tilerelease" ::);
}

#define _tile_loadd(dst,base,stride)		\
  _tile_loadd_internal (dst, base, stride)

#define _tile_loadd_internal(dst,base,stride)				\
  __asm__ volatile							\
  ("{tileloadd\t(%0,%1,1), %%tmm"#dst"|tileloadd\t%%tmm"#dst", [%0+%1*1]}" \
   :: "r" ((const void*) (base)), "r" ((__PTRDIFF_TYPE__) (stride)))

#define _tile_stream_loadd(dst,base,stride)		\
  _tile_stream_loadd_internal (dst, base, stride)

#define _tile_stream_loadd_internal(dst,base,stride)			\
  __asm__ volatile							\
  ("{tileloaddt1\t(%0,%1,1), %%tmm"#dst"|tileloaddt1\t%%tmm"#dst", [%0+%1*1]}" \
   :: "r" ((const void*) (base)), "r" ((__PTRDIFF_TYPE__) (stride)))

#define _tile_stored(dst,base,stride)		\
  _tile_stored_internal (dst, base, stride)

#define _tile_stored_internal(src,base,stride)				\
  __asm__ volatile							\
  ("{tilestored\t%%tmm"#src", (%0,%1,1)|tilestored\t[%0+%1*1], %%tmm"#src"}" \
   :: "r" ((void*) (base)), "r" ((__PTRDIFF_TYPE__) (stride)) \
   : "memory")

#define _tile_zero(dst)				\
  _tile_zero_internal (dst)

#define _tile_zero_internal(dst)		\
  __asm__ volatile				\
  ("tilezero\t%%tmm"#dst ::)

#endif

#ifdef __DISABLE_AMX_TILE__
#undef __DISABLE_AMX_TILE__
#pragma GCC pop_options
#endif /* __DISABLE_AMX_TILE__ */

#endif /* _AMXTILEINTRIN_H_INCLUDED */

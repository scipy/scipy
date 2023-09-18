/* Copyright (C) 2021-2022 Free Software Foundation, Inc.

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

#ifndef _MWAITINTRIN_H_INCLUDED
#define _MWAITINTRIN_H_INCLUDED

#ifndef __MWAIT__
#pragma GCC push_options
#pragma GCC target("mwait")
#define __DISABLE_MWAIT__
#endif /* __MWAIT__ */

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_monitor (void const * __P, unsigned int __E, unsigned int __H)
{
  __builtin_ia32_monitor (__P, __E, __H);
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mwait (unsigned int __E, unsigned int __H)
{
  __builtin_ia32_mwait (__E, __H);
}

#ifdef __DISABLE_MWAIT__
#undef __DISABLE_MWAIT__
#pragma GCC pop_options
#endif /* __DISABLE_MWAIT__ */

#endif /* _MWAITINTRIN_H_INCLUDED */

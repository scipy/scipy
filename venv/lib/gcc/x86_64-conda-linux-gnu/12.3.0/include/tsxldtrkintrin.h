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

#ifndef _X86GPRINTRIN_H_INCLUDED
# error "Never use <tsxldtrkintrin.h> directly; include <x86gprintrin.h> instead."
#endif

#ifndef _TSXLDTRKINTRIN_H_INCLUDED
#define _TSXLDTRKINTRIN_H_INCLUDED

#if !defined(__TSXLDTRK__)
#pragma GCC push_options
#pragma GCC target("tsxldtrk")
#define __DISABLE_TSXLDTRK__
#endif /* __TSXLDTRK__ */

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_xsusldtrk (void)
{
  __builtin_ia32_xsusldtrk ();
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_xresldtrk (void)
{
  __builtin_ia32_xresldtrk ();
}

#ifdef __DISABLE_TSXLDTRK__
#undef __DISABLE_TSXLDTRK__
#pragma GCC pop_options
#endif /* __DISABLE_TSXLDTRK__ */

#endif /* _TSXLDTRKINTRIN_H_INCLUDED */

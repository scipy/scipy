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
# error "Never use <uintrintrin.h> directly; include <x86gprintrin.h> instead."
#endif

#ifndef _UINTRNTRIN_H_INCLUDED
#define _UINTRNTRIN_H_INCLUDED

#ifdef __x86_64__

#ifndef __UINTR__
#pragma GCC push_options
#pragma GCC target ("uintr")
#define __DISABLE_UINTR__
#endif /* __UINTR__ */

struct __uintr_frame
{
  /* RIP of the interrupted user process.  */
  unsigned long long rip;
  /* RFLAGS of the interrupted user process.  */
  unsigned long long rflags;
  /* RSP of the interrupted user process.  */
  unsigned long long rsp;
};

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_clui (void)
{
  __builtin_ia32_clui ();
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_stui (void)
{
  __builtin_ia32_stui ();
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_senduipi (unsigned long long __R)
{
  __builtin_ia32_senduipi (__R);
}

extern __inline unsigned char
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_testui (void)
{
  return __builtin_ia32_testui ();
}

#ifdef __DISABLE_UINTR__
#undef __DISABLE_UINTR__
#pragma GCC pop_options
#endif /* __DISABLE_UINTR__ */

#endif

#endif /* _UINTRNTRIN_H_INCLUDED.  */

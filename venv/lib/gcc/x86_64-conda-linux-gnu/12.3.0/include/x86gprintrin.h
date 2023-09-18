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
#define _X86GPRINTRIN_H_INCLUDED

#if !defined _SOFT_FLOAT || defined __MMX__ || defined __SSE__
#pragma GCC push_options
#pragma GCC target("general-regs-only")
#define __DISABLE_GENERAL_REGS_ONLY__
#endif

#include <ia32intrin.h>

#ifndef __iamcu__

#include <stddef.h>

#include <adxintrin.h>

#include <bmiintrin.h>

#include <bmi2intrin.h>

#include <cetintrin.h>

#include <cldemoteintrin.h>

#include <clflushoptintrin.h>

#include <clwbintrin.h>

#include <clzerointrin.h>

#include <enqcmdintrin.h>

#include <fxsrintrin.h>

#include <lzcntintrin.h>

#include <lwpintrin.h>

#include <movdirintrin.h>

#include <mwaitintrin.h>

#include <mwaitxintrin.h>

#include <pconfigintrin.h>

#include <popcntintrin.h>

#include <pkuintrin.h>

#include <rdseedintrin.h>

#include <rtmintrin.h>

#include <serializeintrin.h>

#include <sgxintrin.h>

#include <tbmintrin.h>

#include <tsxldtrkintrin.h>

#include <uintrintrin.h>

#include <waitpkgintrin.h>

#include <wbnoinvdintrin.h>

#include <xsaveintrin.h>

#include <xsavecintrin.h>

#include <xsaveoptintrin.h>

#include <xsavesintrin.h>

#include <xtestintrin.h>

#include <hresetintrin.h>

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_wbinvd (void)
{
  __builtin_ia32_wbinvd ();
}

#ifndef __RDRND__
#pragma GCC push_options
#pragma GCC target("rdrnd")
#define __DISABLE_RDRND__
#endif /* __RDRND__ */
extern __inline int
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_rdrand16_step (unsigned short *__P)
{
  return __builtin_ia32_rdrand16_step (__P);
}

extern __inline int
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_rdrand32_step (unsigned int *__P)
{
  return __builtin_ia32_rdrand32_step (__P);
}
#ifdef __DISABLE_RDRND__
#undef __DISABLE_RDRND__
#pragma GCC pop_options
#endif /* __DISABLE_RDRND__ */

#ifndef __RDPID__
#pragma GCC push_options
#pragma GCC target("rdpid")
#define __DISABLE_RDPID__
#endif /* __RDPID__ */
extern __inline unsigned int
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_rdpid_u32 (void)
{
  return __builtin_ia32_rdpid ();
}
#ifdef __DISABLE_RDPID__
#undef __DISABLE_RDPID__
#pragma GCC pop_options
#endif /* __DISABLE_RDPID__ */

#ifdef  __x86_64__

#ifndef __FSGSBASE__
#pragma GCC push_options
#pragma GCC target("fsgsbase")
#define __DISABLE_FSGSBASE__
#endif /* __FSGSBASE__ */
extern __inline unsigned int
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_readfsbase_u32 (void)
{
  return __builtin_ia32_rdfsbase32 ();
}

extern __inline unsigned long long
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_readfsbase_u64 (void)
{
  return __builtin_ia32_rdfsbase64 ();
}

extern __inline unsigned int
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_readgsbase_u32 (void)
{
  return __builtin_ia32_rdgsbase32 ();
}

extern __inline unsigned long long
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_readgsbase_u64 (void)
{
  return __builtin_ia32_rdgsbase64 ();
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_writefsbase_u32 (unsigned int __B)
{
  __builtin_ia32_wrfsbase32 (__B);
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_writefsbase_u64 (unsigned long long __B)
{
  __builtin_ia32_wrfsbase64 (__B);
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_writegsbase_u32 (unsigned int __B)
{
  __builtin_ia32_wrgsbase32 (__B);
}

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_writegsbase_u64 (unsigned long long __B)
{
  __builtin_ia32_wrgsbase64 (__B);
}
#ifdef __DISABLE_FSGSBASE__
#undef __DISABLE_FSGSBASE__
#pragma GCC pop_options
#endif /* __DISABLE_FSGSBASE__ */

#ifndef __RDRND__
#pragma GCC push_options
#pragma GCC target("rdrnd")
#define __DISABLE_RDRND__
#endif /* __RDRND__ */
extern __inline int
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_rdrand64_step (unsigned long long *__P)
{
  return __builtin_ia32_rdrand64_step (__P);
}
#ifdef __DISABLE_RDRND__
#undef __DISABLE_RDRND__
#pragma GCC pop_options
#endif /* __DISABLE_RDRND__ */

#endif /* __x86_64__  */

#ifndef __PTWRITE__
#pragma GCC push_options
#pragma GCC target("ptwrite")
#define __DISABLE_PTWRITE__
#endif

#ifdef __x86_64__
extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_ptwrite64 (unsigned long long __B)
{
  __builtin_ia32_ptwrite64 (__B);
}
#endif /* __x86_64__ */

extern __inline void
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
_ptwrite32 (unsigned __B)
{
  __builtin_ia32_ptwrite32 (__B);
}
#ifdef __DISABLE_PTWRITE__
#undef __DISABLE_PTWRITE__
#pragma GCC pop_options
#endif /* __DISABLE_PTWRITE__ */

#endif /* __iamcu__ */

#ifdef __DISABLE_GENERAL_REGS_ONLY__
#undef __DISABLE_GENERAL_REGS_ONLY__
#pragma GCC pop_options
#endif /* __DISABLE_GENERAL_REGS_ONLY__ */

#endif /* _X86GPRINTRIN_H_INCLUDED.  */

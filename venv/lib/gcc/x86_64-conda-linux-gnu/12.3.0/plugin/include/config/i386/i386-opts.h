/* Definitions for option handling for IA-32.
   Copyright (C) 1988-2022 Free Software Foundation, Inc.

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

#ifndef I386_OPTS_H
#define I386_OPTS_H

/* Algorithm to expand string function with.  */
enum stringop_alg
{
#undef DEF_ALG
#define DEF_ALG(alg, name) alg, 

#include "stringop.def"
last_alg

#undef DEF_ALG
};

/* Available call abi.  */
enum calling_abi
{
  SYSV_ABI = 0,
  MS_ABI = 1
};

enum fpmath_unit
{
  FPMATH_387 = 1,
  FPMATH_SSE = 2
};

enum tls_dialect
{
  TLS_DIALECT_GNU,
  TLS_DIALECT_GNU2,
  TLS_DIALECT_SUN
};

enum cmodel {
  CM_32,	/* The traditional 32-bit ABI.  */
  CM_SMALL,	/* Assumes all code and data fits in the low 31 bits.  */
  CM_KERNEL,	/* Assumes all code and data fits in the high 31 bits.  */
  CM_MEDIUM,	/* Assumes code fits in the low 31 bits; data unlimited.  */
  CM_LARGE,	/* No assumptions.  */
  CM_SMALL_PIC,	/* Assumes code+data+got/plt fits in a 31 bit region.  */
  CM_MEDIUM_PIC,/* Assumes code+got/plt fits in a 31 bit region.  */
  CM_LARGE_PIC	/* No assumptions.  */
};

enum pmode {
  PMODE_SI,	/* Pmode == SImode. */
  PMODE_DI 	/* Pmode == DImode. */
};

enum ix86_align_data {
  ix86_align_data_type_compat,
  ix86_align_data_type_abi,
  ix86_align_data_type_cacheline
};

enum asm_dialect {
  ASM_ATT,
  ASM_INTEL
};

enum ix86_veclibabi {
  ix86_veclibabi_type_none,
  ix86_veclibabi_type_svml,
  ix86_veclibabi_type_acml
};

enum stack_protector_guard {
  SSP_TLS,      /* per-thread canary in TLS block */
  SSP_GLOBAL    /* global canary */
};

enum prefer_vector_width {
    PVW_NONE,
    PVW_AVX128,
    PVW_AVX256,
    PVW_AVX512
};

/* This is used to mitigate variant #2 of the speculative execution
   vulnerabilities on x86 processors identified by CVE-2017-5715, aka
   Spectre.  They convert indirect branches and function returns to
   call and return thunks to avoid speculative execution via indirect
   call, jmp and ret.  */
enum indirect_branch {
  indirect_branch_unset = 0,
  indirect_branch_keep,
  indirect_branch_thunk,
  indirect_branch_thunk_inline,
  indirect_branch_thunk_extern
};

enum instrument_return {
  instrument_return_none = 0,
  instrument_return_call,
  instrument_return_nop5
};

enum harden_sls {
  harden_sls_none = 0,
  harden_sls_return = 1 << 0,
  harden_sls_indirect_jmp = 1 << 1,
  harden_sls_all = harden_sls_return | harden_sls_indirect_jmp
};

#endif

/* Get CPU type and Features for x86 processors.
   Copyright (C) 2012-2022 Free Software Foundation, Inc.
   Contributed by Sriraman Tallam (tmsriram@google.com)

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3, or (at your option) any later
version.

GCC is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

Under Section 7 of GPL version 3, you are granted additional
permissions described in the GCC Runtime Library Exception, version
3.1, as published by the Free Software Foundation.

You should have received a copy of the GNU General Public License and
a copy of the GCC Runtime Library Exception along with this program;
see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
<http://www.gnu.org/licenses/>.  */

/* Processor Vendor and Models. */

enum processor_vendor
{
  VENDOR_INTEL = 1,
  VENDOR_AMD,
  VENDOR_OTHER,
  VENDOR_CENTAUR,
  VENDOR_CYRIX,
  VENDOR_NSC,

  /* Maximum values must be at the end of this enum.  */
  VENDOR_MAX,
  BUILTIN_VENDOR_MAX = VENDOR_OTHER
};

/* Any new types or subtypes have to be inserted at the end. */

enum processor_types
{
  INTEL_BONNELL = 1,
  INTEL_CORE2,
  INTEL_COREI7,
  AMDFAM10H,
  AMDFAM15H,
  INTEL_SILVERMONT,
  INTEL_KNL,
  AMD_BTVER1,
  AMD_BTVER2,
  AMDFAM17H,
  INTEL_KNM,
  INTEL_GOLDMONT,
  INTEL_GOLDMONT_PLUS,
  INTEL_TREMONT,
  AMDFAM19H,
  CPU_TYPE_MAX,
  BUILTIN_CPU_TYPE_MAX = CPU_TYPE_MAX
};

enum processor_subtypes
{
  INTEL_COREI7_NEHALEM = 1,
  INTEL_COREI7_WESTMERE,
  INTEL_COREI7_SANDYBRIDGE,
  AMDFAM10H_BARCELONA,
  AMDFAM10H_SHANGHAI,
  AMDFAM10H_ISTANBUL,
  AMDFAM15H_BDVER1,
  AMDFAM15H_BDVER2,
  AMDFAM15H_BDVER3,
  AMDFAM15H_BDVER4,
  AMDFAM17H_ZNVER1,
  INTEL_COREI7_IVYBRIDGE,
  INTEL_COREI7_HASWELL,
  INTEL_COREI7_BROADWELL,
  INTEL_COREI7_SKYLAKE,
  INTEL_COREI7_SKYLAKE_AVX512,
  INTEL_COREI7_CANNONLAKE,
  INTEL_COREI7_ICELAKE_CLIENT,
  INTEL_COREI7_ICELAKE_SERVER,
  AMDFAM17H_ZNVER2,
  INTEL_COREI7_CASCADELAKE,
  INTEL_COREI7_TIGERLAKE,
  INTEL_COREI7_COOPERLAKE,
  INTEL_COREI7_SAPPHIRERAPIDS,
  INTEL_COREI7_ALDERLAKE,
  AMDFAM19H_ZNVER3,
  INTEL_COREI7_ROCKETLAKE,
  AMDFAM19H_ZNVER4,
  CPU_SUBTYPE_MAX
};

/* Priority of i386 features, greater value is higher priority.   This is
   used to decide the order in which function dispatch must happen.  For
   instance, a version specialized for SSE4.2 should be checked for dispatch
   before a version for SSE3, as SSE4.2 implies SSE3.  */
enum feature_priority
{
  P_NONE = 0,
  P_MMX,
  P_SSE,
  P_SSE2,
  P_X86_64_BASELINE,
  P_SSE3,
  P_SSSE3,
  P_PROC_SSSE3,
  P_SSE4_A,
  P_PROC_SSE4_A,
  P_SSE4_1,
  P_SSE4_2,
  P_PROC_SSE4_2,
  P_POPCNT,
  P_X86_64_V2,
  P_AES,
  P_PCLMUL,
  P_AVX,
  P_PROC_AVX,
  P_BMI,
  P_PROC_BMI,
  P_FMA4,
  P_XOP,
  P_PROC_XOP,
  P_FMA,
  P_PROC_FMA,
  P_BMI2,
  P_AVX2,
  P_PROC_AVX2,
  P_X86_64_V3,
  P_AVX512F,
  P_PROC_AVX512F,
  P_X86_64_V4,
  P_PROC_DYNAMIC
};

/* ISA Features supported. New features have to be inserted at the end.  */

enum processor_features
{
  FEATURE_CMOV = 0,
  FEATURE_MMX,
  FEATURE_POPCNT,
  FEATURE_SSE,
  FEATURE_SSE2,
  FEATURE_SSE3,
  FEATURE_SSSE3,
  FEATURE_SSE4_1,
  FEATURE_SSE4_2,
  FEATURE_AVX,
  FEATURE_AVX2,
  FEATURE_SSE4_A,
  FEATURE_FMA4,
  FEATURE_XOP,
  FEATURE_FMA,
  FEATURE_AVX512F,
  FEATURE_BMI,
  FEATURE_BMI2,
  FEATURE_AES,
  FEATURE_PCLMUL,
  FEATURE_AVX512VL,
  FEATURE_AVX512BW,
  FEATURE_AVX512DQ,
  FEATURE_AVX512CD,
  FEATURE_AVX512ER,
  FEATURE_AVX512PF,
  FEATURE_AVX512VBMI,
  FEATURE_AVX512IFMA,
  FEATURE_AVX5124VNNIW,
  FEATURE_AVX5124FMAPS,
  FEATURE_AVX512VPOPCNTDQ,
  FEATURE_AVX512VBMI2,
  FEATURE_GFNI,
  FEATURE_VPCLMULQDQ,
  FEATURE_AVX512VNNI,
  FEATURE_AVX512BITALG,
  FEATURE_AVX512BF16,
  FEATURE_AVX512VP2INTERSECT,
  FEATURE_3DNOW,
  FEATURE_3DNOWP,
  FEATURE_ADX,
  FEATURE_ABM,
  FEATURE_CLDEMOTE,
  FEATURE_CLFLUSHOPT,
  FEATURE_CLWB,
  FEATURE_CLZERO,
  FEATURE_CMPXCHG16B,
  FEATURE_CMPXCHG8B,
  FEATURE_ENQCMD,
  FEATURE_F16C,
  FEATURE_FSGSBASE,
  FEATURE_FXSAVE,
  FEATURE_HLE,
  FEATURE_IBT,
  FEATURE_LAHF_LM,
  FEATURE_LM,
  FEATURE_LWP,
  FEATURE_LZCNT,
  FEATURE_MOVBE,
  FEATURE_MOVDIR64B,
  FEATURE_MOVDIRI,
  FEATURE_MWAITX,
  FEATURE_OSXSAVE,
  FEATURE_PCONFIG,
  FEATURE_PKU,
  FEATURE_PREFETCHWT1,
  FEATURE_PRFCHW,
  FEATURE_PTWRITE,
  FEATURE_RDPID,
  FEATURE_RDRND,
  FEATURE_RDSEED,
  FEATURE_RTM,
  FEATURE_SERIALIZE,
  FEATURE_SGX,
  FEATURE_SHA,
  FEATURE_SHSTK,
  FEATURE_TBM,
  FEATURE_TSXLDTRK,
  FEATURE_VAES,
  FEATURE_WAITPKG,
  FEATURE_WBNOINVD,
  FEATURE_XSAVE,
  FEATURE_XSAVEC,
  FEATURE_XSAVEOPT,
  FEATURE_XSAVES,
  FEATURE_AMX_TILE,
  FEATURE_AMX_INT8,
  FEATURE_AMX_BF16,
  FEATURE_UINTR,
  FEATURE_HRESET,
  FEATURE_KL,
  FEATURE_AESKLE,
  FEATURE_WIDEKL,
  FEATURE_AVXVNNI,
  FEATURE_AVX512FP16,
  FEATURE_X86_64_BASELINE,
  FEATURE_X86_64_V2,
  FEATURE_X86_64_V3,
  FEATURE_X86_64_V4,
  CPU_FEATURE_MAX
};

/* Size of __cpu_features2 array in libgcc/config/i386/cpuinfo.c.  */
#define SIZE_OF_CPU_FEATURES ((CPU_FEATURE_MAX - 1) / 32)

/* These are the values for vendor types, cpu types and subtypes.  Cpu
   types and subtypes should be subtracted by the corresponding start
   value.  */

#define M_CPU_TYPE_START (BUILTIN_VENDOR_MAX)
#define M_CPU_SUBTYPE_START \
  (M_CPU_TYPE_START + BUILTIN_CPU_TYPE_MAX)
#define M_VENDOR(a) (a)
#define M_CPU_TYPE(a) (M_CPU_TYPE_START + a)
#define M_CPU_SUBTYPE(a) (M_CPU_SUBTYPE_START + a)

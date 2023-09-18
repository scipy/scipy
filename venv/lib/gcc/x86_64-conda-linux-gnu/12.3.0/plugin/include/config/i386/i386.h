/* Definitions of target machine for GCC for IA-32.
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

/* The purpose of this file is to define the characteristics of the i386,
   independent of assembler syntax or operating system.

   Three other files build on this one to describe a specific assembler syntax:
   bsd386.h, att386.h, and sun386.h.

   The actual tm.h file for a particular system should include
   this file, and then the file for the appropriate assembler syntax.

   Many macros that specify assembler syntax are omitted entirely from
   this file because they really belong in the files for particular
   assemblers.  These include RP, IP, LPREFIX, PUT_OP_SIZE, USE_STAR,
   ADDR_BEG, ADDR_END, PRINT_IREG, PRINT_SCALE, PRINT_B_I_S, and many
   that start with ASM_ or end in ASM_OP.  */

/* Redefines for option macros.  */

#define TARGET_CMPXCHG16B	TARGET_CX16
#define TARGET_CMPXCHG16B_P(x)	TARGET_CX16_P(x)

#define TARGET_LP64		TARGET_ABI_64
#define TARGET_LP64_P(x)	TARGET_ABI_64_P(x)
#define TARGET_X32		TARGET_ABI_X32
#define TARGET_X32_P(x)		TARGET_ABI_X32_P(x)
#define TARGET_16BIT		TARGET_CODE16
#define TARGET_16BIT_P(x)	TARGET_CODE16_P(x)

#define TARGET_MMX_WITH_SSE	(TARGET_64BIT && TARGET_SSE2)

#include "config/vxworks-dummy.h"

#include "config/i386/i386-opts.h"

#define MAX_STRINGOP_ALGS 4

/* Specify what algorithm to use for stringops on known size.
   When size is unknown, the UNKNOWN_SIZE alg is used.  When size is
   known at compile time or estimated via feedback, the SIZE array
   is walked in order until MAX is greater then the estimate (or -1
   means infinity).  Corresponding ALG is used then.
   When NOALIGN is true the code guaranting the alignment of the memory
   block is skipped.

   For example initializer:
    {{256, loop}, {-1, rep_prefix_4_byte}}
   will use loop for blocks smaller or equal to 256 bytes, rep prefix will
   be used otherwise.  */
struct stringop_algs
{
  const enum stringop_alg unknown_size;
  const struct stringop_strategy {
    /* Several older compilers delete the default constructor because of the
       const entries (see PR100246).  Manually specifying a CTOR works around
       this issue.  Since this header is used by code compiled with the C
       compiler we must guard the addition.  */
#ifdef __cplusplus
    constexpr
    stringop_strategy (int _max = -1, enum stringop_alg _alg = libcall,
		       int _noalign = false)
      : max (_max), alg (_alg), noalign (_noalign) {}
#endif
    const int max;
    const enum stringop_alg alg;
    int noalign;
  } size [MAX_STRINGOP_ALGS];
};

/* Analog of COSTS_N_INSNS when optimizing for size.  */
#ifndef COSTS_N_BYTES
#define COSTS_N_BYTES(N) ((N) * 2)
#endif

/* Define the specific costs for a given cpu.  NB: hard_register is used
   by TARGET_REGISTER_MOVE_COST and TARGET_MEMORY_MOVE_COST to compute
   hard register move costs by register allocator.  Relative costs of
   pseudo register load and store versus pseudo register moves in RTL
   expressions for TARGET_RTX_COSTS can be different from relative
   costs of hard registers to get the most efficient operations with
   pseudo registers.  */

struct processor_costs {
  /* Costs used by register allocator.  integer->integer register move
     cost is 2.  */
  struct
    {
      const int movzbl_load;	/* cost of loading using movzbl */
      const int int_load[3];	/* cost of loading integer registers
				   in QImode, HImode and SImode relative
				   to reg-reg move (2).  */
      const int int_store[3];	/* cost of storing integer register
				   in QImode, HImode and SImode */
      const int fp_move;	/* cost of reg,reg fld/fst */
      const int fp_load[3];	/* cost of loading FP register
				   in SFmode, DFmode and XFmode */
      const int fp_store[3];	/* cost of storing FP register
				   in SFmode, DFmode and XFmode */
      const int mmx_move;	/* cost of moving MMX register.  */
      const int mmx_load[2];	/* cost of loading MMX register
				   in SImode and DImode */
      const int mmx_store[2];	/* cost of storing MMX register
				   in SImode and DImode */
      const int xmm_move;	/* cost of moving XMM register.  */
      const int ymm_move;	/* cost of moving XMM register.  */
      const int zmm_move;	/* cost of moving XMM register.  */
      const int sse_load[5];	/* cost of loading SSE register
				   in 32bit, 64bit, 128bit, 256bit and 512bit */
      const int sse_store[5];	/* cost of storing SSE register
				   in SImode, DImode and TImode.  */
      const int sse_to_integer;	/* cost of moving SSE register to integer.  */
      const int integer_to_sse;	/* cost of moving integer register to SSE. */
      const int mask_to_integer; /* cost of moving mask register to integer.  */
      const int integer_to_mask; /* cost of moving integer register to mask.  */
      const int mask_load[3]; /* cost of loading mask registers
				 in QImode, HImode and SImode.  */
      const int mask_store[3]; /* cost of storing mask register
				  in QImode, HImode and SImode.  */
      const int mask_move; /* cost of moving mask register.  */
    } hard_register;

  const int add;		/* cost of an add instruction */
  const int lea;		/* cost of a lea instruction */
  const int shift_var;		/* variable shift costs */
  const int shift_const;	/* constant shift costs */
  const int mult_init[5];	/* cost of starting a multiply
				   in QImode, HImode, SImode, DImode, TImode*/
  const int mult_bit;		/* cost of multiply per each bit set */
  const int divide[5];		/* cost of a divide/mod
				   in QImode, HImode, SImode, DImode, TImode*/
  int movsx;			/* The cost of movsx operation.  */
  int movzx;			/* The cost of movzx operation.  */
  const int large_insn;		/* insns larger than this cost more */
  const int move_ratio;		/* The threshold of number of scalar
				   memory-to-memory move insns.  */
  const int clear_ratio;	/* The threshold of number of scalar
				   memory clearing insns.  */
  const int int_load[3];	/* cost of loading integer registers
				   in QImode, HImode and SImode relative
				   to reg-reg move (2).  */
  const int int_store[3];	/* cost of storing integer register
				   in QImode, HImode and SImode */
  const int sse_load[5];	/* cost of loading SSE register
				   in 32bit, 64bit, 128bit, 256bit and 512bit */
  const int sse_store[5];	/* cost of storing SSE register
				   in 32bit, 64bit, 128bit, 256bit and 512bit */
  const int sse_unaligned_load[5];/* cost of unaligned load.  */
  const int sse_unaligned_store[5];/* cost of unaligned store.  */
  const int xmm_move, ymm_move, /* cost of moving XMM and YMM register.  */
	    zmm_move;
  const int sse_to_integer;	/* cost of moving SSE register to integer.  */
  const int gather_static, gather_per_elt; /* Cost of gather load is computed
				   as static + per_item * nelts. */
  const int scatter_static, scatter_per_elt; /* Cost of gather store is
				   computed as static + per_item * nelts.  */
  const int l1_cache_size;	/* size of l1 cache, in kilobytes.  */
  const int l2_cache_size;	/* size of l2 cache, in kilobytes.  */
  const int prefetch_block;	/* bytes moved to cache for prefetch.  */
  const int simultaneous_prefetches; /* number of parallel prefetch
				   operations.  */
  const int branch_cost;	/* Default value for BRANCH_COST.  */
  const int fadd;		/* cost of FADD and FSUB instructions.  */
  const int fmul;		/* cost of FMUL instruction.  */
  const int fdiv;		/* cost of FDIV instruction.  */
  const int fabs;		/* cost of FABS instruction.  */
  const int fchs;		/* cost of FCHS instruction.  */
  const int fsqrt;		/* cost of FSQRT instruction.  */
				/* Specify what algorithm
				   to use for stringops on unknown size.  */
  const int sse_op;		/* cost of cheap SSE instruction.  */
  const int addss;		/* cost of ADDSS/SD SUBSS/SD instructions.  */
  const int mulss;		/* cost of MULSS instructions.  */
  const int mulsd;		/* cost of MULSD instructions.  */
  const int fmass;		/* cost of FMASS instructions.  */
  const int fmasd;		/* cost of FMASD instructions.  */
  const int divss;		/* cost of DIVSS instructions.  */
  const int divsd;		/* cost of DIVSD instructions.  */
  const int sqrtss;		/* cost of SQRTSS instructions.  */
  const int sqrtsd;		/* cost of SQRTSD instructions.  */
  const int reassoc_int, reassoc_fp, reassoc_vec_int, reassoc_vec_fp;
				/* Specify reassociation width for integer,
				   fp, vector integer and vector fp
				   operations.  Generally should correspond
				   to number of instructions executed in
				   parallel.  See also
				   ix86_reassociation_width.  */
  struct stringop_algs *memcpy, *memset;
  const int cond_taken_branch_cost;    /* Cost of taken branch for vectorizer
					  cost model.  */
  const int cond_not_taken_branch_cost;/* Cost of not taken branch for
					  vectorizer cost model.  */

  /* The "0:0:8" label alignment specified for some processors generates
     secondary 8-byte alignment only for those label/jump/loop targets
     which have primary alignment.  */
  const char *const align_loop;		/* Loop alignment.  */
  const char *const align_jump;		/* Jump alignment.  */
  const char *const align_label;	/* Label alignment.  */
  const char *const align_func;		/* Function alignment.  */
};

extern const struct processor_costs *ix86_cost;
extern const struct processor_costs ix86_size_cost;

#define ix86_cur_cost() \
  (optimize_insn_for_size_p () ? &ix86_size_cost: ix86_cost)

/* Macros used in the machine description to test the flags.  */

/* configure can arrange to change it.  */

#ifndef TARGET_CPU_DEFAULT
#define TARGET_CPU_DEFAULT PROCESSOR_GENERIC
#endif

#ifndef TARGET_FPMATH_DEFAULT
#define TARGET_FPMATH_DEFAULT \
  (TARGET_64BIT && TARGET_SSE ? FPMATH_SSE : FPMATH_387)
#endif

#ifndef TARGET_FPMATH_DEFAULT_P
#define TARGET_FPMATH_DEFAULT_P(x) \
  (TARGET_64BIT_P(x) && TARGET_SSE_P(x) ? FPMATH_SSE : FPMATH_387)
#endif

/* If the i387 is disabled or -miamcu is used , then do not return
   values in it. */
#define TARGET_FLOAT_RETURNS_IN_80387 \
  (TARGET_FLOAT_RETURNS && TARGET_80387 && !TARGET_IAMCU)
#define TARGET_FLOAT_RETURNS_IN_80387_P(x) \
  (TARGET_FLOAT_RETURNS_P(x) && TARGET_80387_P(x) && !TARGET_IAMCU_P(x))

/* 64bit Sledgehammer mode.  For libgcc2 we make sure this is a
   compile-time constant.  */
#ifdef IN_LIBGCC2
#undef TARGET_64BIT
#ifdef __x86_64__
#define TARGET_64BIT 1
#else
#define TARGET_64BIT 0
#endif
#else
#ifndef TARGET_BI_ARCH
#undef TARGET_64BIT
#undef TARGET_64BIT_P
#if TARGET_64BIT_DEFAULT
#define TARGET_64BIT 1
#define TARGET_64BIT_P(x) 1
#else
#define TARGET_64BIT 0
#define TARGET_64BIT_P(x) 0
#endif
#endif
#endif

#define HAS_LONG_COND_BRANCH 1
#define HAS_LONG_UNCOND_BRANCH 1

#define TARGET_CPU_P(CPU) (ix86_tune == PROCESSOR_ ## CPU)

/* Feature tests against the various tunings.  */
enum ix86_tune_indices {
#undef DEF_TUNE
#define DEF_TUNE(tune, name, selector) tune,
#include "x86-tune.def"
#undef DEF_TUNE
X86_TUNE_LAST
};

extern unsigned char ix86_tune_features[X86_TUNE_LAST];

#define TARGET_USE_LEAVE	ix86_tune_features[X86_TUNE_USE_LEAVE]
#define TARGET_PUSH_MEMORY	ix86_tune_features[X86_TUNE_PUSH_MEMORY]
#define TARGET_ZERO_EXTEND_WITH_AND \
	ix86_tune_features[X86_TUNE_ZERO_EXTEND_WITH_AND]
#define TARGET_UNROLL_STRLEN	ix86_tune_features[X86_TUNE_UNROLL_STRLEN]
#define TARGET_BRANCH_PREDICTION_HINTS \
	ix86_tune_features[X86_TUNE_BRANCH_PREDICTION_HINTS]
#define TARGET_DOUBLE_WITH_ADD	ix86_tune_features[X86_TUNE_DOUBLE_WITH_ADD]
#define TARGET_USE_SAHF		ix86_tune_features[X86_TUNE_USE_SAHF]
#define TARGET_MOVX		ix86_tune_features[X86_TUNE_MOVX]
#define TARGET_PARTIAL_REG_STALL ix86_tune_features[X86_TUNE_PARTIAL_REG_STALL]
#define TARGET_PARTIAL_FLAG_REG_STALL \
	ix86_tune_features[X86_TUNE_PARTIAL_FLAG_REG_STALL]
#define TARGET_LCP_STALL \
	ix86_tune_features[X86_TUNE_LCP_STALL]
#define TARGET_USE_HIMODE_FIOP	ix86_tune_features[X86_TUNE_USE_HIMODE_FIOP]
#define TARGET_USE_SIMODE_FIOP	ix86_tune_features[X86_TUNE_USE_SIMODE_FIOP]
#define TARGET_USE_MOV0		ix86_tune_features[X86_TUNE_USE_MOV0]
#define TARGET_USE_CLTD		ix86_tune_features[X86_TUNE_USE_CLTD]
#define TARGET_USE_XCHGB	ix86_tune_features[X86_TUNE_USE_XCHGB]
#define TARGET_SPLIT_LONG_MOVES	ix86_tune_features[X86_TUNE_SPLIT_LONG_MOVES]
#define TARGET_READ_MODIFY_WRITE ix86_tune_features[X86_TUNE_READ_MODIFY_WRITE]
#define TARGET_READ_MODIFY	ix86_tune_features[X86_TUNE_READ_MODIFY]
#define TARGET_PROMOTE_QImode	ix86_tune_features[X86_TUNE_PROMOTE_QIMODE]
#define TARGET_FAST_PREFIX	ix86_tune_features[X86_TUNE_FAST_PREFIX]
#define TARGET_SINGLE_STRINGOP	ix86_tune_features[X86_TUNE_SINGLE_STRINGOP]
#define TARGET_PREFER_KNOWN_REP_MOVSB_STOSB \
  ix86_tune_features[X86_TUNE_PREFER_KNOWN_REP_MOVSB_STOSB]
#define TARGET_MISALIGNED_MOVE_STRING_PRO_EPILOGUES \
	ix86_tune_features[X86_TUNE_MISALIGNED_MOVE_STRING_PRO_EPILOGUES]
#define TARGET_QIMODE_MATH	ix86_tune_features[X86_TUNE_QIMODE_MATH]
#define TARGET_HIMODE_MATH	ix86_tune_features[X86_TUNE_HIMODE_MATH]
#define TARGET_PROMOTE_QI_REGS	ix86_tune_features[X86_TUNE_PROMOTE_QI_REGS]
#define TARGET_PROMOTE_HI_REGS	ix86_tune_features[X86_TUNE_PROMOTE_HI_REGS]
#define TARGET_SINGLE_POP	ix86_tune_features[X86_TUNE_SINGLE_POP]
#define TARGET_DOUBLE_POP	ix86_tune_features[X86_TUNE_DOUBLE_POP]
#define TARGET_SINGLE_PUSH	ix86_tune_features[X86_TUNE_SINGLE_PUSH]
#define TARGET_DOUBLE_PUSH	ix86_tune_features[X86_TUNE_DOUBLE_PUSH]
#define TARGET_INTEGER_DFMODE_MOVES \
	ix86_tune_features[X86_TUNE_INTEGER_DFMODE_MOVES]
#define TARGET_PARTIAL_REG_DEPENDENCY \
	ix86_tune_features[X86_TUNE_PARTIAL_REG_DEPENDENCY]
#define TARGET_SSE_PARTIAL_REG_DEPENDENCY \
	ix86_tune_features[X86_TUNE_SSE_PARTIAL_REG_DEPENDENCY]
#define TARGET_SSE_PARTIAL_REG_FP_CONVERTS_DEPENDENCY \
	ix86_tune_features[X86_TUNE_SSE_PARTIAL_REG_FP_CONVERTS_DEPENDENCY]
#define TARGET_SSE_PARTIAL_REG_CONVERTS_DEPENDENCY \
	ix86_tune_features[X86_TUNE_SSE_PARTIAL_REG_CONVERTS_DEPENDENCY]
#define TARGET_SSE_UNALIGNED_LOAD_OPTIMAL \
	ix86_tune_features[X86_TUNE_SSE_UNALIGNED_LOAD_OPTIMAL]
#define TARGET_SSE_UNALIGNED_STORE_OPTIMAL \
	ix86_tune_features[X86_TUNE_SSE_UNALIGNED_STORE_OPTIMAL]
#define TARGET_SSE_PACKED_SINGLE_INSN_OPTIMAL \
	ix86_tune_features[X86_TUNE_SSE_PACKED_SINGLE_INSN_OPTIMAL]
#define TARGET_SSE_SPLIT_REGS	ix86_tune_features[X86_TUNE_SSE_SPLIT_REGS]
#define TARGET_SSE_TYPELESS_STORES \
	ix86_tune_features[X86_TUNE_SSE_TYPELESS_STORES]
#define TARGET_SSE_LOAD0_BY_PXOR ix86_tune_features[X86_TUNE_SSE_LOAD0_BY_PXOR]
#define TARGET_MEMORY_MISMATCH_STALL \
	ix86_tune_features[X86_TUNE_MEMORY_MISMATCH_STALL]
#define TARGET_PROLOGUE_USING_MOVE \
	ix86_tune_features[X86_TUNE_PROLOGUE_USING_MOVE]
#define TARGET_EPILOGUE_USING_MOVE \
	ix86_tune_features[X86_TUNE_EPILOGUE_USING_MOVE]
#define TARGET_SHIFT1		ix86_tune_features[X86_TUNE_SHIFT1]
#define TARGET_USE_FFREEP	ix86_tune_features[X86_TUNE_USE_FFREEP]
#define TARGET_INTER_UNIT_MOVES_TO_VEC \
	ix86_tune_features[X86_TUNE_INTER_UNIT_MOVES_TO_VEC]
#define TARGET_INTER_UNIT_MOVES_FROM_VEC \
	ix86_tune_features[X86_TUNE_INTER_UNIT_MOVES_FROM_VEC]
#define TARGET_INTER_UNIT_CONVERSIONS \
	ix86_tune_features[X86_TUNE_INTER_UNIT_CONVERSIONS]
#define TARGET_FOUR_JUMP_LIMIT	ix86_tune_features[X86_TUNE_FOUR_JUMP_LIMIT]
#define TARGET_SCHEDULE		ix86_tune_features[X86_TUNE_SCHEDULE]
#define TARGET_USE_BT		ix86_tune_features[X86_TUNE_USE_BT]
#define TARGET_USE_INCDEC	ix86_tune_features[X86_TUNE_USE_INCDEC]
#define TARGET_PAD_RETURNS	ix86_tune_features[X86_TUNE_PAD_RETURNS]
#define TARGET_PAD_SHORT_FUNCTION \
	ix86_tune_features[X86_TUNE_PAD_SHORT_FUNCTION]
#define TARGET_EXT_80387_CONSTANTS \
	ix86_tune_features[X86_TUNE_EXT_80387_CONSTANTS]
#define TARGET_AVOID_VECTOR_DECODE \
	ix86_tune_features[X86_TUNE_AVOID_VECTOR_DECODE]
#define TARGET_TUNE_PROMOTE_HIMODE_IMUL \
	ix86_tune_features[X86_TUNE_PROMOTE_HIMODE_IMUL]
#define TARGET_SLOW_IMUL_IMM32_MEM \
	ix86_tune_features[X86_TUNE_SLOW_IMUL_IMM32_MEM]
#define TARGET_SLOW_IMUL_IMM8	ix86_tune_features[X86_TUNE_SLOW_IMUL_IMM8]
#define	TARGET_MOVE_M1_VIA_OR	ix86_tune_features[X86_TUNE_MOVE_M1_VIA_OR]
#define TARGET_NOT_UNPAIRABLE	ix86_tune_features[X86_TUNE_NOT_UNPAIRABLE]
#define TARGET_NOT_VECTORMODE	ix86_tune_features[X86_TUNE_NOT_VECTORMODE]
#define TARGET_USE_VECTOR_FP_CONVERTS \
	ix86_tune_features[X86_TUNE_USE_VECTOR_FP_CONVERTS]
#define TARGET_USE_VECTOR_CONVERTS \
	ix86_tune_features[X86_TUNE_USE_VECTOR_CONVERTS]
#define TARGET_SLOW_PSHUFB \
	ix86_tune_features[X86_TUNE_SLOW_PSHUFB]
#define TARGET_AVOID_4BYTE_PREFIXES \
	ix86_tune_features[X86_TUNE_AVOID_4BYTE_PREFIXES]
#define TARGET_USE_GATHER_2PARTS \
	ix86_tune_features[X86_TUNE_USE_GATHER_2PARTS]
#define TARGET_USE_SCATTER_2PARTS \
	ix86_tune_features[X86_TUNE_USE_SCATTER_2PARTS]
#define TARGET_USE_GATHER_4PARTS \
	ix86_tune_features[X86_TUNE_USE_GATHER_4PARTS]
#define TARGET_USE_SCATTER_4PARTS \
	ix86_tune_features[X86_TUNE_USE_SCATTER_4PARTS]
#define TARGET_USE_GATHER \
	ix86_tune_features[X86_TUNE_USE_GATHER]
#define TARGET_USE_SCATTER \
	ix86_tune_features[X86_TUNE_USE_SCATTER]
#define TARGET_FUSE_CMP_AND_BRANCH_32 \
	ix86_tune_features[X86_TUNE_FUSE_CMP_AND_BRANCH_32]
#define TARGET_FUSE_CMP_AND_BRANCH_64 \
	ix86_tune_features[X86_TUNE_FUSE_CMP_AND_BRANCH_64]
#define TARGET_FUSE_CMP_AND_BRANCH \
	(TARGET_64BIT ? TARGET_FUSE_CMP_AND_BRANCH_64 \
	 : TARGET_FUSE_CMP_AND_BRANCH_32)
#define TARGET_FUSE_CMP_AND_BRANCH_SOFLAGS \
	ix86_tune_features[X86_TUNE_FUSE_CMP_AND_BRANCH_SOFLAGS]
#define TARGET_FUSE_ALU_AND_BRANCH \
	ix86_tune_features[X86_TUNE_FUSE_ALU_AND_BRANCH]
#define TARGET_OPT_AGU ix86_tune_features[X86_TUNE_OPT_AGU]
#define TARGET_AVOID_LEA_FOR_ADDR \
	ix86_tune_features[X86_TUNE_AVOID_LEA_FOR_ADDR]
#define TARGET_SOFTWARE_PREFETCHING_BENEFICIAL \
	ix86_tune_features[X86_TUNE_SOFTWARE_PREFETCHING_BENEFICIAL]
#define TARGET_AVX256_SPLIT_REGS \
	ix86_tune_features[X86_TUNE_AVX256_SPLIT_REGS]
#define TARGET_AVX512_SPLIT_REGS \
	ix86_tune_features[X86_TUNE_AVX512_SPLIT_REGS]
#define TARGET_GENERAL_REGS_SSE_SPILL \
	ix86_tune_features[X86_TUNE_GENERAL_REGS_SSE_SPILL]
#define TARGET_AVOID_MEM_OPND_FOR_CMOVE \
	ix86_tune_features[X86_TUNE_AVOID_MEM_OPND_FOR_CMOVE]
#define TARGET_SPLIT_MEM_OPND_FOR_FP_CONVERTS \
	ix86_tune_features[X86_TUNE_SPLIT_MEM_OPND_FOR_FP_CONVERTS]
#define TARGET_ADJUST_UNROLL \
    ix86_tune_features[X86_TUNE_ADJUST_UNROLL]
#define TARGET_AVOID_FALSE_DEP_FOR_BMI \
	ix86_tune_features[X86_TUNE_AVOID_FALSE_DEP_FOR_BMI]
#define TARGET_ONE_IF_CONV_INSN \
	ix86_tune_features[X86_TUNE_ONE_IF_CONV_INSN]
#define TARGET_AVOID_MFENCE ix86_tune_features[X86_TUNE_AVOID_MFENCE]
#define TARGET_EMIT_VZEROUPPER \
	ix86_tune_features[X86_TUNE_EMIT_VZEROUPPER]
#define TARGET_EXPAND_ABS \
	ix86_tune_features[X86_TUNE_EXPAND_ABS]
#define TARGET_V2DF_REDUCTION_PREFER_HADDPD \
	ix86_tune_features[X86_TUNE_V2DF_REDUCTION_PREFER_HADDPD]
#define TARGET_DEST_FALSE_DEP_FOR_GLC \
	ix86_tune_features[X86_TUNE_DEST_FALSE_DEP_FOR_GLC]

/* Feature tests against the various architecture variations.  */
enum ix86_arch_indices {
  X86_ARCH_CMOV,
  X86_ARCH_CMPXCHG,
  X86_ARCH_CMPXCHG8B,
  X86_ARCH_XADD,
  X86_ARCH_BSWAP,

  X86_ARCH_LAST
};

extern unsigned char ix86_arch_features[X86_ARCH_LAST];

#define TARGET_CMOV		ix86_arch_features[X86_ARCH_CMOV]
#define TARGET_CMPXCHG		ix86_arch_features[X86_ARCH_CMPXCHG]
#define TARGET_CMPXCHG8B	ix86_arch_features[X86_ARCH_CMPXCHG8B]
#define TARGET_XADD		ix86_arch_features[X86_ARCH_XADD]
#define TARGET_BSWAP		ix86_arch_features[X86_ARCH_BSWAP]

/* For sane SSE instruction set generation we need fcomi instruction.
   It is safe to enable all CMOVE instructions.  Also, RDRAND intrinsic
   expands to a sequence that includes conditional move. */
#define TARGET_CMOVE		(TARGET_CMOV || TARGET_SSE || TARGET_RDRND)

#define TARGET_FISTTP		(TARGET_SSE3 && TARGET_80387)

extern unsigned char ix86_prefetch_sse;
#define TARGET_PREFETCH_SSE	ix86_prefetch_sse

#define ASSEMBLER_DIALECT	(ix86_asm_dialect)

#define TARGET_SSE_MATH		((ix86_fpmath & FPMATH_SSE) != 0)
#define TARGET_MIX_SSE_I387 \
 ((ix86_fpmath & (FPMATH_SSE | FPMATH_387)) == (FPMATH_SSE | FPMATH_387))

#define TARGET_HARD_SF_REGS	(TARGET_80387 || TARGET_MMX || TARGET_SSE)
#define TARGET_HARD_DF_REGS	(TARGET_80387 || TARGET_SSE)
#define TARGET_HARD_XF_REGS	(TARGET_80387)

#define TARGET_GNU_TLS		(ix86_tls_dialect == TLS_DIALECT_GNU)
#define TARGET_GNU2_TLS		(ix86_tls_dialect == TLS_DIALECT_GNU2)
#define TARGET_ANY_GNU_TLS	(TARGET_GNU_TLS || TARGET_GNU2_TLS)
#define TARGET_SUN_TLS		0

#ifndef TARGET_64BIT_DEFAULT
#define TARGET_64BIT_DEFAULT 0
#endif
#ifndef TARGET_TLS_DIRECT_SEG_REFS_DEFAULT
#define TARGET_TLS_DIRECT_SEG_REFS_DEFAULT 0
#endif

#define TARGET_SSP_GLOBAL_GUARD (ix86_stack_protector_guard == SSP_GLOBAL)
#define TARGET_SSP_TLS_GUARD    (ix86_stack_protector_guard == SSP_TLS)

/* Fence to use after loop using storent.  */

extern GTY(()) tree x86_mfence;
#define FENCE_FOLLOWING_MOVNT x86_mfence

/* Once GDB has been enhanced to deal with functions without frame
   pointers, we can change this to allow for elimination of
   the frame pointer in leaf functions.  */
#define TARGET_DEFAULT 0

/* Extra bits to force.  */
#define TARGET_SUBTARGET_DEFAULT 0
#define TARGET_SUBTARGET_ISA_DEFAULT 0

/* Extra bits to force on w/ 32-bit mode.  */
#define TARGET_SUBTARGET32_DEFAULT 0
#define TARGET_SUBTARGET32_ISA_DEFAULT 0

/* Extra bits to force on w/ 64-bit mode.  */
#define TARGET_SUBTARGET64_DEFAULT 0
/* Enable MMX, SSE and SSE2 by default.  */
#define TARGET_SUBTARGET64_ISA_DEFAULT \
  (OPTION_MASK_ISA_MMX | OPTION_MASK_ISA_SSE | OPTION_MASK_ISA_SSE2)

/* Replace MACH-O, ifdefs by in-line tests, where possible. 
   (a) Macros defined in config/i386/darwin.h  */
#define TARGET_MACHO 0
#define TARGET_MACHO_SYMBOL_STUBS 0
#define MACHOPIC_ATT_STUB 0
/* (b) Macros defined in config/darwin.h  */
#define MACHO_DYNAMIC_NO_PIC_P 0
#define MACHOPIC_INDIRECT 0
#define MACHOPIC_PURE 0

/* For the RDOS  */
#define TARGET_RDOS 0

/* For the Windows 64-bit ABI.  */
#define TARGET_64BIT_MS_ABI (TARGET_64BIT && ix86_cfun_abi () == MS_ABI)

/* For the Windows 32-bit ABI.  */
#define TARGET_32BIT_MS_ABI (!TARGET_64BIT && ix86_cfun_abi () == MS_ABI)

/* This is re-defined by cygming.h.  */
#define TARGET_SEH 0

/* The default abi used by target.  */
#define DEFAULT_ABI SYSV_ABI

/* The default TLS segment register used by target.  */
#define DEFAULT_TLS_SEG_REG \
  (TARGET_64BIT ? ADDR_SPACE_SEG_FS : ADDR_SPACE_SEG_GS)

/* Subtargets may reset this to 1 in order to enable 96-bit long double
   with the rounding mode forced to 53 bits.  */
#define TARGET_96_ROUND_53_LONG_DOUBLE 0

#ifndef SUBTARGET_DRIVER_SELF_SPECS
# define SUBTARGET_DRIVER_SELF_SPECS ""
#endif

#define DRIVER_SELF_SPECS SUBTARGET_DRIVER_SELF_SPECS

/* -march=native handling only makes sense with compiler running on
   an x86 or x86_64 chip.  If changing this condition, also change
   the condition in driver-i386.cc.  */
#if defined(__i386__) || defined(__x86_64__)
/* In driver-i386.cc.  */
extern const char *host_detect_local_cpu (int argc, const char **argv);
#define EXTRA_SPEC_FUNCTIONS \
  { "local_cpu_detect", host_detect_local_cpu },
#define HAVE_LOCAL_CPU_DETECT
#endif

#if TARGET_64BIT_DEFAULT
#define OPT_ARCH64 "!m32"
#define OPT_ARCH32 "m32"
#else
#define OPT_ARCH64 "m64|mx32"
#define OPT_ARCH32 "m64|mx32:;"
#endif

/* Support for configure-time defaults of some command line options.
   The order here is important so that -march doesn't squash the
   tune or cpu values.  */
#define OPTION_DEFAULT_SPECS					   \
  {"tune", "%{!mtune=*:%{!mcpu=*:%{!march=*:-mtune=%(VALUE)}}}" }, \
  {"tune_32", "%{" OPT_ARCH32 ":%{!mtune=*:%{!mcpu=*:%{!march=*:-mtune=%(VALUE)}}}}" }, \
  {"tune_64", "%{" OPT_ARCH64 ":%{!mtune=*:%{!mcpu=*:%{!march=*:-mtune=%(VALUE)}}}}" }, \
  {"cpu", "%{!mtune=*:%{!mcpu=*:%{!march=*:-mtune=%(VALUE)}}}" },  \
  {"cpu_32", "%{" OPT_ARCH32 ":%{!mtune=*:%{!mcpu=*:%{!march=*:-mtune=%(VALUE)}}}}" }, \
  {"cpu_64", "%{" OPT_ARCH64 ":%{!mtune=*:%{!mcpu=*:%{!march=*:-mtune=%(VALUE)}}}}" }, \
  {"arch", "%{!march=*:-march=%(VALUE)}"},			   \
  {"arch_32", "%{" OPT_ARCH32 ":%{!march=*:-march=%(VALUE)}}"},	   \
  {"arch_64", "%{" OPT_ARCH64 ":%{!march=*:-march=%(VALUE)}}"},

/* Specs for the compiler proper */

#ifndef CC1_CPU_SPEC
#define CC1_CPU_SPEC_1 ""

#ifndef HAVE_LOCAL_CPU_DETECT
#define CC1_CPU_SPEC CC1_CPU_SPEC_1
#else
#define ARCH_ARG "%{" OPT_ARCH64 ":64;:32}"
#define CC1_CPU_SPEC CC1_CPU_SPEC_1 \
"%{march=native:%>march=native %:local_cpu_detect(arch " ARCH_ARG ") \
  %{!mtune=*:%>mtune=native %:local_cpu_detect(tune " ARCH_ARG ")}} \
%{mtune=native:%>mtune=native %:local_cpu_detect(tune " ARCH_ARG ")}"
#endif
#endif

/* Target CPU builtins.  */
#define TARGET_CPU_CPP_BUILTINS() ix86_target_macros ()

/* Target Pragmas.  */
#define REGISTER_TARGET_PRAGMAS() ix86_register_pragmas ()

#ifndef CC1_SPEC
#define CC1_SPEC "%(cc1_cpu) "
#endif

/* This macro defines names of additional specifications to put in the
   specs that can be used in various specifications like CC1_SPEC.  Its
   definition is an initializer with a subgrouping for each command option.

   Each subgrouping contains a string constant, that defines the
   specification name, and a string constant that used by the GCC driver
   program.

   Do not define this macro if it does not need to do anything.  */

#ifndef SUBTARGET_EXTRA_SPECS
#define SUBTARGET_EXTRA_SPECS
#endif

#define EXTRA_SPECS							\
  { "cc1_cpu",  CC1_CPU_SPEC },						\
  SUBTARGET_EXTRA_SPECS


/* Whether to allow x87 floating-point arithmetic on MODE (one of
   SFmode, DFmode and XFmode) in the current excess precision
   configuration.  */
#define X87_ENABLE_ARITH(MODE)				\
  (ix86_unsafe_math_optimizations			\
   || ix86_excess_precision == EXCESS_PRECISION_FAST	\
   || (MODE) == XFmode)

/* Likewise, whether to allow direct conversions from integer mode
   IMODE (HImode, SImode or DImode) to MODE.  */
#define X87_ENABLE_FLOAT(MODE, IMODE)			\
  (ix86_unsafe_math_optimizations			\
   || ix86_excess_precision == EXCESS_PRECISION_FAST	\
   || (MODE) == XFmode					\
   || ((MODE) == DFmode && (IMODE) == SImode)		\
   || (IMODE) == HImode)

/* target machine storage layout */

#define SHORT_TYPE_SIZE 16
#define INT_TYPE_SIZE 32
#define LONG_TYPE_SIZE (TARGET_X32 ? 32 : BITS_PER_WORD)
#define POINTER_SIZE (TARGET_X32 ? 32 : BITS_PER_WORD)
#define LONG_LONG_TYPE_SIZE 64
#define FLOAT_TYPE_SIZE 32
#define DOUBLE_TYPE_SIZE 64
#define LONG_DOUBLE_TYPE_SIZE \
  (TARGET_LONG_DOUBLE_64 ? 64 : (TARGET_LONG_DOUBLE_128 ? 128 : 80))

#define WIDEST_HARDWARE_FP_SIZE 80

#if defined (TARGET_BI_ARCH) || TARGET_64BIT_DEFAULT
#define MAX_BITS_PER_WORD 64
#else
#define MAX_BITS_PER_WORD 32
#endif

/* Define this if most significant byte of a word is the lowest numbered.  */
/* That is true on the 80386.  */

#define BITS_BIG_ENDIAN 0

/* Define this if most significant byte of a word is the lowest numbered.  */
/* That is not true on the 80386.  */
#define BYTES_BIG_ENDIAN 0

/* Define this if most significant word of a multiword number is the lowest
   numbered.  */
/* Not true for 80386 */
#define WORDS_BIG_ENDIAN 0

/* Width of a word, in units (bytes).  */
#define UNITS_PER_WORD		(TARGET_64BIT ? 8 : 4)

#ifndef IN_LIBGCC2
#define MIN_UNITS_PER_WORD	4
#endif

/* Allocation boundary (in *bits*) for storing arguments in argument list.  */
#define PARM_BOUNDARY BITS_PER_WORD

/* Boundary (in *bits*) on which stack pointer should be aligned.  */
#define STACK_BOUNDARY (TARGET_64BIT_MS_ABI ? 128 : BITS_PER_WORD)

/* Stack boundary of the main function guaranteed by OS.  */
#define MAIN_STACK_BOUNDARY (TARGET_64BIT ? 128 : 32)

/* Minimum stack boundary.  */
#define MIN_STACK_BOUNDARY BITS_PER_WORD

/* Boundary (in *bits*) on which the stack pointer prefers to be
   aligned; the compiler cannot rely on having this alignment.  */
#define PREFERRED_STACK_BOUNDARY ix86_preferred_stack_boundary

/* It should be MIN_STACK_BOUNDARY.  But we set it to 128 bits for
   both 32bit and 64bit, to support codes that need 128 bit stack
   alignment for SSE instructions, but can't realign the stack.  */
#define PREFERRED_STACK_BOUNDARY_DEFAULT \
  (TARGET_IAMCU ? MIN_STACK_BOUNDARY : 128)

/* 1 if -mstackrealign should be turned on by default.  It will
   generate an alternate prologue and epilogue that realigns the
   runtime stack if nessary.  This supports mixing codes that keep a
   4-byte aligned stack, as specified by i386 psABI, with codes that
   need a 16-byte aligned stack, as required by SSE instructions.  */
#define STACK_REALIGN_DEFAULT 0

/* Boundary (in *bits*) on which the incoming stack is aligned.  */
#define INCOMING_STACK_BOUNDARY ix86_incoming_stack_boundary

/* According to Windows x64 software convention, the maximum stack allocatable
   in the prologue is 4G - 8 bytes.  Furthermore, there is a limited set of
   instructions allowed to adjust the stack pointer in the epilog, forcing the
   use of frame pointer for frames larger than 2 GB.  This theorical limit
   is reduced by 256, an over-estimated upper bound for the stack use by the
   prologue.
   We define only one threshold for both the prolog and the epilog.  When the
   frame size is larger than this threshold, we allocate the area to save SSE
   regs, then save them, and then allocate the remaining.  There is no SEH
   unwind info for this later allocation.  */
#define SEH_MAX_FRAME_SIZE ((2U << 30) - 256)

/* Target OS keeps a vector-aligned (128-bit, 16-byte) stack.  This is
   mandatory for the 64-bit ABI, and may or may not be true for other
   operating systems.  */
#define TARGET_KEEPS_VECTOR_ALIGNED_STACK TARGET_64BIT

/* Minimum allocation boundary for the code of a function.  */
#define FUNCTION_BOUNDARY 8

/* C++ stores the virtual bit in the lowest bit of function pointers.  */
#define TARGET_PTRMEMFUNC_VBIT_LOCATION ptrmemfunc_vbit_in_pfn

/* Minimum size in bits of the largest boundary to which any
   and all fundamental data types supported by the hardware
   might need to be aligned. No data type wants to be aligned
   rounder than this.

   Pentium+ prefers DFmode values to be aligned to 64 bit boundary
   and Pentium Pro XFmode values at 128 bit boundaries.

   When increasing the maximum, also update
   TARGET_ABSOLUTE_BIGGEST_ALIGNMENT.  */

#define BIGGEST_ALIGNMENT \
  (TARGET_IAMCU ? 32 : (TARGET_AVX512F ? 512 : (TARGET_AVX ? 256 : 128)))

/* Maximum stack alignment.  */
#define MAX_STACK_ALIGNMENT MAX_OFILE_ALIGNMENT

/* Alignment value for attribute ((aligned)).  It is a constant since
   it is the part of the ABI.  We shouldn't change it with -mavx.  */
#define ATTRIBUTE_ALIGNED_VALUE (TARGET_IAMCU ? 32 : 128)

/* Decide whether a variable of mode MODE should be 128 bit aligned.  */
#define ALIGN_MODE_128(MODE) \
 ((MODE) == XFmode || SSE_REG_MODE_P (MODE))

/* The published ABIs say that doubles should be aligned on word
   boundaries, so lower the alignment for structure fields unless
   -malign-double is set.  */

/* ??? Blah -- this macro is used directly by libobjc.  Since it
   supports no vector modes, cut out the complexity and fall back
   on BIGGEST_FIELD_ALIGNMENT.  */
#ifdef IN_TARGET_LIBS
#ifdef __x86_64__
#define BIGGEST_FIELD_ALIGNMENT 128
#else
#define BIGGEST_FIELD_ALIGNMENT 32
#endif
#else
#define ADJUST_FIELD_ALIGN(FIELD, TYPE, COMPUTED) \
  x86_field_alignment ((TYPE), (COMPUTED))
#endif

/* If defined, a C expression to compute the alignment for a static
   variable.  TYPE is the data type, and ALIGN is the alignment that
   the object would ordinarily have.  The value of this macro is used
   instead of that alignment to align the object.

   If this macro is not defined, then ALIGN is used.

   One use of this macro is to increase alignment of medium-size
   data to make it all fit in fewer cache lines.  Another is to
   cause character arrays to be word-aligned so that `strcpy' calls
   that copy constants to character arrays can be done inline.  */

#define DATA_ALIGNMENT(TYPE, ALIGN) \
  ix86_data_alignment ((TYPE), (ALIGN), true)

/* Similar to DATA_ALIGNMENT, but for the cases where the ABI mandates
   some alignment increase, instead of optimization only purposes.  E.g.
   AMD x86-64 psABI says that variables with array type larger than 15 bytes
   must be aligned to 16 byte boundaries.

   If this macro is not defined, then ALIGN is used.  */

#define DATA_ABI_ALIGNMENT(TYPE, ALIGN) \
  ix86_data_alignment ((TYPE), (ALIGN), false)

/* If defined, a C expression to compute the alignment for a local
   variable.  TYPE is the data type, and ALIGN is the alignment that
   the object would ordinarily have.  The value of this macro is used
   instead of that alignment to align the object.

   If this macro is not defined, then ALIGN is used.

   One use of this macro is to increase alignment of medium-size
   data to make it all fit in fewer cache lines.  */

#define LOCAL_ALIGNMENT(TYPE, ALIGN) \
  ix86_local_alignment ((TYPE), VOIDmode, (ALIGN))

/* If defined, a C expression to compute the alignment for stack slot.
   TYPE is the data type, MODE is the widest mode available, and ALIGN
   is the alignment that the slot would ordinarily have.  The value of
   this macro is used instead of that alignment to align the slot.

   If this macro is not defined, then ALIGN is used when TYPE is NULL,
   Otherwise, LOCAL_ALIGNMENT will be used.

   One use of this macro is to set alignment of stack slot to the
   maximum alignment of all possible modes which the slot may have.  */

#define STACK_SLOT_ALIGNMENT(TYPE, MODE, ALIGN) \
  ix86_local_alignment ((TYPE), (MODE), (ALIGN))

/* If defined, a C expression to compute the alignment for a local
   variable DECL.

   If this macro is not defined, then
   LOCAL_ALIGNMENT (TREE_TYPE (DECL), DECL_ALIGN (DECL)) will be used.

   One use of this macro is to increase alignment of medium-size
   data to make it all fit in fewer cache lines.  */

#define LOCAL_DECL_ALIGNMENT(DECL) \
  ix86_local_alignment ((DECL), VOIDmode, DECL_ALIGN (DECL))

/* If defined, a C expression to compute the minimum required alignment
   for dynamic stack realignment purposes for EXP (a TYPE or DECL),
   MODE, assuming normal alignment ALIGN.

   If this macro is not defined, then (ALIGN) will be used.  */

#define MINIMUM_ALIGNMENT(EXP, MODE, ALIGN) \
  ix86_minimum_alignment ((EXP), (MODE), (ALIGN))


/* Set this nonzero if move instructions will actually fail to work
   when given unaligned data.  */
#define STRICT_ALIGNMENT 0

/* If bit field type is int, don't let it cross an int,
   and give entire struct the alignment of an int.  */
/* Required on the 386 since it doesn't have bit-field insns.  */
#define PCC_BITFIELD_TYPE_MATTERS 1

/* Standard register usage.  */

/* This processor has special stack-like registers.  See reg-stack.cc
   for details.  */

#define STACK_REGS

#define IS_STACK_MODE(MODE)				\
  (X87_FLOAT_MODE_P (MODE)				\
   && (!(SSE_FLOAT_MODE_P (MODE) && TARGET_SSE_MATH)	\
       || TARGET_MIX_SSE_I387))

/* Number of actual hardware registers.
   The hardware registers are assigned numbers for the compiler
   from 0 to just below FIRST_PSEUDO_REGISTER.
   All registers that the compiler knows about must be given numbers,
   even those that are not normally considered general registers.

   In the 80386 we give the 8 general purpose registers the numbers 0-7.
   We number the floating point registers 8-15.
   Note that registers 0-7 can be accessed as a  short or int,
   while only 0-3 may be used with byte `mov' instructions.

   Reg 16 does not correspond to any hardware register, but instead
   appears in the RTL as an argument pointer prior to reload, and is
   eliminated during reloading in favor of either the stack or frame
   pointer.  */

#define FIRST_PSEUDO_REGISTER FIRST_PSEUDO_REG

/* Number of hardware registers that go into the DWARF-2 unwind info.
   If not defined, equals FIRST_PSEUDO_REGISTER.  */

#define DWARF_FRAME_REGISTERS 17

/* 1 for registers that have pervasive standard uses
   and are not available for the register allocator.
   On the 80386, the stack pointer is such, as is the arg pointer.

   REX registers are disabled for 32bit targets in
   TARGET_CONDITIONAL_REGISTER_USAGE.  */

#define FIXED_REGISTERS						\
/*ax,dx,cx,bx,si,di,bp,sp,st,st1,st2,st3,st4,st5,st6,st7*/	\
{  0, 0, 0, 0, 0, 0, 0, 1, 0,  0,  0,  0,  0,  0,  0,  0,	\
/*arg,flags,fpsr,frame*/					\
    1,    1,   1,    1,						\
/*xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7*/			\
     0,   0,   0,   0,   0,   0,   0,   0,			\
/* mm0, mm1, mm2, mm3, mm4, mm5, mm6, mm7*/			\
     0,   0,   0,   0,   0,   0,   0,   0,			\
/*  r8,  r9, r10, r11, r12, r13, r14, r15*/			\
     0,   0,   0,   0,   0,   0,   0,   0,			\
/*xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15*/		\
     0,   0,    0,    0,    0,    0,    0,    0,		\
/*xmm16,xmm17,xmm18,xmm19,xmm20,xmm21,xmm22,xmm23*/		\
     0,   0,    0,    0,    0,    0,    0,    0,		\
/*xmm24,xmm25,xmm26,xmm27,xmm28,xmm29,xmm30,xmm31*/		\
     0,   0,    0,    0,    0,    0,    0,    0,		\
/*  k0,  k1, k2, k3, k4, k5, k6, k7*/				\
     0,  0,   0,  0,  0,  0,  0,  0 }

/* 1 for registers not available across function calls.
   These must include the FIXED_REGISTERS and also any
   registers that can be used without being saved.
   The latter must include the registers where values are returned
   and the register where structure-value addresses are passed.
   Aside from that, you can include as many other registers as you like.

   Value is set to 1 if the register is call used unconditionally.
   Bit one is set if the register is call used on TARGET_32BIT ABI.
   Bit two is set if the register is call used on TARGET_64BIT ABI.
   Bit three is set if the register is call used on TARGET_64BIT_MS_ABI.

   Proper values are computed in TARGET_CONDITIONAL_REGISTER_USAGE.  */

#define CALL_USED_REGISTERS_MASK(IS_64BIT_MS_ABI) \
  ((IS_64BIT_MS_ABI) ? (1 << 3) : TARGET_64BIT ? (1 << 2) : (1 << 1))

#define CALL_USED_REGISTERS					\
/*ax,dx,cx,bx,si,di,bp,sp,st,st1,st2,st3,st4,st5,st6,st7*/	\
{  1, 1, 1, 0, 4, 4, 0, 1, 1,  1,  1,  1,  1,  1,  1,  1,	\
/*arg,flags,fpsr,frame*/					\
    1,   1,    1,    1,						\
/*xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7*/			\
     1,   1,   1,   1,   1,   1,   6,   6,			\
/* mm0, mm1, mm2, mm3, mm4, mm5, mm6, mm7*/			\
     1,   1,   1,   1,   1,   1,   1,   1,			\
/*  r8,  r9, r10, r11, r12, r13, r14, r15*/			\
     1,   1,   1,   1,   2,   2,   2,   2,			\
/*xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15*/		\
     6,   6,    6,    6,    6,    6,    6,    6,		\
/*xmm16,xmm17,xmm18,xmm19,xmm20,xmm21,xmm22,xmm23*/		\
     1,    1,     1,    1,    1,    1,    1,    1,		\
/*xmm24,xmm25,xmm26,xmm27,xmm28,xmm29,xmm30,xmm31*/		\
     1,    1,     1,    1,    1,    1,    1,    1,		\
 /* k0,  k1,  k2,  k3,  k4,  k5,  k6,  k7*/			\
     1,   1,   1,   1,   1,   1,   1,   1 }

/* Order in which to allocate registers.  Each register must be
   listed once, even those in FIXED_REGISTERS.  List frame pointer
   late and fixed registers last.  Note that, in general, we prefer
   registers listed in CALL_USED_REGISTERS, keeping the others
   available for storage of persistent values.

   The ADJUST_REG_ALLOC_ORDER actually overwrite the order,
   so this is just empty initializer for array.  */

#define REG_ALLOC_ORDER							\
{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,			\
  16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,	\
  32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,	\
  48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,	\
  64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75 }

/* ADJUST_REG_ALLOC_ORDER is a macro which permits reg_alloc_order
   to be rearranged based on a particular function.  When using sse math,
   we want to allocate SSE before x87 registers and vice versa.  */

#define ADJUST_REG_ALLOC_ORDER x86_order_regs_for_local_alloc ()


#define OVERRIDE_ABI_FORMAT(FNDECL) ix86_call_abi_override (FNDECL)

#define HARD_REGNO_NREGS_HAS_PADDING(REGNO, MODE)			\
  (TARGET_128BIT_LONG_DOUBLE && !TARGET_64BIT				\
   && GENERAL_REGNO_P (REGNO)						\
   && ((MODE) == XFmode || (MODE) == XCmode))

#define HARD_REGNO_NREGS_WITH_PADDING(REGNO, MODE) ((MODE) == XFmode ? 4 : 8)

#define REGMODE_NATURAL_SIZE(MODE) ix86_regmode_natural_size (MODE)

#define VALID_AVX256_REG_MODE(MODE)					\
  ((MODE) == V32QImode || (MODE) == V16HImode || (MODE) == V8SImode	\
   || (MODE) == V4DImode || (MODE) == V2TImode || (MODE) == V8SFmode	\
   || (MODE) == V4DFmode || (MODE) == V16HFmode)

#define VALID_AVX256_REG_OR_OI_MODE(MODE)		\
  (VALID_AVX256_REG_MODE (MODE) || (MODE) == OImode)

#define VALID_AVX512F_SCALAR_MODE(MODE)					\
  ((MODE) == DImode || (MODE) == DFmode || (MODE) == SImode		\
   || (MODE) == SFmode)

#define VALID_AVX512FP16_SCALAR_MODE(MODE)	\
  ((MODE) == HImode || (MODE) == HFmode)

#define VALID_AVX512F_REG_MODE(MODE)					\
  ((MODE) == V8DImode || (MODE) == V8DFmode || (MODE) == V64QImode	\
   || (MODE) == V16SImode || (MODE) == V16SFmode || (MODE) == V32HImode \
   || (MODE) == V4TImode || (MODE) == V32HFmode)

#define VALID_AVX512F_REG_OR_XI_MODE(MODE)				\
  (VALID_AVX512F_REG_MODE (MODE) || (MODE) == XImode)

#define VALID_AVX512VL_128_REG_MODE(MODE)				\
  ((MODE) == V2DImode || (MODE) == V2DFmode || (MODE) == V16QImode	\
   || (MODE) == V4SImode || (MODE) == V4SFmode || (MODE) == V8HImode	\
   || (MODE) == TFmode || (MODE) == V1TImode || (MODE) == V8HFmode	\
   || (MODE) == TImode)

#define VALID_AVX512FP16_REG_MODE(MODE)					\
  ((MODE) == V8HFmode || (MODE) == V16HFmode || (MODE) == V32HFmode	\
   || (MODE) == V2HFmode)

#define VALID_SSE2_REG_MODE(MODE)					\
  ((MODE) == V16QImode || (MODE) == V8HImode || (MODE) == V2DFmode	\
   || (MODE) == V8HFmode || (MODE) == V4HFmode || (MODE) == V2HFmode	\
   || (MODE) == V4QImode || (MODE) == V2HImode || (MODE) == V1SImode	\
   || (MODE) == V2DImode || (MODE) == V2QImode || (MODE) == DFmode	\
   || (MODE) == HFmode)

#define VALID_SSE_REG_MODE(MODE)					\
  ((MODE) == V1TImode || (MODE) == TImode				\
   || (MODE) == V4SFmode || (MODE) == V4SImode				\
   || (MODE) == SFmode || (MODE) == TFmode || (MODE) == TDmode)

#define VALID_MMX_REG_MODE_3DNOW(MODE) \
  ((MODE) == V2SFmode || (MODE) == SFmode)

/* To match ia32 psABI, V4HFmode should be added here.  */
#define VALID_MMX_REG_MODE(MODE)					\
  ((MODE) == V1DImode || (MODE) == DImode				\
   || (MODE) == V2SImode || (MODE) == SImode				\
   || (MODE) == V4HImode || (MODE) == V8QImode				\
   || (MODE) == V4HFmode)

#define VALID_MASK_REG_MODE(MODE) ((MODE) == HImode || (MODE) == QImode)

#define VALID_MASK_AVX512BW_MODE(MODE) ((MODE) == SImode || (MODE) == DImode)

#define VALID_FP_MODE_P(MODE)						\
  ((MODE) == SFmode || (MODE) == DFmode || (MODE) == XFmode		\
   || (MODE) == SCmode || (MODE) == DCmode || (MODE) == XCmode)

#define VALID_INT_MODE_P(MODE)						\
  ((MODE) == QImode || (MODE) == HImode					\
   || (MODE) == SImode || (MODE) == DImode				\
   || (MODE) == CQImode || (MODE) == CHImode				\
   || (MODE) == CSImode || (MODE) == CDImode				\
   || (MODE) == SDmode || (MODE) == DDmode				\
   || (MODE) == HFmode || (MODE) == HCmode				\
   || (MODE) == V2HImode || (MODE) == V2HFmode				\
   || (MODE) == V1SImode || (MODE) == V4QImode || (MODE) == V2QImode	\
   || (TARGET_64BIT							\
       && ((MODE) == TImode || (MODE) == CTImode			\
	   || (MODE) == TFmode || (MODE) == TCmode			\
	   || (MODE) == V8QImode || (MODE) == V4HImode			\
	   || (MODE) == V2SImode || (MODE) == TDmode)))

/* Return true for modes passed in SSE registers.  */
#define SSE_REG_MODE_P(MODE)						\
  ((MODE) == V1TImode || (MODE) == TImode || (MODE) == V16QImode	\
   || (MODE) == TFmode || (MODE) == V8HImode || (MODE) == V2DFmode	\
   || (MODE) == V2DImode || (MODE) == V4SFmode || (MODE) == V4SImode	\
   || (MODE) == V32QImode || (MODE) == V16HImode || (MODE) == V8SImode	\
   || (MODE) == V4DImode || (MODE) == V8SFmode || (MODE) == V4DFmode	\
   || (MODE) == V2TImode || (MODE) == V8DImode || (MODE) == V64QImode	\
   || (MODE) == V16SImode || (MODE) == V32HImode || (MODE) == V8DFmode	\
   || (MODE) == V16SFmode || (MODE) == V32HFmode || (MODE) == V16HFmode \
   || (MODE) == V8HFmode)

#define X87_FLOAT_MODE_P(MODE)	\
  (TARGET_80387 && ((MODE) == SFmode || (MODE) == DFmode || (MODE) == XFmode))

#define SSE_FLOAT_MODE_P(MODE) \
  ((TARGET_SSE && (MODE) == SFmode) || (TARGET_SSE2 && (MODE) == DFmode))

#define SSE_FLOAT_MODE_SSEMATH_OR_HF_P(MODE)				\
  ((SSE_FLOAT_MODE_P (MODE) && TARGET_SSE_MATH)				\
   || (TARGET_AVX512FP16 && (MODE) == HFmode))

#define FMA4_VEC_FLOAT_MODE_P(MODE) \
  (TARGET_FMA4 && ((MODE) == V4SFmode || (MODE) == V2DFmode \
		  || (MODE) == V8SFmode || (MODE) == V4DFmode))

#define VALID_BCST_MODE_P(MODE)			\
  ((MODE) == SFmode || (MODE) == DFmode		\
   || (MODE) == SImode || (MODE) == DImode	\
   || (MODE) == HFmode)

/* It is possible to write patterns to move flags; but until someone
   does it,  */
#define AVOID_CCMODE_COPIES

/* Specify the modes required to caller save a given hard regno.
   We do this on i386 to prevent flags from being saved at all.

   Kill any attempts to combine saving of modes.  */

#define HARD_REGNO_CALLER_SAVE_MODE(REGNO, NREGS, MODE)			\
  (CC_REGNO_P (REGNO) ? VOIDmode					\
   : (MODE) == VOIDmode && (NREGS) != 1 ? VOIDmode			\
   : (MODE) == VOIDmode ? choose_hard_reg_mode ((REGNO), (NREGS), NULL)	\
   : (MODE) == HImode && !((GENERAL_REGNO_P (REGNO)			\
			    && TARGET_PARTIAL_REG_STALL)		\
			   || MASK_REGNO_P (REGNO)) ? SImode		\
   : (MODE) == QImode && !(ANY_QI_REGNO_P (REGNO)			\
			   || MASK_REGNO_P (REGNO)) ? SImode		\
   : (MODE))

/* Specify the registers used for certain standard purposes.
   The values of these macros are register numbers.  */

/* on the 386 the pc register is %eip, and is not usable as a general
   register.  The ordinary mov instructions won't work */
/* #define PC_REGNUM  */

/* Base register for access to arguments of the function.  */
#define ARG_POINTER_REGNUM ARGP_REG

/* Register to use for pushing function arguments.  */
#define STACK_POINTER_REGNUM SP_REG

/* Base register for access to local variables of the function.  */
#define FRAME_POINTER_REGNUM FRAME_REG
#define HARD_FRAME_POINTER_REGNUM BP_REG

#define FIRST_INT_REG AX_REG
#define LAST_INT_REG  SP_REG

#define FIRST_QI_REG AX_REG
#define LAST_QI_REG  BX_REG

/* First & last stack-like regs */
#define FIRST_STACK_REG ST0_REG
#define LAST_STACK_REG  ST7_REG

#define FIRST_SSE_REG XMM0_REG
#define LAST_SSE_REG  XMM7_REG

#define FIRST_MMX_REG  MM0_REG
#define LAST_MMX_REG   MM7_REG

#define FIRST_REX_INT_REG  R8_REG
#define LAST_REX_INT_REG   R15_REG

#define FIRST_REX_SSE_REG  XMM8_REG
#define LAST_REX_SSE_REG   XMM15_REG

#define FIRST_EXT_REX_SSE_REG  XMM16_REG
#define LAST_EXT_REX_SSE_REG   XMM31_REG

#define FIRST_MASK_REG  MASK0_REG
#define LAST_MASK_REG   MASK7_REG

/* Override this in other tm.h files to cope with various OS lossage
   requiring a frame pointer.  */
#ifndef SUBTARGET_FRAME_POINTER_REQUIRED
#define SUBTARGET_FRAME_POINTER_REQUIRED 0
#endif

/* Define the shadow offset for asan. Other OS's can override in the
   respective tm.h files.  */
#ifndef SUBTARGET_SHADOW_OFFSET
#define SUBTARGET_SHADOW_OFFSET	    \
  (TARGET_LP64 ? HOST_WIDE_INT_C (0x7fff8000) : HOST_WIDE_INT_1 << 29)
#endif

/* Make sure we can access arbitrary call frames.  */
#define SETUP_FRAME_ADDRESSES()  ix86_setup_frame_addresses ()

/* Register to hold the addressing base for position independent
   code access to data items.  We don't use PIC pointer for 64bit
   mode.  Define the regnum to dummy value to prevent gcc from
   pessimizing code dealing with EBX.

   To avoid clobbering a call-saved register unnecessarily, we renumber
   the pic register when possible.  The change is visible after the
   prologue has been emitted.  */

#define REAL_PIC_OFFSET_TABLE_REGNUM  (TARGET_64BIT ? R15_REG : BX_REG)

#define PIC_OFFSET_TABLE_REGNUM						\
  (ix86_use_pseudo_pic_reg ()						\
   ? (pic_offset_table_rtx						\
      ? INVALID_REGNUM							\
      : REAL_PIC_OFFSET_TABLE_REGNUM)					\
   : INVALID_REGNUM)

#define GOT_SYMBOL_NAME "_GLOBAL_OFFSET_TABLE_"

/* This is overridden by <cygwin.h>.  */
#define MS_AGGREGATE_RETURN 0

#define KEEP_AGGREGATE_RETURN_POINTER 0

/* Define the classes of registers for register constraints in the
   machine description.  Also define ranges of constants.

   One of the classes must always be named ALL_REGS and include all hard regs.
   If there is more than one class, another class must be named NO_REGS
   and contain no registers.

   The name GENERAL_REGS must be the name of a class (or an alias for
   another name such as ALL_REGS).  This is the class of registers
   that is allowed by "g" or "r" in a register constraint.
   Also, registers outside this class are allocated only when
   instructions express preferences for them.

   The classes must be numbered in nondecreasing order; that is,
   a larger-numbered class must never be contained completely
   in a smaller-numbered class.  This is why CLOBBERED_REGS class
   is listed early, even though in 64-bit mode it contains more
   registers than just %eax, %ecx, %edx.

   For any two classes, it is very desirable that there be another
   class that represents their union.

   The flags and fpsr registers are in no class.  */

enum reg_class
{
  NO_REGS,
  AREG, DREG, CREG, BREG, SIREG, DIREG,
  AD_REGS,			/* %eax/%edx for DImode */
  CLOBBERED_REGS,		/* call-clobbered integer registers */
  Q_REGS,			/* %eax %ebx %ecx %edx */
  NON_Q_REGS,			/* %esi %edi %ebp %esp */
  TLS_GOTBASE_REGS,		/* %ebx %ecx %edx %esi %edi %ebp */
  INDEX_REGS,			/* %eax %ebx %ecx %edx %esi %edi %ebp */
  LEGACY_REGS,			/* %eax %ebx %ecx %edx %esi %edi %ebp %esp */
  GENERAL_REGS,			/* %eax %ebx %ecx %edx %esi %edi %ebp %esp
				   %r8 %r9 %r10 %r11 %r12 %r13 %r14 %r15 */
  FP_TOP_REG, FP_SECOND_REG,	/* %st(0) %st(1) */
  FLOAT_REGS,
  SSE_FIRST_REG,
  NO_REX_SSE_REGS,
  SSE_REGS,
  ALL_SSE_REGS,
  MMX_REGS,
  FLOAT_SSE_REGS,
  FLOAT_INT_REGS,
  INT_SSE_REGS,
  FLOAT_INT_SSE_REGS,
  MASK_REGS,
  ALL_MASK_REGS,
  INT_MASK_REGS,
  ALL_REGS,
  LIM_REG_CLASSES
};

#define N_REG_CLASSES ((int) LIM_REG_CLASSES)

#define INTEGER_CLASS_P(CLASS) \
  reg_class_subset_p ((CLASS), GENERAL_REGS)
#define FLOAT_CLASS_P(CLASS) \
  reg_class_subset_p ((CLASS), FLOAT_REGS)
#define SSE_CLASS_P(CLASS) \
  reg_class_subset_p ((CLASS), ALL_SSE_REGS)
#define INT_SSE_CLASS_P(CLASS) \
  reg_class_subset_p ((CLASS), INT_SSE_REGS)
#define MMX_CLASS_P(CLASS) \
  ((CLASS) == MMX_REGS)
#define MASK_CLASS_P(CLASS) \
  reg_class_subset_p ((CLASS), ALL_MASK_REGS)
#define MAYBE_INTEGER_CLASS_P(CLASS) \
  reg_classes_intersect_p ((CLASS), GENERAL_REGS)
#define MAYBE_FLOAT_CLASS_P(CLASS) \
  reg_classes_intersect_p ((CLASS), FLOAT_REGS)
#define MAYBE_SSE_CLASS_P(CLASS) \
  reg_classes_intersect_p ((CLASS), ALL_SSE_REGS)
#define MAYBE_MMX_CLASS_P(CLASS) \
  reg_classes_intersect_p ((CLASS), MMX_REGS)
#define MAYBE_MASK_CLASS_P(CLASS) \
  reg_classes_intersect_p ((CLASS), ALL_MASK_REGS)

#define Q_CLASS_P(CLASS) \
  reg_class_subset_p ((CLASS), Q_REGS)

#define MAYBE_NON_Q_CLASS_P(CLASS) \
  reg_classes_intersect_p ((CLASS), NON_Q_REGS)

/* Give names of register classes as strings for dump file.  */

#define REG_CLASS_NAMES \
{  "NO_REGS",				\
   "AREG", "DREG", "CREG", "BREG",	\
   "SIREG", "DIREG",			\
   "AD_REGS",				\
   "CLOBBERED_REGS",			\
   "Q_REGS", "NON_Q_REGS",		\
   "TLS_GOTBASE_REGS",			\
   "INDEX_REGS",			\
   "LEGACY_REGS",			\
   "GENERAL_REGS",			\
   "FP_TOP_REG", "FP_SECOND_REG",	\
   "FLOAT_REGS",			\
   "SSE_FIRST_REG",			\
   "NO_REX_SSE_REGS",			\
   "SSE_REGS",				\
   "ALL_SSE_REGS",			\
   "MMX_REGS",				\
   "FLOAT_SSE_REGS",			\
   "FLOAT_INT_REGS",			\
   "INT_SSE_REGS",			\
   "FLOAT_INT_SSE_REGS",		\
   "MASK_REGS",				\
   "ALL_MASK_REGS",			\
   "INT_MASK_REGS",			\
   "ALL_REGS" }

/* Define which registers fit in which classes.  This is an initializer
   for a vector of HARD_REG_SET of length N_REG_CLASSES.

   Note that CLOBBERED_REGS are calculated by
   TARGET_CONDITIONAL_REGISTER_USAGE.  */

#define REG_CLASS_CONTENTS						\
{      { 0x0,        0x0,   0x0 },	/* NO_REGS */			\
      { 0x01,        0x0,   0x0 },	/* AREG */			\
      { 0x02,        0x0,   0x0 },	/* DREG */			\
      { 0x04,        0x0,   0x0 },	/* CREG */			\
      { 0x08,        0x0,   0x0 },	/* BREG */			\
      { 0x10,        0x0,   0x0 },	/* SIREG */			\
      { 0x20,        0x0,   0x0 },	/* DIREG */			\
      { 0x03,        0x0,   0x0 },	/* AD_REGS */			\
      { 0x07,        0x0,   0x0 },	/* CLOBBERED_REGS */		\
      { 0x0f,        0x0,   0x0 },	/* Q_REGS */			\
   { 0x900f0,        0x0,   0x0 },	/* NON_Q_REGS */		\
      { 0x7e,      0xff0,   0x0 },	/* TLS_GOTBASE_REGS */		\
      { 0x7f,      0xff0,   0x0 },	/* INDEX_REGS */		\
   { 0x900ff,        0x0,   0x0 },	/* LEGACY_REGS */		\
   { 0x900ff,      0xff0,   0x0 },	/* GENERAL_REGS */		\
     { 0x100,        0x0,   0x0 },	/* FP_TOP_REG */		\
     { 0x200,        0x0,   0x0 },	/* FP_SECOND_REG */		\
    { 0xff00,        0x0,   0x0 },	/* FLOAT_REGS */		\
  { 0x100000,        0x0,   0x0 },	/* SSE_FIRST_REG */		\
 { 0xff00000,        0x0,   0x0 },	/* NO_REX_SSE_REGS */		\
 { 0xff00000,    0xff000,   0x0 },	/* SSE_REGS */			\
 { 0xff00000, 0xfffff000,   0xf },	/* ALL_SSE_REGS */		\
{ 0xf0000000,        0xf,   0x0 },	/* MMX_REGS */			\
 { 0xff0ff00, 0xfffff000,   0xf },	/* FLOAT_SSE_REGS */		\
 {   0x9ffff,      0xff0,   0x0 },	/* FLOAT_INT_REGS */		\
 { 0xff900ff, 0xfffffff0,   0xf },	/* INT_SSE_REGS */		\
 { 0xff9ffff, 0xfffffff0,   0xf },	/* FLOAT_INT_SSE_REGS */	\
       { 0x0,        0x0, 0xfe0 },	/* MASK_REGS */			\
       { 0x0,        0x0, 0xff0 },	/* ALL_MASK_REGS */		\
   { 0x900ff,      0xff0, 0xff0 },	/* INT_MASK_REGS */	\
{ 0xffffffff, 0xffffffff, 0xfff }	/* ALL_REGS  */			\
}

/* The same information, inverted:
   Return the class number of the smallest class containing
   reg number REGNO.  This could be a conditional expression
   or could index an array.  */

#define REGNO_REG_CLASS(REGNO) (regclass_map[(REGNO)])

/* When this hook returns true for MODE, the compiler allows
   registers explicitly used in the rtl to be used as spill registers
   but prevents the compiler from extending the lifetime of these
   registers.  */
#define TARGET_SMALL_REGISTER_CLASSES_FOR_MODE_P hook_bool_mode_true

#define QI_REG_P(X) (REG_P (X) && QI_REGNO_P (REGNO (X)))
#define QI_REGNO_P(N) IN_RANGE ((N), FIRST_QI_REG, LAST_QI_REG)

#define LEGACY_INT_REG_P(X) (REG_P (X) && LEGACY_INT_REGNO_P (REGNO (X)))
#define LEGACY_INT_REGNO_P(N) (IN_RANGE ((N), FIRST_INT_REG, LAST_INT_REG))

#define REX_INT_REG_P(X) (REG_P (X) && REX_INT_REGNO_P (REGNO (X)))
#define REX_INT_REGNO_P(N) \
  IN_RANGE ((N), FIRST_REX_INT_REG, LAST_REX_INT_REG)

#define GENERAL_REG_P(X) (REG_P (X) && GENERAL_REGNO_P (REGNO (X)))
#define GENERAL_REGNO_P(N) \
  (LEGACY_INT_REGNO_P (N) || REX_INT_REGNO_P (N))

#define ANY_QI_REG_P(X) (REG_P (X) && ANY_QI_REGNO_P (REGNO (X)))
#define ANY_QI_REGNO_P(N) \
  (TARGET_64BIT ? GENERAL_REGNO_P (N) : QI_REGNO_P (N))

#define STACK_REG_P(X) (REG_P (X) && STACK_REGNO_P (REGNO (X)))
#define STACK_REGNO_P(N) IN_RANGE ((N), FIRST_STACK_REG, LAST_STACK_REG)

#define SSE_REG_P(X) (REG_P (X) && SSE_REGNO_P (REGNO (X)))
#define SSE_REGNO_P(N)						\
  (LEGACY_SSE_REGNO_P (N)					\
   || REX_SSE_REGNO_P (N)					\
   || EXT_REX_SSE_REGNO_P (N))

#define LEGACY_SSE_REGNO_P(N) \
  IN_RANGE ((N), FIRST_SSE_REG, LAST_SSE_REG)

#define REX_SSE_REGNO_P(N) \
  IN_RANGE ((N), FIRST_REX_SSE_REG, LAST_REX_SSE_REG)

#define EXT_REX_SSE_REG_P(X) (REG_P (X) && EXT_REX_SSE_REGNO_P (REGNO (X)))

#define EXT_REX_SSE_REGNO_P(N) \
  IN_RANGE ((N), FIRST_EXT_REX_SSE_REG, LAST_EXT_REX_SSE_REG)

#define ANY_FP_REG_P(X) (REG_P (X) && ANY_FP_REGNO_P (REGNO (X)))
#define ANY_FP_REGNO_P(N) (STACK_REGNO_P (N) || SSE_REGNO_P (N))

#define MASK_REG_P(X) (REG_P (X) && MASK_REGNO_P (REGNO (X)))
#define MASK_REGNO_P(N) IN_RANGE ((N), FIRST_MASK_REG, LAST_MASK_REG)
#define MASK_PAIR_REGNO_P(N) ((((N) - FIRST_MASK_REG) & 1) == 0)

#define MMX_REG_P(X) (REG_P (X) && MMX_REGNO_P (REGNO (X)))
#define MMX_REGNO_P(N) IN_RANGE ((N), FIRST_MMX_REG, LAST_MMX_REG)

#define CC_REG_P(X) (REG_P (X) && CC_REGNO_P (REGNO (X)))
#define CC_REGNO_P(X) ((X) == FLAGS_REG)

#define MOD4_SSE_REG_P(X) (REG_P (X) && MOD4_SSE_REGNO_P (REGNO (X)))
#define MOD4_SSE_REGNO_P(N) ((N) == XMM0_REG  \
			     || (N) == XMM4_REG  \
			     || (N) == XMM8_REG  \
			     || (N) == XMM12_REG \
			     || (N) == XMM16_REG \
			     || (N) == XMM20_REG \
			     || (N) == XMM24_REG \
			     || (N) == XMM28_REG)

/* First floating point reg */
#define FIRST_FLOAT_REG FIRST_STACK_REG
#define STACK_TOP_P(X) (REG_P (X) && REGNO (X) == FIRST_FLOAT_REG)

#define GET_SSE_REGNO(N)			\
  ((N) < 8 ? FIRST_SSE_REG + (N)		\
   : (N) < 16 ? FIRST_REX_SSE_REG + (N) - 8	\
   : FIRST_EXT_REX_SSE_REG + (N) - 16)

/* The class value for index registers, and the one for base regs.  */

#define INDEX_REG_CLASS INDEX_REGS
#define BASE_REG_CLASS GENERAL_REGS

/* Stack layout; function entry, exit and calling.  */

/* Define this if pushing a word on the stack
   makes the stack pointer a smaller address.  */
#define STACK_GROWS_DOWNWARD 1

/* Define this to nonzero if the nominal address of the stack frame
   is at the high-address end of the local variables;
   that is, each additional local variable allocated
   goes at a more negative offset in the frame.  */
#define FRAME_GROWS_DOWNWARD 1

#define PUSH_ROUNDING(BYTES) ix86_push_rounding (BYTES)

/* If defined, the maximum amount of space required for outgoing arguments
   will be computed and placed into the variable `crtl->outgoing_args_size'.
   No space will be pushed onto the stack for each call; instead, the
   function prologue should increase the stack frame size by this amount.  

   In 32bit mode enabling argument accumulation results in about 5% code size
   growth because move instructions are less compact than push.  In 64bit
   mode the difference is less drastic but visible.  

   FIXME: Unlike earlier implementations, the size of unwind info seems to
   actually grow with accumulation.  Is that because accumulated args
   unwind info became unnecesarily bloated?

   With the 64-bit MS ABI, we can generate correct code with or without
   accumulated args, but because of OUTGOING_REG_PARM_STACK_SPACE the code
   generated without accumulated args is terrible.

   If stack probes are required, the space used for large function
   arguments on the stack must also be probed, so enable
   -maccumulate-outgoing-args so this happens in the prologue.

   We must use argument accumulation in interrupt function if stack
   may be realigned to avoid DRAP.  */

#define ACCUMULATE_OUTGOING_ARGS \
  ((TARGET_ACCUMULATE_OUTGOING_ARGS \
    && optimize_function_for_speed_p (cfun)) \
   || (cfun->machine->func_type != TYPE_NORMAL \
       && crtl->stack_realign_needed) \
   || TARGET_STACK_PROBE \
   || TARGET_64BIT_MS_ABI \
   || (TARGET_MACHO && crtl->profile))

/* We want the stack and args grow in opposite directions, even if
   targetm.calls.push_argument returns false.  */
#define PUSH_ARGS_REVERSED 1

/* Offset of first parameter from the argument pointer register value.  */
#define FIRST_PARM_OFFSET(FNDECL) 0

/* Define this macro if functions should assume that stack space has been
   allocated for arguments even when their values are passed in registers.

   The value of this macro is the size, in bytes, of the area reserved for
   arguments passed in registers for the function represented by FNDECL.

   This space can be allocated by the caller, or be a part of the
   machine-dependent stack frame: `OUTGOING_REG_PARM_STACK_SPACE' says
   which.  */
#define REG_PARM_STACK_SPACE(FNDECL) ix86_reg_parm_stack_space (FNDECL)

#define OUTGOING_REG_PARM_STACK_SPACE(FNTYPE) \
  (TARGET_64BIT && ix86_function_type_abi (FNTYPE) == MS_ABI)

/* Define how to find the value returned by a library function
   assuming the value has mode MODE.  */

#define LIBCALL_VALUE(MODE) ix86_libcall_value (MODE)

/* Define the size of the result block used for communication between
   untyped_call and untyped_return.  The block contains a DImode value
   followed by the block used by fnsave and frstor.  */

#define APPLY_RESULT_SIZE (8+108)

/* 1 if N is a possible register number for function argument passing.  */
#define FUNCTION_ARG_REGNO_P(N) ix86_function_arg_regno_p (N)

/* Define a data type for recording info about an argument list
   during the scan of that argument list.  This data type should
   hold all necessary information about the function itself
   and about the args processed so far, enough to enable macros
   such as FUNCTION_ARG to determine where the next arg should go.  */

typedef struct ix86_args {
  int words;			/* # words passed so far */
  int nregs;			/* # registers available for passing */
  int regno;			/* next available register number */
  int fastcall;			/* fastcall or thiscall calling convention
				   is used */
  int sse_words;		/* # sse words passed so far */
  int sse_nregs;		/* # sse registers available for passing */
  int warn_avx512f;		/* True when we want to warn
				   about AVX512F ABI.  */
  int warn_avx;			/* True when we want to warn about AVX ABI.  */
  int warn_sse;			/* True when we want to warn about SSE ABI.  */
  int warn_mmx;			/* True when we want to warn about MMX ABI.  */
  int warn_empty;		/* True when we want to warn about empty classes
				   passing ABI change.  */
  int sse_regno;		/* next available sse register number */
  int mmx_words;		/* # mmx words passed so far */
  int mmx_nregs;		/* # mmx registers available for passing */
  int mmx_regno;		/* next available mmx register number */
  int maybe_vaarg;		/* true for calls to possibly vardic fncts.  */
  int caller;			/* true if it is caller.  */
  int float_in_sse;		/* Set to 1 or 2 for 32bit targets if
				   SFmode/DFmode arguments should be passed
				   in SSE registers.  Otherwise 0.  */
  int stdarg;                   /* Set to 1 if function is stdarg.  */
  enum calling_abi call_abi;	/* Set to SYSV_ABI for sysv abi. Otherwise
 				   MS_ABI for ms abi.  */
  tree decl;			/* Callee decl.  */
} CUMULATIVE_ARGS;

/* Initialize a variable CUM of type CUMULATIVE_ARGS
   for a call to a function whose data type is FNTYPE.
   For a library call, FNTYPE is 0.  */

#define INIT_CUMULATIVE_ARGS(CUM, FNTYPE, LIBNAME, FNDECL, N_NAMED_ARGS) \
  init_cumulative_args (&(CUM), (FNTYPE), (LIBNAME), (FNDECL), \
			(N_NAMED_ARGS) != -1)

/* Output assembler code to FILE to increment profiler label # LABELNO
   for profiling a function entry.  */

#define FUNCTION_PROFILER(FILE, LABELNO) \
  x86_function_profiler ((FILE), (LABELNO))

#define MCOUNT_NAME "_mcount"

#define MCOUNT_NAME_BEFORE_PROLOGUE "__fentry__"

#define PROFILE_COUNT_REGISTER "edx"

/* EXIT_IGNORE_STACK should be nonzero if, when returning from a function,
   the stack pointer does not matter.  The value is tested only in
   functions that have frame pointers.
   No definition is equivalent to always zero.  */
/* Note on the 386 it might be more efficient not to define this since
   we have to restore it ourselves from the frame pointer, in order to
   use pop */

#define EXIT_IGNORE_STACK 1

/* Define this macro as a C expression that is nonzero for registers
   used by the epilogue or the `return' pattern.  */

#define EPILOGUE_USES(REGNO) ix86_epilogue_uses (REGNO)

/* Output assembler code for a block containing the constant parts
   of a trampoline, leaving space for the variable parts.  */

/* On the 386, the trampoline contains two instructions:
     mov #STATIC,ecx
     jmp FUNCTION
   The trampoline is generated entirely at runtime.  The operand of JMP
   is the address of FUNCTION relative to the instruction following the
   JMP (which is 5 bytes long).  */

/* Length in units of the trampoline for entering a nested function.  */

#define TRAMPOLINE_SIZE (TARGET_64BIT ? 28 : 14)

/* Definitions for register eliminations.

   This is an array of structures.  Each structure initializes one pair
   of eliminable registers.  The "from" register number is given first,
   followed by "to".  Eliminations of the same "from" register are listed
   in order of preference.

   There are two registers that can always be eliminated on the i386.
   The frame pointer and the arg pointer can be replaced by either the
   hard frame pointer or to the stack pointer, depending upon the
   circumstances.  The hard frame pointer is not used before reload and
   so it is not eligible for elimination.  */

#define ELIMINABLE_REGS					\
{{ ARG_POINTER_REGNUM, STACK_POINTER_REGNUM},		\
 { ARG_POINTER_REGNUM, HARD_FRAME_POINTER_REGNUM},	\
 { FRAME_POINTER_REGNUM, STACK_POINTER_REGNUM},		\
 { FRAME_POINTER_REGNUM, HARD_FRAME_POINTER_REGNUM}}	\

/* Define the offset between two registers, one to be eliminated, and the other
   its replacement, at the start of a routine.  */

#define INITIAL_ELIMINATION_OFFSET(FROM, TO, OFFSET) \
  ((OFFSET) = ix86_initial_elimination_offset ((FROM), (TO)))

/* Addressing modes, and classification of registers for them.  */

/* Macros to check register numbers against specific register classes.  */

/* These assume that REGNO is a hard or pseudo reg number.
   They give nonzero only if REGNO is a hard reg of the suitable class
   or a pseudo reg currently allocated to a suitable hard reg.
   Since they use reg_renumber, they are safe only once reg_renumber
   has been allocated, which happens in reginfo.cc during register
   allocation.  */

#define REGNO_OK_FOR_INDEX_P(REGNO) 					\
  ((REGNO) < STACK_POINTER_REGNUM 					\
   || REX_INT_REGNO_P (REGNO)						\
   || (unsigned) reg_renumber[(REGNO)] < STACK_POINTER_REGNUM		\
   || REX_INT_REGNO_P ((unsigned) reg_renumber[(REGNO)]))

#define REGNO_OK_FOR_BASE_P(REGNO) 					\
  (GENERAL_REGNO_P (REGNO)						\
   || (REGNO) == ARG_POINTER_REGNUM 					\
   || (REGNO) == FRAME_POINTER_REGNUM 					\
   || GENERAL_REGNO_P ((unsigned) reg_renumber[(REGNO)]))

/* The macros REG_OK_FOR..._P assume that the arg is a REG rtx
   and check its validity for a certain class.
   We have two alternate definitions for each of them.
   The usual definition accepts all pseudo regs; the other rejects
   them unless they have been allocated suitable hard regs.
   The symbol REG_OK_STRICT causes the latter definition to be used.

   Most source files want to accept pseudo regs in the hope that
   they will get allocated to the class that the insn wants them to be in.
   Source files for reload pass need to be strict.
   After reload, it makes no difference, since pseudo regs have
   been eliminated by then.  */


/* Non strict versions, pseudos are ok.  */
#define REG_OK_FOR_INDEX_NONSTRICT_P(X)					\
  (REGNO (X) < STACK_POINTER_REGNUM					\
   || REX_INT_REGNO_P (REGNO (X))					\
   || REGNO (X) >= FIRST_PSEUDO_REGISTER)

#define REG_OK_FOR_BASE_NONSTRICT_P(X)					\
  (GENERAL_REGNO_P (REGNO (X))						\
   || REGNO (X) == ARG_POINTER_REGNUM					\
   || REGNO (X) == FRAME_POINTER_REGNUM 				\
   || REGNO (X) >= FIRST_PSEUDO_REGISTER)

/* Strict versions, hard registers only */
#define REG_OK_FOR_INDEX_STRICT_P(X) REGNO_OK_FOR_INDEX_P (REGNO (X))
#define REG_OK_FOR_BASE_STRICT_P(X)  REGNO_OK_FOR_BASE_P (REGNO (X))

#ifndef REG_OK_STRICT
#define REG_OK_FOR_INDEX_P(X)  REG_OK_FOR_INDEX_NONSTRICT_P (X)
#define REG_OK_FOR_BASE_P(X)   REG_OK_FOR_BASE_NONSTRICT_P (X)

#else
#define REG_OK_FOR_INDEX_P(X)  REG_OK_FOR_INDEX_STRICT_P (X)
#define REG_OK_FOR_BASE_P(X)   REG_OK_FOR_BASE_STRICT_P (X)
#endif

/* TARGET_LEGITIMATE_ADDRESS_P recognizes an RTL expression
   that is a valid memory address for an instruction.
   The MODE argument is the machine mode for the MEM expression
   that wants to use this address.

   The other macros defined here are used only in TARGET_LEGITIMATE_ADDRESS_P,
   except for CONSTANT_ADDRESS_P which is usually machine-independent.

   See legitimize_pic_address in i386.cc for details as to what
   constitutes a legitimate address when -fpic is used.  */

#define MAX_REGS_PER_ADDRESS 2

#define CONSTANT_ADDRESS_P(X)  constant_address_p (X)

/* If defined, a C expression to determine the base term of address X.
   This macro is used in only one place: `find_base_term' in alias.cc.

   It is always safe for this macro to not be defined.  It exists so
   that alias analysis can understand machine-dependent addresses.

   The typical use of this macro is to handle addresses containing
   a label_ref or symbol_ref within an UNSPEC.  */

#define FIND_BASE_TERM(X) ix86_find_base_term (X)

/* Nonzero if the constant value X is a legitimate general operand
   when generating PIC code.  It is given that flag_pic is on and
   that X satisfies CONSTANT_P or is a CONST_DOUBLE.  */

#define LEGITIMATE_PIC_OPERAND_P(X) legitimate_pic_operand_p (X)

#define STRIP_UNARY(X) (UNARY_P (X) ? XEXP (X, 0) : X)

#define SYMBOLIC_CONST(X)	\
  (GET_CODE (X) == SYMBOL_REF						\
   || GET_CODE (X) == LABEL_REF						\
   || (GET_CODE (X) == CONST && symbolic_reference_mentioned_p (X)))

/* Max number of args passed in registers.  If this is more than 3, we will
   have problems with ebx (register #4), since it is a caller save register and
   is also used as the pic register in ELF.  So for now, don't allow more than
   3 registers to be passed in registers.  */

/* Abi specific values for REGPARM_MAX and SSE_REGPARM_MAX */
#define X86_64_REGPARM_MAX 6
#define X86_64_MS_REGPARM_MAX 4

#define X86_32_REGPARM_MAX 3

#define REGPARM_MAX							\
  (TARGET_64BIT								\
   ? (TARGET_64BIT_MS_ABI						\
      ? X86_64_MS_REGPARM_MAX						\
      : X86_64_REGPARM_MAX)						\
   : X86_32_REGPARM_MAX)

#define X86_64_SSE_REGPARM_MAX 8
#define X86_64_MS_SSE_REGPARM_MAX 4

#define X86_32_SSE_REGPARM_MAX (TARGET_SSE ? (TARGET_MACHO ? 4 : 3) : 0)

#define SSE_REGPARM_MAX							\
  (TARGET_64BIT								\
   ? (TARGET_64BIT_MS_ABI						\
      ? X86_64_MS_SSE_REGPARM_MAX					\
      : X86_64_SSE_REGPARM_MAX)						\
   : X86_32_SSE_REGPARM_MAX)

#define X86_32_MMX_REGPARM_MAX (TARGET_MMX ? (TARGET_MACHO ? 0 : 3) : 0)

#define MMX_REGPARM_MAX (TARGET_64BIT ? 0 : X86_32_MMX_REGPARM_MAX)

/* Specify the machine mode that this machine uses
   for the index in the tablejump instruction.  */
#define CASE_VECTOR_MODE \
 (!TARGET_LP64 || (flag_pic && ix86_cmodel != CM_LARGE_PIC) ? SImode : DImode)

/* Define this as 1 if `char' should by default be signed; else as 0.  */
#define DEFAULT_SIGNED_CHAR 1

/* The constant maximum number of bytes that a single instruction can
   move quickly between memory and registers or between two memory
   locations.  */
#define MAX_MOVE_MAX 64

/* Max number of bytes we can move from memory to memory in one
   reasonably fast instruction, as opposed to MOVE_MAX_PIECES which
   is the number of bytes at a time which we can move efficiently.
   MOVE_MAX_PIECES defaults to MOVE_MAX.  */

#define MOVE_MAX \
  ((TARGET_AVX512F \
    && (ix86_move_max == PVW_AVX512 \
	|| ix86_store_max == PVW_AVX512)) \
   ? 64 \
   : ((TARGET_AVX \
       && (ix86_move_max >= PVW_AVX256 \
	   || ix86_store_max >= PVW_AVX256)) \
      ? 32 \
      : ((TARGET_SSE2 \
	  && TARGET_SSE_UNALIGNED_LOAD_OPTIMAL \
	  && TARGET_SSE_UNALIGNED_STORE_OPTIMAL) \
	 ? 16 : UNITS_PER_WORD)))

/* STORE_MAX_PIECES is the number of bytes at a time that we can store
   efficiently.  Allow 16/32/64 bytes only if inter-unit move is enabled
   since vec_duplicate enabled by inter-unit move is used to implement
   store_by_pieces of 16/32/64 bytes.  */
#define STORE_MAX_PIECES \
  (TARGET_INTER_UNIT_MOVES_TO_VEC \
   ? ((TARGET_AVX512F && ix86_store_max == PVW_AVX512) \
      ? 64 \
      : ((TARGET_AVX \
	  && ix86_store_max >= PVW_AVX256) \
	  ? 32 \
	  : ((TARGET_SSE2 \
	      && TARGET_SSE_UNALIGNED_STORE_OPTIMAL) \
	      ? 16 : UNITS_PER_WORD))) \
   : UNITS_PER_WORD)

/* If a memory-to-memory move would take MOVE_RATIO or more simple
   move-instruction pairs, we will do a cpymem or libcall instead.
   Increasing the value will always make code faster, but eventually
   incurs high cost in increased code size.

   If you don't define this, a reasonable default is used.  */

#define MOVE_RATIO(speed) ((speed) ? ix86_cost->move_ratio : 3)

/* If a clear memory operation would take CLEAR_RATIO or more simple
   move-instruction sequences, we will do a clrmem or libcall instead.  */

#define CLEAR_RATIO(speed) ((speed) ? ix86_cost->clear_ratio : 2)

/* Define if shifts truncate the shift count which implies one can
   omit a sign-extension or zero-extension of a shift count.

   On i386, shifts do truncate the count.  But bit test instructions
   take the modulo of the bit offset operand.  */

/* #define SHIFT_COUNT_TRUNCATED */

/* A macro to update M and UNSIGNEDP when an object whose type is
   TYPE and which has the specified mode and signedness is to be
   stored in a register.  This macro is only called when TYPE is a
   scalar type.

   On i386 it is sometimes useful to promote HImode and QImode
   quantities to SImode.  The choice depends on target type.  */

#define PROMOTE_MODE(MODE, UNSIGNEDP, TYPE) 		\
do {							\
  if (((MODE) == HImode && TARGET_PROMOTE_HI_REGS)	\
      || ((MODE) == QImode && TARGET_PROMOTE_QI_REGS))	\
    (MODE) = SImode;					\
} while (0)

/* Specify the machine mode that pointers have.
   After generation of rtl, the compiler makes no further distinction
   between pointers and any other objects of this machine mode.  */
#define Pmode (ix86_pmode == PMODE_DI ? DImode : SImode)

/* Supply a definition of STACK_SAVEAREA_MODE for emit_stack_save.
   NONLOCAL needs space to save both shadow stack and stack pointers.

   FIXME: We only need to save and restore stack pointer in ptr_mode.
   But expand_builtin_setjmp_setup and expand_builtin_longjmp use Pmode
   to save and restore stack pointer.  See
   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=84150
 */
#define STACK_SAVEAREA_MODE(LEVEL)			\
  ((LEVEL) == SAVE_NONLOCAL ? (TARGET_64BIT ? TImode : DImode) : Pmode)

/* Specify the machine_mode of the size increment
   operand of an 'allocate_stack' named pattern.  */
#define STACK_SIZE_MODE Pmode

/* A C expression whose value is zero if pointers that need to be extended
   from being `POINTER_SIZE' bits wide to `Pmode' are sign-extended and
   greater then zero if they are zero-extended and less then zero if the
   ptr_extend instruction should be used.  */

#define POINTERS_EXTEND_UNSIGNED 1

/* A function address in a call instruction
   is a byte address (for indexing purposes)
   so give the MEM rtx a byte's mode.  */
#define FUNCTION_MODE QImode


/* A C expression for the cost of a branch instruction.  A value of 1
   is the default; other values are interpreted relative to that.  */

#define BRANCH_COST(speed_p, predictable_p) \
  (!(speed_p) ? 2 : (predictable_p) ? 0 : ix86_branch_cost)

/* An integer expression for the size in bits of the largest integer machine
   mode that should actually be used.  We allow pairs of registers.  */
#define MAX_FIXED_MODE_SIZE GET_MODE_BITSIZE (TARGET_64BIT ? TImode : DImode)

/* Define this macro as a C expression which is nonzero if accessing
   less than a word of memory (i.e. a `char' or a `short') is no
   faster than accessing a word of memory, i.e., if such access
   require more than one instruction or if there is no difference in
   cost between byte and (aligned) word loads.

   When this macro is not defined, the compiler will access a field by
   finding the smallest containing object; when it is defined, a
   fullword load will be used if alignment permits.  Unless bytes
   accesses are faster than word accesses, using word accesses is
   preferable since it may eliminate subsequent memory access if
   subsequent accesses occur to other fields in the same word of the
   structure, but to different bytes.  */

#define SLOW_BYTE_ACCESS 0

/* Nonzero if access to memory by shorts is slow and undesirable.  */
#define SLOW_SHORT_ACCESS 0

/* Define this macro if it is as good or better to call a constant
   function address than to call an address kept in a register.

   Desirable on the 386 because a CALL with a constant address is
   faster than one with a register address.  */

#define NO_FUNCTION_CSE 1

/* Given a comparison code (EQ, NE, etc.) and the first operand of a COMPARE,
   return the mode to be used for the comparison.

   For floating-point equality comparisons, CCFPEQmode should be used.
   VOIDmode should be used in all other cases.

   For integer comparisons against zero, reduce to CCNOmode or CCZmode if
   possible, to allow for more combinations.  */

#define SELECT_CC_MODE(OP, X, Y) ix86_cc_mode ((OP), (X), (Y))

/* Return nonzero if MODE implies a floating point inequality can be
   reversed.  */

#define REVERSIBLE_CC_MODE(MODE) 1

/* A C expression whose value is reversed condition code of the CODE for
   comparison done in CC_MODE mode.  */
#define REVERSE_CONDITION(CODE, MODE) ix86_reverse_condition ((CODE), (MODE))


/* Control the assembler format that we output, to the extent
   this does not vary between assemblers.  */

/* How to refer to registers in assembler output.
   This sequence is indexed by compiler's hard-register-number (see above).  */

/* In order to refer to the first 8 regs as 32-bit regs, prefix an "e".
   For non floating point regs, the following are the HImode names.

   For float regs, the stack top is sometimes referred to as "%st(0)"
   instead of just "%st".  TARGET_PRINT_OPERAND handles this with the
   "y" code.  */

#define HI_REGISTER_NAMES						\
{"ax","dx","cx","bx","si","di","bp","sp",				\
 "st","st(1)","st(2)","st(3)","st(4)","st(5)","st(6)","st(7)",		\
 "argp", "flags", "fpsr", "frame",					\
 "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7",		\
 "mm0", "mm1", "mm2", "mm3", "mm4", "mm5", "mm6", "mm7",		\
 "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15",			\
 "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15",	\
 "xmm16", "xmm17", "xmm18", "xmm19",					\
 "xmm20", "xmm21", "xmm22", "xmm23",					\
 "xmm24", "xmm25", "xmm26", "xmm27",					\
 "xmm28", "xmm29", "xmm30", "xmm31",					\
 "k0", "k1", "k2", "k3", "k4", "k5", "k6", "k7" }

#define REGISTER_NAMES HI_REGISTER_NAMES

#define QI_REGISTER_NAMES \
{"al", "dl", "cl", "bl", "sil", "dil", "bpl", "spl"}

#define QI_HIGH_REGISTER_NAMES \
{"ah", "dh", "ch", "bh"}

/* Table of additional register names to use in user input.  */

#define ADDITIONAL_REGISTER_NAMES						\
{										\
  { "eax", AX_REG }, { "edx", DX_REG }, { "ecx", CX_REG }, { "ebx", BX_REG },	\
  { "esi", SI_REG }, { "edi", DI_REG }, { "ebp", BP_REG }, { "esp", SP_REG },	\
  { "rax", AX_REG }, { "rdx", DX_REG }, { "rcx", CX_REG }, { "rbx", BX_REG },	\
  { "rsi", SI_REG }, { "rdi", DI_REG }, { "rbp", BP_REG }, { "rsp", SP_REG },	\
  { "al", AX_REG }, { "dl", DX_REG }, { "cl", CX_REG }, { "bl", BX_REG },	\
  { "sil", SI_REG }, { "dil", DI_REG }, { "bpl", BP_REG }, { "spl", SP_REG },	\
  { "ah", AX_REG }, { "dh", DX_REG }, { "ch", CX_REG }, { "bh", BX_REG },	\
  { "ymm0", XMM0_REG }, { "ymm1", XMM1_REG }, { "ymm2", XMM2_REG }, { "ymm3", XMM3_REG }, \
  { "ymm4", XMM4_REG }, { "ymm5", XMM5_REG }, { "ymm6", XMM6_REG }, { "ymm7", XMM7_REG }, \
  { "ymm8", XMM8_REG }, { "ymm9", XMM9_REG }, { "ymm10", XMM10_REG }, { "ymm11", XMM11_REG }, \
  { "ymm12", XMM12_REG }, { "ymm13", XMM13_REG }, { "ymm14", XMM14_REG }, { "ymm15", XMM15_REG }, \
  { "ymm16", XMM16_REG }, { "ymm17", XMM17_REG }, { "ymm18", XMM18_REG }, { "ymm19", XMM19_REG }, \
  { "ymm20", XMM20_REG }, { "ymm21", XMM21_REG }, { "ymm22", XMM22_REG }, { "ymm23", XMM23_REG }, \
  { "ymm24", XMM24_REG }, { "ymm25", XMM25_REG }, { "ymm26", XMM26_REG }, { "ymm27", XMM27_REG }, \
  { "ymm28", XMM28_REG }, { "ymm29", XMM29_REG }, { "ymm30", XMM30_REG }, { "ymm31", XMM31_REG }, \
  { "zmm0", XMM0_REG }, { "zmm1", XMM1_REG }, { "zmm2", XMM2_REG }, { "zmm3", XMM3_REG }, \
  { "zmm4", XMM4_REG }, { "zmm5", XMM5_REG }, { "zmm6", XMM6_REG }, { "zmm7", XMM7_REG }, \
  { "zmm8", XMM8_REG }, { "zmm9", XMM9_REG }, { "zmm10", XMM10_REG }, { "zmm11", XMM11_REG }, \
  { "zmm12", XMM12_REG }, { "zmm13", XMM13_REG }, { "zmm14", XMM14_REG }, { "zmm15", XMM15_REG }, \
  { "zmm16", XMM16_REG }, { "zmm17", XMM17_REG }, { "zmm18", XMM18_REG }, { "zmm19", XMM19_REG }, \
  { "zmm20", XMM20_REG }, { "zmm21", XMM21_REG }, { "zmm22", XMM22_REG }, { "zmm23", XMM23_REG }, \
  { "zmm24", XMM24_REG }, { "zmm25", XMM25_REG }, { "zmm26", XMM26_REG }, { "zmm27", XMM27_REG }, \
  { "zmm28", XMM28_REG }, { "zmm29", XMM29_REG }, { "zmm30", XMM30_REG }, { "zmm31", XMM31_REG }  \
}

/* How to renumber registers for dbx and gdb.  */

#define DBX_REGISTER_NUMBER(N) \
  (TARGET_64BIT ? dbx64_register_map[(N)] : dbx_register_map[(N)])

extern int const dbx_register_map[FIRST_PSEUDO_REGISTER];
extern int const dbx64_register_map[FIRST_PSEUDO_REGISTER];
extern int const svr4_dbx_register_map[FIRST_PSEUDO_REGISTER];

/* Before the prologue, RA is at 0(%esp).  */
#define INCOMING_RETURN_ADDR_RTX \
  gen_rtx_MEM (Pmode, stack_pointer_rtx)

/* After the prologue, RA is at -4(AP) in the current frame.  */
#define RETURN_ADDR_RTX(COUNT, FRAME)					\
  ((COUNT) == 0								\
   ? gen_rtx_MEM (Pmode, plus_constant (Pmode, arg_pointer_rtx,		\
					-UNITS_PER_WORD))		\
   : gen_rtx_MEM (Pmode, plus_constant (Pmode, (FRAME), UNITS_PER_WORD)))

/* PC is dbx register 8; let's use that column for RA.  */
#define DWARF_FRAME_RETURN_COLUMN 	(TARGET_64BIT ? 16 : 8)

/* Before the prologue, there are return address and error code for
   exception handler on the top of the frame.  */
#define INCOMING_FRAME_SP_OFFSET \
  (cfun->machine->func_type == TYPE_EXCEPTION \
   ? 2 * UNITS_PER_WORD : UNITS_PER_WORD)

/* The value of INCOMING_FRAME_SP_OFFSET the assembler assumes in
   .cfi_startproc.  */
#define DEFAULT_INCOMING_FRAME_SP_OFFSET UNITS_PER_WORD

/* Describe how we implement __builtin_eh_return.  */
#define EH_RETURN_DATA_REGNO(N)	((N) <= DX_REG ? (N) : INVALID_REGNUM)
#define EH_RETURN_STACKADJ_RTX	gen_rtx_REG (Pmode, CX_REG)


/* Select a format to encode pointers in exception handling data.  CODE
   is 0 for data, 1 for code labels, 2 for function pointers.  GLOBAL is
   true if the symbol may be affected by dynamic relocations.

   ??? All x86 object file formats are capable of representing this.
   After all, the relocation needed is the same as for the call insn.
   Whether or not a particular assembler allows us to enter such, I
   guess we'll have to see.  */
#define ASM_PREFERRED_EH_DATA_FORMAT(CODE, GLOBAL)       		\
  asm_preferred_eh_data_format ((CODE), (GLOBAL))

/* These are a couple of extensions to the formats accepted
   by asm_fprintf:
     %z prints out opcode suffix for word-mode instruction
     %r prints out word-mode name for reg_names[arg]  */
#define ASM_FPRINTF_EXTENSIONS(FILE, ARGS, P)		\
  case 'z':						\
    fputc (TARGET_64BIT ? 'q' : 'l', (FILE));		\
    break;						\
							\
  case 'r':						\
    {							\
      unsigned int regno = va_arg ((ARGS), int);	\
      if (LEGACY_INT_REGNO_P (regno))			\
	fputc (TARGET_64BIT ? 'r' : 'e', (FILE));	\
      fputs (reg_names[regno], (FILE));			\
      break;						\
    }

/* This is how to output an insn to push a register on the stack.  */

#define ASM_OUTPUT_REG_PUSH(FILE, REGNO)		\
  asm_fprintf ((FILE), "\tpush%z\t%%%r\n", (REGNO))

/* This is how to output an insn to pop a register from the stack.  */

#define ASM_OUTPUT_REG_POP(FILE, REGNO)  \
  asm_fprintf ((FILE), "\tpop%z\t%%%r\n", (REGNO))

/* This is how to output an element of a case-vector that is absolute.  */

#define ASM_OUTPUT_ADDR_VEC_ELT(FILE, VALUE)  \
  ix86_output_addr_vec_elt ((FILE), (VALUE))

/* This is how to output an element of a case-vector that is relative.  */

#define ASM_OUTPUT_ADDR_DIFF_ELT(FILE, BODY, VALUE, REL) \
  ix86_output_addr_diff_elt ((FILE), (VALUE), (REL))

/* When we see %v, we will print the 'v' prefix if TARGET_AVX is true.  */

#define ASM_OUTPUT_AVX_PREFIX(STREAM, PTR)	\
{						\
  if ((PTR)[0] == '%' && (PTR)[1] == 'v')	\
    (PTR) += TARGET_AVX ? 1 : 2;		\
}

/* A C statement or statements which output an assembler instruction
   opcode to the stdio stream STREAM.  The macro-operand PTR is a
   variable of type `char *' which points to the opcode name in
   its "internal" form--the form that is written in the machine
   description.  */

#define ASM_OUTPUT_OPCODE(STREAM, PTR) \
  ASM_OUTPUT_AVX_PREFIX ((STREAM), (PTR))

/* A C statement to output to the stdio stream FILE an assembler
   command to pad the location counter to a multiple of 1<<LOG
   bytes if it is within MAX_SKIP bytes.  */

#ifdef HAVE_GAS_MAX_SKIP_P2ALIGN
# define ASM_OUTPUT_MAX_SKIP_ALIGN(FILE,LOG,MAX_SKIP)			\
  do {									\
    if ((LOG) != 0) {							\
      if ((MAX_SKIP) == 0 || (MAX_SKIP) >= (1 << (LOG)) - 1)		\
	fprintf ((FILE), "\t.p2align %d\n", (LOG));			\
      else								\
	fprintf ((FILE), "\t.p2align %d,,%d\n", (LOG), (MAX_SKIP));	\
    }									\
  } while (0)
#endif

/* Write the extra assembler code needed to declare a function
   properly.  */

#undef ASM_OUTPUT_FUNCTION_LABEL
#define ASM_OUTPUT_FUNCTION_LABEL(FILE, NAME, DECL) \
  ix86_asm_output_function_label ((FILE), (NAME), (DECL))

/* A C statement (sans semicolon) to output a reference to SYMBOL_REF SYM.
   If not defined, assemble_name will be used to output the name of the
   symbol.  This macro may be used to modify the way a symbol is referenced
   depending on information encoded by TARGET_ENCODE_SECTION_INFO.  */

#ifndef ASM_OUTPUT_SYMBOL_REF
#define ASM_OUTPUT_SYMBOL_REF(FILE, SYM) \
  do {							\
    const char *name					\
      = assemble_name_resolve (XSTR (x, 0));		\
    /* In -masm=att wrap identifiers that start with $	\
       into parens.  */					\
    if (ASSEMBLER_DIALECT == ASM_ATT			\
	&& name[0] == '$'				\
	&& user_label_prefix[0] == '\0')		\
      {							\
	fputc ('(', (FILE));				\
	assemble_name_raw ((FILE), name);		\
	fputc (')', (FILE));				\
      }							\
    else						\
      assemble_name_raw ((FILE), name);			\
  } while (0)
#endif

/* Under some conditions we need jump tables in the text section,
   because the assembler cannot handle label differences between
   sections.  */

#define JUMP_TABLES_IN_TEXT_SECTION \
  (flag_pic && !(TARGET_64BIT || HAVE_AS_GOTOFF_IN_DATA))

/* Switch to init or fini section via SECTION_OP, emit a call to FUNC,
   and switch back.  For x86 we do this only to save a few bytes that
   would otherwise be unused in the text section.  */
#define CRT_MKSTR2(VAL) #VAL
#define CRT_MKSTR(x) CRT_MKSTR2(x)

#define CRT_CALL_STATIC_FUNCTION(SECTION_OP, FUNC)		\
   asm (SECTION_OP "\n\t"					\
	"call " CRT_MKSTR(__USER_LABEL_PREFIX__) #FUNC "\n"	\
	TEXT_SECTION_ASM_OP);

/* Default threshold for putting data in large sections
   with x86-64 medium memory model */
#define DEFAULT_LARGE_SECTION_THRESHOLD 65536

/* Which processor to tune code generation for.  These must be in sync
   with processor_target_table in i386.cc.  */ 

enum processor_type
{
  PROCESSOR_GENERIC = 0,
  PROCESSOR_I386,			/* 80386 */
  PROCESSOR_I486,			/* 80486DX, 80486SX, 80486DX[24] */
  PROCESSOR_PENTIUM,
  PROCESSOR_LAKEMONT,
  PROCESSOR_PENTIUMPRO,
  PROCESSOR_PENTIUM4,
  PROCESSOR_NOCONA,
  PROCESSOR_CORE2,
  PROCESSOR_NEHALEM,
  PROCESSOR_SANDYBRIDGE,
  PROCESSOR_HASWELL,
  PROCESSOR_BONNELL,
  PROCESSOR_SILVERMONT,
  PROCESSOR_GOLDMONT,
  PROCESSOR_GOLDMONT_PLUS,
  PROCESSOR_TREMONT,
  PROCESSOR_KNL,
  PROCESSOR_KNM,
  PROCESSOR_SKYLAKE,
  PROCESSOR_SKYLAKE_AVX512,
  PROCESSOR_CANNONLAKE,
  PROCESSOR_ICELAKE_CLIENT,
  PROCESSOR_ICELAKE_SERVER,
  PROCESSOR_CASCADELAKE,
  PROCESSOR_TIGERLAKE,
  PROCESSOR_COOPERLAKE,
  PROCESSOR_SAPPHIRERAPIDS,
  PROCESSOR_ALDERLAKE,
  PROCESSOR_ROCKETLAKE,
  PROCESSOR_INTEL,
  PROCESSOR_GEODE,
  PROCESSOR_K6,
  PROCESSOR_ATHLON,
  PROCESSOR_K8,
  PROCESSOR_AMDFAM10,
  PROCESSOR_BDVER1,
  PROCESSOR_BDVER2,
  PROCESSOR_BDVER3,
  PROCESSOR_BDVER4,
  PROCESSOR_BTVER1,
  PROCESSOR_BTVER2,
  PROCESSOR_ZNVER1,
  PROCESSOR_ZNVER2,
  PROCESSOR_ZNVER3,
  PROCESSOR_ZNVER4,
  PROCESSOR_max
};

#if !defined(IN_LIBGCC2) && !defined(IN_TARGET_LIBS) && !defined(IN_RTS)
extern const char *const processor_names[];

#include "wide-int-bitmask.h"

enum pta_flag
{
#define DEF_PTA(NAME) _ ## NAME,
#include "i386-isa.def"
#undef DEF_PTA
  END_PTA
};

/* wide_int_bitmask can handle only 128 flags.  */
STATIC_ASSERT (END_PTA <= 128);

#define WIDE_INT_BITMASK_FROM_NTH(N) (N < 64 ? wide_int_bitmask (0, 1ULL << N) \
				      : wide_int_bitmask (1ULL << (N - 64), 0))

#define DEF_PTA(NAME) constexpr wide_int_bitmask PTA_ ## NAME \
  = WIDE_INT_BITMASK_FROM_NTH ((pta_flag) _ ## NAME);
#include "i386-isa.def"
#undef DEF_PTA

constexpr wide_int_bitmask PTA_X86_64_BASELINE = PTA_64BIT | PTA_MMX | PTA_SSE
  | PTA_SSE2 | PTA_NO_SAHF | PTA_FXSR;
constexpr wide_int_bitmask PTA_X86_64_V2 = (PTA_X86_64_BASELINE
					    & (~PTA_NO_SAHF))
  | PTA_CX16 | PTA_POPCNT | PTA_SSE3 | PTA_SSE4_1 | PTA_SSE4_2 | PTA_SSSE3;
constexpr wide_int_bitmask PTA_X86_64_V3 = PTA_X86_64_V2
  | PTA_AVX | PTA_AVX2 | PTA_BMI | PTA_BMI2 | PTA_F16C | PTA_FMA | PTA_LZCNT
  | PTA_MOVBE | PTA_XSAVE;
constexpr wide_int_bitmask PTA_X86_64_V4 = PTA_X86_64_V3
  | PTA_AVX512F | PTA_AVX512BW | PTA_AVX512CD | PTA_AVX512DQ | PTA_AVX512VL;

constexpr wide_int_bitmask PTA_CORE2 = PTA_64BIT | PTA_MMX | PTA_SSE | PTA_SSE2
  | PTA_SSE3 | PTA_SSSE3 | PTA_CX16 | PTA_FXSR;
constexpr wide_int_bitmask PTA_NEHALEM = PTA_CORE2 | PTA_SSE4_1 | PTA_SSE4_2
  | PTA_POPCNT;
constexpr wide_int_bitmask PTA_WESTMERE = PTA_NEHALEM | PTA_PCLMUL;
constexpr wide_int_bitmask PTA_SANDYBRIDGE = PTA_WESTMERE | PTA_AVX | PTA_XSAVE
  | PTA_XSAVEOPT;
constexpr wide_int_bitmask PTA_IVYBRIDGE = PTA_SANDYBRIDGE | PTA_FSGSBASE
  | PTA_RDRND | PTA_F16C;
constexpr wide_int_bitmask PTA_HASWELL = PTA_IVYBRIDGE | PTA_AVX2 | PTA_BMI
  | PTA_BMI2 | PTA_LZCNT | PTA_FMA | PTA_MOVBE | PTA_HLE;
constexpr wide_int_bitmask PTA_BROADWELL = PTA_HASWELL | PTA_ADX | PTA_RDSEED
  | PTA_PRFCHW;
constexpr wide_int_bitmask PTA_SKYLAKE = PTA_BROADWELL | PTA_AES
  | PTA_CLFLUSHOPT | PTA_XSAVEC | PTA_XSAVES | PTA_SGX;
constexpr wide_int_bitmask PTA_SKYLAKE_AVX512 = PTA_SKYLAKE | PTA_AVX512F
  | PTA_AVX512CD | PTA_AVX512VL | PTA_AVX512BW | PTA_AVX512DQ | PTA_PKU
  | PTA_CLWB;
constexpr wide_int_bitmask PTA_CASCADELAKE = PTA_SKYLAKE_AVX512
  | PTA_AVX512VNNI;
constexpr wide_int_bitmask PTA_COOPERLAKE = PTA_CASCADELAKE | PTA_AVX512BF16;
constexpr wide_int_bitmask PTA_CANNONLAKE = PTA_SKYLAKE | PTA_AVX512F
  | PTA_AVX512CD | PTA_AVX512VL | PTA_AVX512BW | PTA_AVX512DQ | PTA_PKU
  | PTA_AVX512VBMI | PTA_AVX512IFMA | PTA_SHA;
constexpr wide_int_bitmask PTA_ICELAKE_CLIENT = PTA_CANNONLAKE | PTA_AVX512VNNI
  | PTA_GFNI | PTA_VAES | PTA_AVX512VBMI2 | PTA_VPCLMULQDQ | PTA_AVX512BITALG
  | PTA_RDPID | PTA_AVX512VPOPCNTDQ;
constexpr wide_int_bitmask PTA_ROCKETLAKE = PTA_ICELAKE_CLIENT & ~PTA_SGX;
constexpr wide_int_bitmask PTA_ICELAKE_SERVER = PTA_ICELAKE_CLIENT
  | PTA_PCONFIG | PTA_WBNOINVD | PTA_CLWB;
constexpr wide_int_bitmask PTA_TIGERLAKE = PTA_ICELAKE_CLIENT | PTA_MOVDIRI
  | PTA_MOVDIR64B | PTA_CLWB | PTA_AVX512VP2INTERSECT | PTA_KL | PTA_WIDEKL;
constexpr wide_int_bitmask PTA_SAPPHIRERAPIDS = PTA_ICELAKE_SERVER | PTA_MOVDIRI
  | PTA_MOVDIR64B | PTA_ENQCMD | PTA_CLDEMOTE | PTA_PTWRITE | PTA_WAITPKG
  | PTA_SERIALIZE | PTA_TSXLDTRK | PTA_AMX_TILE | PTA_AMX_INT8 | PTA_AMX_BF16
  | PTA_UINTR | PTA_AVXVNNI | PTA_AVX512FP16 | PTA_AVX512BF16;
constexpr wide_int_bitmask PTA_KNL = PTA_BROADWELL | PTA_AVX512PF
  | PTA_AVX512ER | PTA_AVX512F | PTA_AVX512CD | PTA_PREFETCHWT1;
constexpr wide_int_bitmask PTA_BONNELL = PTA_CORE2 | PTA_MOVBE;
constexpr wide_int_bitmask PTA_SILVERMONT = PTA_WESTMERE | PTA_MOVBE
  | PTA_RDRND | PTA_PRFCHW;
constexpr wide_int_bitmask PTA_GOLDMONT = PTA_SILVERMONT | PTA_AES | PTA_SHA
  | PTA_XSAVE | PTA_RDSEED | PTA_XSAVEC | PTA_XSAVES | PTA_CLFLUSHOPT
  | PTA_XSAVEOPT | PTA_FSGSBASE;
constexpr wide_int_bitmask PTA_GOLDMONT_PLUS = PTA_GOLDMONT | PTA_RDPID
  | PTA_SGX | PTA_PTWRITE;
constexpr wide_int_bitmask PTA_TREMONT = PTA_GOLDMONT_PLUS | PTA_CLWB
  | PTA_GFNI | PTA_MOVDIRI | PTA_MOVDIR64B | PTA_CLDEMOTE | PTA_WAITPKG;
constexpr wide_int_bitmask PTA_ALDERLAKE = PTA_TREMONT | PTA_ADX | PTA_AVX
  | PTA_AVX2 | PTA_BMI | PTA_BMI2 | PTA_F16C | PTA_FMA | PTA_LZCNT
  | PTA_PCONFIG | PTA_PKU | PTA_VAES | PTA_VPCLMULQDQ | PTA_SERIALIZE
  | PTA_HRESET | PTA_KL | PTA_WIDEKL | PTA_AVXVNNI;
constexpr wide_int_bitmask PTA_KNM = PTA_KNL | PTA_AVX5124VNNIW
  | PTA_AVX5124FMAPS | PTA_AVX512VPOPCNTDQ;
constexpr wide_int_bitmask PTA_ZNVER1 = PTA_64BIT | PTA_MMX | PTA_SSE | PTA_SSE2
  | PTA_SSE3 | PTA_SSE4A | PTA_CX16 | PTA_ABM | PTA_SSSE3 | PTA_SSE4_1
  | PTA_SSE4_2 | PTA_AES | PTA_PCLMUL | PTA_AVX | PTA_AVX2 | PTA_BMI | PTA_BMI2
  | PTA_F16C | PTA_FMA | PTA_PRFCHW | PTA_FXSR | PTA_XSAVE | PTA_XSAVEOPT
  | PTA_FSGSBASE | PTA_RDRND | PTA_MOVBE | PTA_MWAITX | PTA_ADX | PTA_RDSEED
  | PTA_CLZERO | PTA_CLFLUSHOPT | PTA_XSAVEC | PTA_XSAVES | PTA_SHA | PTA_LZCNT
  | PTA_POPCNT;
constexpr wide_int_bitmask PTA_ZNVER2 = PTA_ZNVER1 | PTA_CLWB | PTA_RDPID
  | PTA_WBNOINVD;
constexpr wide_int_bitmask PTA_ZNVER3 = PTA_ZNVER2 | PTA_VAES | PTA_VPCLMULQDQ
  | PTA_PKU;
constexpr wide_int_bitmask PTA_ZNVER4 = PTA_ZNVER3 | PTA_AVX512F | PTA_AVX512DQ
  | PTA_AVX512IFMA | PTA_AVX512CD | PTA_AVX512BW | PTA_AVX512VL
  | PTA_AVX512BF16 | PTA_AVX512VBMI | PTA_AVX512VBMI2 | PTA_GFNI
  | PTA_AVX512VNNI | PTA_AVX512BITALG | PTA_AVX512VPOPCNTDQ;

#ifndef GENERATOR_FILE

#include "insn-attr-common.h"

#include "common/config/i386/i386-cpuinfo.h"

class pta
{
public:
  const char *const name;		/* processor name or nickname.  */
  const enum processor_type processor;
  const enum attr_cpu schedule;
  const wide_int_bitmask flags;
  const int model;
  const enum feature_priority priority;
};

extern const pta processor_alias_table[];
extern unsigned int const pta_size;
extern unsigned int const num_arch_names;
#endif

#endif

extern enum processor_type ix86_tune;
extern enum processor_type ix86_arch;

/* Size of the RED_ZONE area.  */
#define RED_ZONE_SIZE 128
/* Reserved area of the red zone for temporaries.  */
#define RED_ZONE_RESERVE 8

extern unsigned int ix86_preferred_stack_boundary;
extern unsigned int ix86_incoming_stack_boundary;

/* Smallest class containing REGNO.  */
extern enum reg_class const regclass_map[FIRST_PSEUDO_REGISTER];

enum ix86_fpcmp_strategy {
  IX86_FPCMP_SAHF,
  IX86_FPCMP_COMI,
  IX86_FPCMP_ARITH
};

/* To properly truncate FP values into integers, we need to set i387 control
   word.  We can't emit proper mode switching code before reload, as spills
   generated by reload may truncate values incorrectly, but we still can avoid
   redundant computation of new control word by the mode switching pass.
   The fldcw instructions are still emitted redundantly, but this is probably
   not going to be noticeable problem, as most CPUs do have fast path for
   the sequence.

   The machinery is to emit simple truncation instructions and split them
   before reload to instructions having USEs of two memory locations that
   are filled by this code to old and new control word.

   Post-reload pass may be later used to eliminate the redundant fildcw if
   needed.  */

enum ix86_stack_slot
{
  SLOT_TEMP = 0,
  SLOT_CW_STORED,
  SLOT_CW_ROUNDEVEN,
  SLOT_CW_TRUNC,
  SLOT_CW_FLOOR,
  SLOT_CW_CEIL,
  SLOT_STV_TEMP,
  SLOT_FLOATxFDI_387,
  MAX_386_STACK_LOCALS
};

enum ix86_entity
{
  X86_DIRFLAG = 0,
  AVX_U128,
  I387_ROUNDEVEN,
  I387_TRUNC,
  I387_FLOOR,
  I387_CEIL,
  MAX_386_ENTITIES
};

enum x86_dirflag_state
{
  X86_DIRFLAG_RESET,
  X86_DIRFLAG_ANY
};

enum avx_u128_state
{
  AVX_U128_CLEAN,
  AVX_U128_DIRTY,
  AVX_U128_ANY
};

/* Define this macro if the port needs extra instructions inserted
   for mode switching in an optimizing compilation.  */

#define OPTIMIZE_MODE_SWITCHING(ENTITY) \
   ix86_optimize_mode_switching[(ENTITY)]

/* If you define `OPTIMIZE_MODE_SWITCHING', you have to define this as
   initializer for an array of integers.  Each initializer element N
   refers to an entity that needs mode switching, and specifies the
   number of different modes that might need to be set for this
   entity.  The position of the initializer in the initializer -
   starting counting at zero - determines the integer that is used to
   refer to the mode-switched entity in question.  */

#define NUM_MODES_FOR_MODE_SWITCHING			\
  { X86_DIRFLAG_ANY, AVX_U128_ANY,			\
    I387_CW_ANY, I387_CW_ANY, I387_CW_ANY, I387_CW_ANY  }


/* Avoid renaming of stack registers, as doing so in combination with
   scheduling just increases amount of live registers at time and in
   the turn amount of fxch instructions needed.

   ??? Maybe Pentium chips benefits from renaming, someone can try....

   Don't rename evex to non-evex sse registers.  */

#define HARD_REGNO_RENAME_OK(SRC, TARGET)				\
  (!STACK_REGNO_P (SRC)							\
   && EXT_REX_SSE_REGNO_P (SRC) == EXT_REX_SSE_REGNO_P (TARGET))


#define FASTCALL_PREFIX '@'

#ifndef USED_FOR_TARGET
/* Structure describing stack frame layout.
   Stack grows downward:

   [arguments]
					<- ARG_POINTER
   saved pc

   saved static chain			if ix86_static_chain_on_stack

   saved frame pointer			if frame_pointer_needed
					<- HARD_FRAME_POINTER
   [saved regs]
					<- reg_save_offset
   [padding0]
					<- stack_realign_offset
   [saved SSE regs]
	OR
   [stub-saved registers for ms x64 --> sysv clobbers
			<- Start of out-of-line, stub-saved/restored regs
			   (see libgcc/config/i386/(sav|res)ms64*.S)
     [XMM6-15]
     [RSI]
     [RDI]
     [?RBX]		only if RBX is clobbered
     [?RBP]		only if RBP and RBX are clobbered
     [?R12]		only if R12 and all previous regs are clobbered
     [?R13]		only if R13 and all previous regs are clobbered
     [?R14]		only if R14 and all previous regs are clobbered
     [?R15]		only if R15 and all previous regs are clobbered
			<- end of stub-saved/restored regs
     [padding1]
   ]
					<- sse_reg_save_offset
   [padding2]
		       |		<- FRAME_POINTER
   [va_arg registers]  |
		       |
   [frame]	       |
		       |
   [padding2]	       | = to_allocate
					<- STACK_POINTER
  */
struct GTY(()) ix86_frame
{
  int nsseregs;
  int nregs;
  int va_arg_size;
  int red_zone_size;
  int outgoing_arguments_size;

  /* The offsets relative to ARG_POINTER.  */
  HOST_WIDE_INT frame_pointer_offset;
  HOST_WIDE_INT hard_frame_pointer_offset;
  HOST_WIDE_INT stack_pointer_offset;
  HOST_WIDE_INT hfp_save_offset;
  HOST_WIDE_INT reg_save_offset;
  HOST_WIDE_INT stack_realign_allocate;
  HOST_WIDE_INT stack_realign_offset;
  HOST_WIDE_INT sse_reg_save_offset;

  /* When save_regs_using_mov is set, emit prologue using
     move instead of push instructions.  */
  bool save_regs_using_mov;

  /* Assume without checking that:
       EXPENSIVE_P = expensive_function_p (EXPENSIVE_COUNT).  */
  bool expensive_p;
  int expensive_count;
};

/* Machine specific frame tracking during prologue/epilogue generation.  All
   values are positive, but since the x86 stack grows downward, are subtratced
   from the CFA to produce a valid address.  */

struct GTY(()) machine_frame_state
{
  /* This pair tracks the currently active CFA as reg+offset.  When reg
     is drap_reg, we don't bother trying to record here the real CFA when
     it might really be a DW_CFA_def_cfa_expression.  */
  rtx cfa_reg;
  HOST_WIDE_INT cfa_offset;

  /* The current offset (canonically from the CFA) of ESP and EBP.
     When stack frame re-alignment is active, these may not be relative
     to the CFA.  However, in all cases they are relative to the offsets
     of the saved registers stored in ix86_frame.  */
  HOST_WIDE_INT sp_offset;
  HOST_WIDE_INT fp_offset;

  /* The size of the red-zone that may be assumed for the purposes of
     eliding register restore notes in the epilogue.  This may be zero
     if no red-zone is in effect, or may be reduced from the real
     red-zone value by a maximum runtime stack re-alignment value.  */
  int red_zone_offset;

  /* Indicate whether each of ESP, EBP or DRAP currently holds a valid
     value within the frame.  If false then the offset above should be
     ignored.  Note that DRAP, if valid, *always* points to the CFA and
     thus has an offset of zero.  */
  BOOL_BITFIELD sp_valid : 1;
  BOOL_BITFIELD fp_valid : 1;
  BOOL_BITFIELD drap_valid : 1;

  /* Indicate whether the local stack frame has been re-aligned.  When
     set, the SP/FP offsets above are relative to the aligned frame
     and not the CFA.  */
  BOOL_BITFIELD realigned : 1;

  /* Indicates whether the stack pointer has been re-aligned.  When set,
     SP/FP continue to be relative to the CFA, but the stack pointer
     should only be used for offsets > sp_realigned_offset, while
     the frame pointer should be used for offsets <= sp_realigned_fp_last.
     The flags realigned and sp_realigned are mutually exclusive.  */
  BOOL_BITFIELD sp_realigned : 1;

  /* If sp_realigned is set, this is the last valid offset from the CFA
     that can be used for access with the frame pointer.  */
  HOST_WIDE_INT sp_realigned_fp_last;

  /* If sp_realigned is set, this is the offset from the CFA that the stack
     pointer was realigned, and may or may not be equal to sp_realigned_fp_last.
     Access via the stack pointer is only valid for offsets that are greater than
     this value.  */
  HOST_WIDE_INT sp_realigned_offset;
};

/* Private to winnt.cc.  */
struct seh_frame_state;

enum function_type
{
  TYPE_UNKNOWN = 0,
  TYPE_NORMAL,
  /* The current function is an interrupt service routine with a
     pointer argument as specified by the "interrupt" attribute.  */
  TYPE_INTERRUPT,
  /* The current function is an interrupt service routine with a
     pointer argument and an integer argument as specified by the
     "interrupt" attribute.  */
  TYPE_EXCEPTION
};

enum queued_insn_type
{
  TYPE_NONE = 0,
  TYPE_ENDBR,
  TYPE_PATCHABLE_AREA
};

struct GTY(()) machine_function {
  struct stack_local_entry *stack_locals;
  int varargs_gpr_size;
  int varargs_fpr_size;
  int optimize_mode_switching[MAX_386_ENTITIES];

  /* Cached initial frame layout for the current function.  */
  struct ix86_frame frame;

  /* For -fsplit-stack support: A stack local which holds a pointer to
     the stack arguments for a function with a variable number of
     arguments.  This is set at the start of the function and is used
     to initialize the overflow_arg_area field of the va_list
     structure.  */
  rtx split_stack_varargs_pointer;

  /* This value is used for amd64 targets and specifies the current abi
     to be used. MS_ABI means ms abi. Otherwise SYSV_ABI means sysv abi.  */
  ENUM_BITFIELD(calling_abi) call_abi : 8;

  /* Nonzero if the function accesses a previous frame.  */
  BOOL_BITFIELD accesses_prev_frame : 1;

  /* Set by ix86_compute_frame_layout and used by prologue/epilogue
     expander to determine the style used.  */
  BOOL_BITFIELD use_fast_prologue_epilogue : 1;

  /* Nonzero if the current function calls pc thunk and
     must not use the red zone.  */
  BOOL_BITFIELD pc_thunk_call_expanded : 1;

  /* If true, the current function needs the default PIC register, not
     an alternate register (on x86) and must not use the red zone (on
     x86_64), even if it's a leaf function.  We don't want the
     function to be regarded as non-leaf because TLS calls need not
     affect register allocation.  This flag is set when a TLS call
     instruction is expanded within a function, and never reset, even
     if all such instructions are optimized away.  Use the
     ix86_current_function_calls_tls_descriptor macro for a better
     approximation.  */
  BOOL_BITFIELD tls_descriptor_call_expanded_p : 1;

  /* If true, the current function has a STATIC_CHAIN is placed on the
     stack below the return address.  */
  BOOL_BITFIELD static_chain_on_stack : 1;

  /* If true, it is safe to not save/restore DRAP register.  */
  BOOL_BITFIELD no_drap_save_restore : 1;

  /* Function type.  */
  ENUM_BITFIELD(function_type) func_type : 2;

  /* How to generate indirec branch.  */
  ENUM_BITFIELD(indirect_branch) indirect_branch_type : 3;

  /* If true, the current function has local indirect jumps, like
     "indirect_jump" or "tablejump".  */
  BOOL_BITFIELD has_local_indirect_jump : 1;

  /* How to generate function return.  */
  ENUM_BITFIELD(indirect_branch) function_return_type : 3;

  /* If true, the current function is a function specified with
     the "interrupt" or "no_caller_saved_registers" attribute.  */
  BOOL_BITFIELD no_caller_saved_registers : 1;

  /* If true, there is register available for argument passing.  This
     is used only in ix86_function_ok_for_sibcall by 32-bit to determine
     if there is scratch register available for indirect sibcall.  In
     64-bit, rax, r10 and r11 are scratch registers which aren't used to
     pass arguments and can be used for indirect sibcall.  */
  BOOL_BITFIELD arg_reg_available : 1;

  /* If true, we're out-of-lining reg save/restore for regs clobbered
     by 64-bit ms_abi functions calling a sysv_abi function.  */
  BOOL_BITFIELD call_ms2sysv : 1;

  /* If true, the incoming 16-byte aligned stack has an offset (of 8) and
     needs padding prior to out-of-line stub save/restore area.  */
  BOOL_BITFIELD call_ms2sysv_pad_in : 1;

  /* This is the number of extra registers saved by stub (valid range is
     0-6). Each additional register is only saved/restored by the stubs
     if all successive ones are. (Will always be zero when using a hard
     frame pointer.) */
  unsigned int call_ms2sysv_extra_regs:3;

  /* Nonzero if the function places outgoing arguments on stack.  */
  BOOL_BITFIELD outgoing_args_on_stack : 1;

  /* If true, ENDBR or patchable area is queued at function entrance.  */
  ENUM_BITFIELD(queued_insn_type) insn_queued_at_entrance : 2;

  /* If true, the function label has been emitted.  */
  BOOL_BITFIELD function_label_emitted : 1;

  /* True if the function needs a stack frame.  */
  BOOL_BITFIELD stack_frame_required : 1;

  /* True if we should act silently, rather than raise an error for
     invalid calls.  */
  BOOL_BITFIELD silent_p : 1;

  /* True if red zone is used.  */
  BOOL_BITFIELD red_zone_used : 1;

  /* The largest alignment, in bytes, of stack slot actually used.  */
  unsigned int max_used_stack_alignment;

  /* During prologue/epilogue generation, the current frame state.
     Otherwise, the frame state at the end of the prologue.  */
  struct machine_frame_state fs;

  /* During SEH output, this is non-null.  */
  struct seh_frame_state * GTY((skip(""))) seh;
};

extern GTY(()) tree sysv_va_list_type_node;
extern GTY(()) tree ms_va_list_type_node;
#endif

#define ix86_stack_locals (cfun->machine->stack_locals)
#define ix86_varargs_gpr_size (cfun->machine->varargs_gpr_size)
#define ix86_varargs_fpr_size (cfun->machine->varargs_fpr_size)
#define ix86_optimize_mode_switching (cfun->machine->optimize_mode_switching)
#define ix86_pc_thunk_call_expanded (cfun->machine->pc_thunk_call_expanded)
#define ix86_tls_descriptor_calls_expanded_in_cfun \
  (cfun->machine->tls_descriptor_call_expanded_p)
/* Since tls_descriptor_call_expanded is not cleared, even if all TLS
   calls are optimized away, we try to detect cases in which it was
   optimized away.  Since such instructions (use (reg REG_SP)), we can
   verify whether there's any such instruction live by testing that
   REG_SP is live.  */
#define ix86_current_function_calls_tls_descriptor \
  (ix86_tls_descriptor_calls_expanded_in_cfun && df_regs_ever_live_p (SP_REG))
#define ix86_static_chain_on_stack (cfun->machine->static_chain_on_stack)
#define ix86_red_zone_used (cfun->machine->red_zone_used)

/* Control behavior of x86_file_start.  */
#define X86_FILE_START_VERSION_DIRECTIVE false
#define X86_FILE_START_FLTUSED false

/* Flag to mark data that is in the large address area.  */
#define SYMBOL_FLAG_FAR_ADDR		(SYMBOL_FLAG_MACH_DEP << 0)
#define SYMBOL_REF_FAR_ADDR_P(X)	\
	((SYMBOL_REF_FLAGS (X) & SYMBOL_FLAG_FAR_ADDR) != 0)

/* Flags to mark dllimport/dllexport.  Used by PE ports, but handy to
   have defined always, to avoid ifdefing.  */
#define SYMBOL_FLAG_DLLIMPORT		(SYMBOL_FLAG_MACH_DEP << 1)
#define SYMBOL_REF_DLLIMPORT_P(X) \
	((SYMBOL_REF_FLAGS (X) & SYMBOL_FLAG_DLLIMPORT) != 0)

#define SYMBOL_FLAG_DLLEXPORT		(SYMBOL_FLAG_MACH_DEP << 2)
#define SYMBOL_REF_DLLEXPORT_P(X) \
	((SYMBOL_REF_FLAGS (X) & SYMBOL_FLAG_DLLEXPORT) != 0)

#define SYMBOL_FLAG_STUBVAR	(SYMBOL_FLAG_MACH_DEP << 4)
#define SYMBOL_REF_STUBVAR_P(X) \
	((SYMBOL_REF_FLAGS (X) & SYMBOL_FLAG_STUBVAR) != 0)

extern void debug_ready_dispatch (void);
extern void debug_dispatch_window (int);

/* The value at zero is only defined for the BMI instructions
   LZCNT and TZCNT, not the BSR/BSF insns in the original isa.  */
#define CTZ_DEFINED_VALUE_AT_ZERO(MODE, VALUE) \
	((VALUE) = GET_MODE_BITSIZE (MODE), TARGET_BMI ? 2 : 0)
#define CLZ_DEFINED_VALUE_AT_ZERO(MODE, VALUE) \
	((VALUE) = GET_MODE_BITSIZE (MODE), TARGET_LZCNT ? 2 : 0)


/* Flags returned by ix86_get_callcvt ().  */
#define IX86_CALLCVT_CDECL	0x1
#define IX86_CALLCVT_STDCALL	0x2
#define IX86_CALLCVT_FASTCALL	0x4
#define IX86_CALLCVT_THISCALL	0x8
#define IX86_CALLCVT_REGPARM	0x10
#define IX86_CALLCVT_SSEREGPARM	0x20

#define IX86_BASE_CALLCVT(FLAGS) \
	((FLAGS) & (IX86_CALLCVT_CDECL | IX86_CALLCVT_STDCALL \
		    | IX86_CALLCVT_FASTCALL | IX86_CALLCVT_THISCALL))

#define RECIP_MASK_NONE		0x00
#define RECIP_MASK_DIV		0x01
#define RECIP_MASK_SQRT		0x02
#define RECIP_MASK_VEC_DIV	0x04
#define RECIP_MASK_VEC_SQRT	0x08
#define RECIP_MASK_ALL	(RECIP_MASK_DIV | RECIP_MASK_SQRT \
			 | RECIP_MASK_VEC_DIV | RECIP_MASK_VEC_SQRT)
#define RECIP_MASK_DEFAULT (RECIP_MASK_VEC_DIV | RECIP_MASK_VEC_SQRT)

#define TARGET_RECIP_DIV	((recip_mask & RECIP_MASK_DIV) != 0)
#define TARGET_RECIP_SQRT	((recip_mask & RECIP_MASK_SQRT) != 0)
#define TARGET_RECIP_VEC_DIV	((recip_mask & RECIP_MASK_VEC_DIV) != 0)
#define TARGET_RECIP_VEC_SQRT	((recip_mask & RECIP_MASK_VEC_SQRT) != 0)

/* Use 128-bit AVX instructions in the auto-vectorizer.  */
#define TARGET_PREFER_AVX128	(prefer_vector_width_type == PVW_AVX128)
/* Use 256-bit AVX instructions in the auto-vectorizer.  */
#define TARGET_PREFER_AVX256	(TARGET_PREFER_AVX128 \
				 || prefer_vector_width_type == PVW_AVX256)

#define TARGET_INDIRECT_BRANCH_REGISTER \
  (ix86_indirect_branch_register \
   || cfun->machine->indirect_branch_type != indirect_branch_keep)

#define IX86_HLE_ACQUIRE (1 << 16)
#define IX86_HLE_RELEASE (1 << 17)

/* For switching between functions with different target attributes.  */
#define SWITCHABLE_TARGET 1

#define TARGET_SUPPORTS_WIDE_INT 1

#if !defined(GENERATOR_FILE) && !defined(IN_LIBGCC2)
extern enum attr_cpu ix86_schedule;

#define NUM_X86_64_MS_CLOBBERED_REGS 12
#endif

/* __builtin_eh_return can't handle stack realignment, so disable MMX/SSE
   in 32-bit libgcc functions that call it.  */
#ifndef __x86_64__
#define LIBGCC2_UNWIND_ATTRIBUTE __attribute__((target ("no-mmx,no-sse")))
#endif

/*
Local variables:
version-control: t
End:
*/

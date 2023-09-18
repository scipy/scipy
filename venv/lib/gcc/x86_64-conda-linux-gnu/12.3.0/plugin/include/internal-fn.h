/* Internal functions.
   Copyright (C) 2011-2022 Free Software Foundation, Inc.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3, or (at your option) any later
version.

GCC is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

#ifndef GCC_INTERNAL_FN_H
#define GCC_INTERNAL_FN_H

/* INTEGER_CST values for IFN_UNIQUE function arg-0.

   UNSPEC: Undifferentiated UNIQUE.

   FORK and JOIN mark the points at which OpenACC partitioned
   execution is entered or exited.
      DEP_VAR = UNIQUE ({FORK,JOIN}, DEP_VAR, AXIS)

   HEAD_MARK and TAIL_MARK are used to demark the sequence entering
   or leaving partitioned execution.
      DEP_VAR = UNIQUE ({HEAD,TAIL}_MARK, REMAINING_MARKS, ...PRIMARY_FLAGS)

   The PRIMARY_FLAGS only occur on the first HEAD_MARK of a sequence.

   PRIVATE captures variables to be made private at the surrounding parallelism
   level.  */
#define IFN_UNIQUE_CODES				  \
  DEF(UNSPEC),	\
    DEF(OACC_FORK), DEF(OACC_JOIN),		\
    DEF(OACC_HEAD_MARK), DEF(OACC_TAIL_MARK),	\
    DEF(OACC_PRIVATE)

enum ifn_unique_kind {
#define DEF(X) IFN_UNIQUE_##X
  IFN_UNIQUE_CODES
#undef DEF
};

/* INTEGER_CST values for IFN_GOACC_LOOP arg-0.  Allows the precise
   stepping of the compute geometry over the loop iterations to be
   deferred until it is known which compiler is generating the code.
   The action is encoded in a constant first argument.

     CHUNK_MAX = LOOP (CODE_CHUNKS, DIR, RANGE, STEP, CHUNK_SIZE, MASK)
     STEP = LOOP (CODE_STEP, DIR, RANGE, STEP, CHUNK_SIZE, MASK)
     OFFSET = LOOP (CODE_OFFSET, DIR, RANGE, STEP, CHUNK_SIZE, MASK, CHUNK_NO)
     BOUND = LOOP (CODE_BOUND, DIR, RANGE, STEP, CHUNK_SIZE, MASK, OFFSET)

     DIR - +1 for up loop, -1 for down loop
     RANGE - Range of loop (END - BASE)
     STEP - iteration step size
     CHUNKING - size of chunking, (constant zero for no chunking)
     CHUNK_NO - chunk number
     MASK - partitioning mask.  */

#define IFN_GOACC_LOOP_CODES \
  DEF(CHUNKS), DEF(STEP), DEF(OFFSET), DEF(BOUND)
enum ifn_goacc_loop_kind {
#define DEF(X) IFN_GOACC_LOOP_##X
  IFN_GOACC_LOOP_CODES
#undef DEF
};

/* The GOACC_REDUCTION function defines a generic interface to support
   gang, worker and vector reductions.  All calls are of the following
   form:

     V = REDUCTION (CODE, REF_TO_RES, LOCAL_VAR, LEVEL, OP, OFFSET)

   REF_TO_RES - is a reference to the original reduction varl, may be NULL
   LOCAL_VAR is the intermediate reduction variable
   LEVEL corresponds to the GOMP_DIM of the reduction
   OP is the tree code of the reduction operation
   OFFSET may be used as an offset into a reduction array for the
          reductions occuring at this level.
   In general the return value is LOCAL_VAR, which creates a data
   dependency between calls operating on the same reduction.  */

#define IFN_GOACC_REDUCTION_CODES \
  DEF(SETUP), DEF(INIT), DEF(FINI), DEF(TEARDOWN)
enum ifn_goacc_reduction_kind {
#define DEF(X) IFN_GOACC_REDUCTION_##X
  IFN_GOACC_REDUCTION_CODES
#undef DEF
};

/* Initialize internal function tables.  */

extern void init_internal_fns ();

/* Return the name of internal function FN.  The name is only meaningful
   for dumps; it has no linkage.  */

extern const char *const internal_fn_name_array[];

static inline const char *
internal_fn_name (enum internal_fn fn)
{
  return internal_fn_name_array[(int) fn];
}

extern internal_fn lookup_internal_fn (const char *);

/* Return the ECF_* flags for function FN.  */

extern const int internal_fn_flags_array[];

static inline int
internal_fn_flags (enum internal_fn fn)
{
  return internal_fn_flags_array[(int) fn];
}

/* Return fnspec for function FN.  */

extern GTY(()) const_tree internal_fn_fnspec_array[IFN_LAST + 1];

static inline const_tree
internal_fn_fnspec (enum internal_fn fn)
{
  return internal_fn_fnspec_array[(int) fn];
}

/* Describes an internal function that maps directly to an optab.  */
struct direct_internal_fn_info
{
  /* optabs can be parameterized by one or two modes.  These fields describe
     how to select those modes from the types of the return value and
     arguments.  A value of -1 says that the mode is determined by the
     return type while a value N >= 0 says that the mode is determined by
     the type of argument N.  A value of -2 says that this internal
     function isn't directly mapped to an optab.  */
  signed int type0 : 8;
  signed int type1 : 8;
  /* True if the function is pointwise, so that it can be vectorized by
     converting the return type and all argument types to vectors of the
     same number of elements.  E.g. we can vectorize an IFN_SQRT on
     floats as an IFN_SQRT on vectors of N floats.

     This only needs 1 bit, but occupies the full 16 to ensure a nice
     layout.  */
  unsigned int vectorizable : 16;
};

extern const direct_internal_fn_info direct_internal_fn_array[IFN_LAST + 1];

/* Return true if FN is mapped directly to an optab.  */

inline bool
direct_internal_fn_p (internal_fn fn)
{
  return direct_internal_fn_array[fn].type0 >= -1;
}

/* Return true if FN is a direct internal function that can be vectorized by
   converting the return type and all argument types to vectors of the same
   number of elements.  E.g. we can vectorize an IFN_SQRT on floats as an
   IFN_SQRT on vectors of N floats.  */

inline bool
vectorizable_internal_fn_p (internal_fn fn)
{
  return direct_internal_fn_array[fn].vectorizable;
}

/* Return optab information about internal function FN.  Only meaningful
   if direct_internal_fn_p (FN).  */

inline const direct_internal_fn_info &
direct_internal_fn (internal_fn fn)
{
  gcc_checking_assert (direct_internal_fn_p (fn));
  return direct_internal_fn_array[fn];
}

extern tree_pair direct_internal_fn_types (internal_fn, tree, tree *);
extern tree_pair direct_internal_fn_types (internal_fn, gcall *);
extern bool direct_internal_fn_supported_p (internal_fn, tree_pair,
					    optimization_type);
extern bool direct_internal_fn_supported_p (internal_fn, tree,
					    optimization_type);
extern bool direct_internal_fn_supported_p (gcall *, optimization_type);

/* Return true if FN is supported for types TYPE0 and TYPE1 when the
   optimization type is OPT_TYPE.  The types are those associated with
   the "type0" and "type1" fields of FN's direct_internal_fn_info
   structure.  */

inline bool
direct_internal_fn_supported_p (internal_fn fn, tree type0, tree type1,
				optimization_type opt_type)
{
  return direct_internal_fn_supported_p (fn, tree_pair (type0, type1),
					 opt_type);
}

extern bool commutative_binary_fn_p (internal_fn);
extern bool commutative_ternary_fn_p (internal_fn);
extern int first_commutative_argument (internal_fn);
extern bool associative_binary_fn_p (internal_fn);

extern bool set_edom_supported_p (void);

extern internal_fn get_conditional_internal_fn (tree_code);
extern internal_fn get_conditional_internal_fn (internal_fn);
extern tree_code conditional_internal_fn_code (internal_fn);
extern internal_fn get_unconditional_internal_fn (internal_fn);
extern bool can_interpret_as_conditional_op_p (gimple *, tree *,
					       tree_code *, tree (&)[3],
					       tree *);

extern bool internal_load_fn_p (internal_fn);
extern bool internal_store_fn_p (internal_fn);
extern bool internal_gather_scatter_fn_p (internal_fn);
extern int internal_fn_mask_index (internal_fn);
extern int internal_fn_stored_value_index (internal_fn);
extern bool internal_gather_scatter_fn_supported_p (internal_fn, tree,
						    tree, tree, int);
extern bool internal_check_ptrs_fn_supported_p (internal_fn, tree,
						poly_uint64, unsigned int);
#define VECT_PARTIAL_BIAS_UNSUPPORTED 127

extern signed char internal_len_load_store_bias (internal_fn ifn,
						 machine_mode);

extern void expand_addsub_overflow (location_t, tree_code, tree, tree, tree,
				    bool, bool, bool, bool, tree *);
extern void expand_internal_call (gcall *);
extern void expand_internal_call (internal_fn, gcall *);
extern void expand_PHI (internal_fn, gcall *);
extern void expand_SHUFFLEVECTOR (internal_fn, gcall *);
extern void expand_SPACESHIP (internal_fn, gcall *);

extern bool vectorized_internal_fn_supported_p (internal_fn, tree);

enum {
  ATOMIC_OP_FETCH_CMP_0_EQ = 0,
  ATOMIC_OP_FETCH_CMP_0_NE = 1,
  ATOMIC_OP_FETCH_CMP_0_LT = 2,
  ATOMIC_OP_FETCH_CMP_0_LE = 3,
  ATOMIC_OP_FETCH_CMP_0_GT = 4,
  ATOMIC_OP_FETCH_CMP_0_GE = 5
};

#endif

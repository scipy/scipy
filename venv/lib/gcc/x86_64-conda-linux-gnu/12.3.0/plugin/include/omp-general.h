/* General types and functions that are uselful for processing of OpenMP,
   OpenACC and similar directivers at various stages of compilation.

   Copyright (C) 2005-2022 Free Software Foundation, Inc.

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

#ifndef GCC_OMP_GENERAL_H
#define GCC_OMP_GENERAL_H

#include "gomp-constants.h"

/*  Flags for an OpenACC loop.  */

enum oacc_loop_flags {
  OLF_SEQ	= 1u << 0,  /* Explicitly sequential  */
  OLF_AUTO	= 1u << 1,	/* Compiler chooses axes.  */
  OLF_INDEPENDENT = 1u << 2,	/* Iterations are known independent.  */
  OLF_GANG_STATIC = 1u << 3,	/* Gang partitioning is static (has op). */
  OLF_TILE	= 1u << 4,	/* Tiled loop. */
  OLF_REDUCTION = 1u << 5,	/* Reduction loop.  */
  
  /* Explicitly specified loop axes.  */
  OLF_DIM_BASE = 6,
  OLF_DIM_GANG   = 1u << (OLF_DIM_BASE + GOMP_DIM_GANG),
  OLF_DIM_WORKER = 1u << (OLF_DIM_BASE + GOMP_DIM_WORKER),
  OLF_DIM_VECTOR = 1u << (OLF_DIM_BASE + GOMP_DIM_VECTOR),

  OLF_MAX = OLF_DIM_BASE + GOMP_DIM_MAX
};

/* A structure holding the elements of:
   for (V = N1; V cond N2; V += STEP) [...]
   or for non-rectangular loops:
   for (V = M1 * W + N1; V cond M2 * W + N2; V += STEP;
   where W is V of the OUTER-th loop (e.g. for OUTER 1 it is the
   the index of the immediately surrounding loop).
   NON_RECT_REFERENCED is true for loops referenced by loops
   with non-NULL M1 or M2.  */

struct omp_for_data_loop
{
  tree v, n1, n2, step, m1, m2;
  enum tree_code cond_code;
  int outer;
  bool non_rect_referenced;
};

/* A structure describing the main elements of a parallel loop.  */

struct omp_for_data
{
  struct omp_for_data_loop loop;
  tree chunk_size;
  gomp_for *for_stmt;
  tree pre, iter_type;
  tree tiling;  /* Tiling values (if non null).  */
  int collapse;  /* Collapsed loops, 1 for a non-collapsed loop.  */
  int ordered;
  int first_nonrect, last_nonrect;
  bool have_nowait, have_ordered, simd_schedule, have_reductemp;
  bool have_pointer_condtemp, have_scantemp, have_nonctrl_scantemp;
  bool non_rect;
  int lastprivate_conditional;
  unsigned char sched_modifiers;
  enum omp_clause_schedule_kind sched_kind;
  struct omp_for_data_loop *loops;
  /* The following are relevant only for non-rectangular loops
     where only a single loop depends on an outer loop iterator.  */
  tree first_inner_iterations; /* Number of iterations of the inner
				  loop with the first outer iterator
				  (or adjn1, if that is non-NULL).  */
  tree factor; /* (m2 - m1) * outer_step / inner_step.  */
  /* Adjusted n1 of the outer loop in such loop nests (if needed).  */
  tree adjn1;
};

#define OACC_FN_ATTRIB "oacc function"

extern tree omp_find_clause (tree clauses, enum omp_clause_code kind);
extern bool omp_is_allocatable_or_ptr (tree decl);
extern tree omp_check_optional_argument (tree decl, bool for_present_check);
extern bool omp_privatize_by_reference (tree decl);
extern void omp_adjust_for_condition (location_t loc, enum tree_code *cond_code,
				      tree *n2, tree v, tree step);
extern tree omp_get_for_step_from_incr (location_t loc, tree incr);
extern void omp_extract_for_data (gomp_for *for_stmt, struct omp_for_data *fd,
				  struct omp_for_data_loop *loops);
extern gimple *omp_build_barrier (tree lhs);
extern tree find_combined_omp_for (tree *, int *, void *);
extern poly_uint64 omp_max_vf (void);
extern int omp_max_simt_vf (void);
extern int omp_constructor_traits_to_codes (tree, enum tree_code *);
extern tree omp_check_context_selector (location_t loc, tree ctx);
extern void omp_mark_declare_variant (location_t loc, tree variant,
				      tree construct);
extern int omp_context_selector_matches (tree);
extern int omp_context_selector_set_compare (const char *, tree, tree);
extern tree omp_get_context_selector (tree, const char *, const char *);
extern tree omp_resolve_declare_variant (tree);
extern tree oacc_launch_pack (unsigned code, tree device, unsigned op);
extern tree oacc_replace_fn_attrib_attr (tree attribs, tree dims);
extern void oacc_replace_fn_attrib (tree fn, tree dims);
extern void oacc_set_fn_attrib (tree fn, tree clauses, vec<tree> *args);
extern int oacc_verify_routine_clauses (tree, tree *, location_t,
					const char *);
extern tree oacc_build_routine_dims (tree clauses);
extern tree oacc_get_fn_attrib (tree fn);
extern bool offloading_function_p (tree fn);
extern int oacc_get_fn_dim_size (tree fn, int axis);
extern int oacc_get_ifn_dim_arg (const gimple *stmt);

enum omp_requires {
  OMP_REQUIRES_ATOMIC_DEFAULT_MEM_ORDER = 0xf,
  OMP_REQUIRES_UNIFIED_ADDRESS = 0x10,
  OMP_REQUIRES_UNIFIED_SHARED_MEMORY = 0x20,
  OMP_REQUIRES_DYNAMIC_ALLOCATORS = 0x40,
  OMP_REQUIRES_REVERSE_OFFLOAD = 0x80,
  OMP_REQUIRES_ATOMIC_DEFAULT_MEM_ORDER_USED = 0x100,
  OMP_REQUIRES_TARGET_USED = 0x200
};

extern GTY(()) enum omp_requires omp_requires_mask;

static inline dump_flags_t
get_openacc_privatization_dump_flags ()
{
  dump_flags_t l_dump_flags = MSG_NOTE;

  /* For '--param=openacc-privatization=quiet', diagnostics only go to dump
     files.  */
  if (param_openacc_privatization == OPENACC_PRIVATIZATION_QUIET)
    l_dump_flags |= MSG_PRIORITY_INTERNALS;

  return l_dump_flags;
}

extern tree omp_build_component_ref (tree obj, tree field);

#endif /* GCC_OMP_GENERAL_H */

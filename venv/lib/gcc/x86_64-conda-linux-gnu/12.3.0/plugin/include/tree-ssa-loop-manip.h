/* Header file for High-level loop manipulation functions.
   Copyright (C) 2013-2022 Free Software Foundation, Inc.

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

#ifndef GCC_TREE_SSA_LOOP_MANIP_H
#define GCC_TREE_SSA_LOOP_MANIP_H

typedef void (*transform_callback)(class loop *, void *);

extern void create_iv (tree, tree, tree, class loop *, gimple_stmt_iterator *,
		       bool, tree *, tree *);
extern void rewrite_into_loop_closed_ssa_1 (bitmap, unsigned, int,
					    class loop *);
extern void rewrite_into_loop_closed_ssa (bitmap, unsigned);
extern void rewrite_virtuals_into_loop_closed_ssa (class loop *);
extern void verify_loop_closed_ssa (bool, class loop * = NULL);

static inline void
checking_verify_loop_closed_ssa (bool verify_ssa_p, class loop *loop = NULL)
{
  if (flag_checking)
    verify_loop_closed_ssa (verify_ssa_p, loop);
}

extern basic_block split_loop_exit_edge (edge, bool = false);
extern basic_block ip_end_pos (class loop *);
extern basic_block ip_normal_pos (class loop *);
extern void standard_iv_increment_position (class loop *,
					    gimple_stmt_iterator *, bool *);
extern bool
gimple_duplicate_loop_body_to_header_edge (class loop *, edge, unsigned int,
					   sbitmap, edge, vec<edge> *, int);
extern bool can_unroll_loop_p (class loop *loop, unsigned factor,
			       class tree_niter_desc *niter);
extern gcov_type niter_for_unrolled_loop (class loop *, unsigned);
extern void tree_transform_and_unroll_loop (class loop *, unsigned,
					    tree_niter_desc *,
					    transform_callback, void *);
extern void tree_unroll_loop (class loop *, unsigned, tree_niter_desc *);
extern tree canonicalize_loop_ivs (class loop *, tree *, bool);
extern unsigned int loop_invariant_motion_in_fun (function *, bool);


#endif /* GCC_TREE_SSA_LOOP_MANIP_H */

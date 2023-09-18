/* Loop manipulation header.
   Copyright (C) 2014-2022 Free Software Foundation, Inc.

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

#ifndef GCC_CFGLOOPMANIP_H
#define GCC_CFGLOOPMANIP_H

enum
{
  CP_SIMPLE_PREHEADERS = 1,
  CP_FALLTHRU_PREHEADERS = 2
};

#define DLTHE_FLAG_UPDATE_FREQ	1	/* Update frequencies in
					   duplicate_loop_to_header_edge.  */
#define DLTHE_RECORD_COPY_NUMBER 2	/* Record copy number in the aux
					   field of newly create BB.  */
#define DLTHE_FLAG_COMPLETTE_PEEL 4	/* Update frequencies expecting
					   a complete peeling.  */
extern edge mfb_kj_edge;

extern bool remove_path (edge, bool * = NULL, bitmap = NULL);
extern void place_new_loop (struct function *, class loop *);
extern void add_loop (class loop *, class loop *);
extern void scale_loop_frequencies (class loop *, profile_probability);
extern void scale_loop_profile (class loop *, profile_probability, gcov_type);
extern edge create_empty_if_region_on_edge (edge, tree);
extern class loop *create_empty_loop_on_edge (edge, tree, tree, tree, tree,
					       tree *, tree *, class loop *);
extern void unloop (class loop *, bool *, bitmap);
extern void copy_loop_info (class loop *loop, class loop *target);
extern class loop * duplicate_loop (class loop *, class loop *,
				     class loop * = NULL);
extern void duplicate_subloops (class loop *, class loop *);
extern bool can_duplicate_loop_p (const class loop *loop);
extern bool
duplicate_loop_body_to_header_edge (class loop *, edge, unsigned, sbitmap, edge,
				    vec<edge> *, int);
extern bool mfb_keep_just (edge);
basic_block create_preheader (class loop *, int);
extern void create_preheaders (int);
extern void force_single_succ_latches (void);
class loop * loop_version (class loop *, void *,
			    basic_block *,
			    profile_probability, profile_probability,
			    profile_probability, profile_probability, bool);

#endif /* GCC_CFGLOOPMANIP_H */

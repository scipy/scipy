/* Control flow graph manipulation code header file.
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

#ifndef GCC_CFG_H
#define GCC_CFG_H

#include "dominance.h"

/* What sort of profiling information we have.  */
enum profile_status_d
{
  PROFILE_ABSENT,
  PROFILE_GUESSED,
  PROFILE_READ,
  PROFILE_LAST	/* Last value, used by profile streaming.  */
};

/* A structure to group all the per-function control flow graph data.
   The x_* prefixing is necessary because otherwise references to the
   fields of this struct are interpreted as the defines for backward
   source compatibility following the definition of this struct.  */
struct GTY(()) control_flow_graph {
  /* Block pointers for the exit and entry of a function.
     These are always the head and tail of the basic block list.  */
  basic_block x_entry_block_ptr;
  basic_block x_exit_block_ptr;

  /* Index by basic block number, get basic block struct info.  */
  vec<basic_block, va_gc> *x_basic_block_info;

  /* Number of basic blocks in this flow graph.  */
  int x_n_basic_blocks;

  /* Number of edges in this flow graph.  */
  int x_n_edges;

  /* The first free basic block number.  */
  int x_last_basic_block;

  /* UIDs for LABEL_DECLs.  */
  int last_label_uid;

  /* Mapping of labels to their associated blocks.  At present
     only used for the gimple CFG.  */
  vec<basic_block, va_gc> *x_label_to_block_map;

  enum profile_status_d x_profile_status;

  /* Whether the dominators and the postdominators are available.  */
  enum dom_state x_dom_computed[2];

  /* Number of basic blocks in the dominance tree.  */
  unsigned x_n_bbs_in_dom_tree[2];

  /* Maximal number of entities in the single jumptable.  Used to estimate
     final flowgraph size.  */
  int max_jumptable_ents;

  /* Maximal count of BB in function.  */
  profile_count count_max;

  /* Dynamically allocated edge/bb flags.  */
  int edge_flags_allocated;
  int bb_flags_allocated;
};


extern void init_flow (function *);
extern void free_cfg (function *);
extern basic_block alloc_block (void);
extern void link_block (basic_block, basic_block);
extern void unlink_block (basic_block);
extern void compact_blocks (void);
extern void expunge_block (basic_block);
extern edge unchecked_make_edge (basic_block, basic_block, int);
extern edge cached_make_edge (sbitmap, basic_block, basic_block, int);
extern edge make_edge (basic_block, basic_block, int);
extern edge make_single_succ_edge (basic_block, basic_block, int);
extern void remove_edge_raw (edge);
extern void redirect_edge_succ (edge, basic_block);
extern void redirect_edge_pred (edge, basic_block);
extern void clear_bb_flags (void);
extern void dump_edge_info (FILE *, edge, dump_flags_t, int);
extern void debug (edge_def &ref);
extern void debug (edge_def *ptr);
extern void alloc_aux_for_blocks (int);
extern void clear_aux_for_blocks (void);
extern void free_aux_for_blocks (void);
extern void alloc_aux_for_edge (edge, int);
extern void alloc_aux_for_edges (int);
extern void clear_aux_for_edges (void);
extern void free_aux_for_edges (void);
extern void debug_bb (basic_block);
extern basic_block debug_bb_n (int);
extern void debug_bb (basic_block, dump_flags_t);
extern basic_block debug_bb_n (int, dump_flags_t);
extern void dump_bb_info (FILE *, basic_block, int, dump_flags_t, bool, bool);
extern void brief_dump_cfg (FILE *, dump_flags_t);
extern void update_bb_profile_for_threading (basic_block, profile_count, edge);
extern void scale_bbs_frequencies_profile_count (basic_block *, int,
					     profile_count, profile_count);
extern void scale_bbs_frequencies (basic_block *, int, profile_probability);
extern void initialize_original_copy_tables (void);
extern void reset_original_copy_tables (void);
extern void free_original_copy_tables (void);
extern bool original_copy_tables_initialized_p (void);
extern void set_bb_original (basic_block, basic_block);
extern basic_block get_bb_original (basic_block);
extern void set_bb_copy (basic_block, basic_block);
extern basic_block get_bb_copy (basic_block);
void set_loop_copy (class loop *, class loop *);
class loop *get_loop_copy (class loop *);

/* Generic RAII class to allocate a bit from storage of integer type T.
   The allocated bit is accessible as mask with the single bit set
   via the conversion operator to T.  */

template <class T>
class auto_flag
{
public:
  /* static assert T is integer type of max HOST_WIDE_INT precision.  */
  auto_flag (T *sptr)
    {
      m_sptr = sptr;
      int free_bit = ffs_hwi (~*sptr);
      /* If there are no unset bits... */
      if (free_bit == 0)
	gcc_unreachable ();
      m_flag = HOST_WIDE_INT_1U << (free_bit - 1);
      /* ...or if T is signed and thus the complement is sign-extended,
         check if we ran out of bits.  We could spare us this bit
	 if we could use C++11 std::make_unsigned<T>::type to pass
	 ~*sptr to ffs_hwi.  */
      if (m_flag == 0)
	gcc_unreachable ();
      gcc_checking_assert ((*sptr & m_flag) == 0);
      *sptr |= m_flag;
    }
  ~auto_flag ()
    {
      gcc_checking_assert ((*m_sptr & m_flag) == m_flag);
      *m_sptr &= ~m_flag;
    }
  operator T () const { return m_flag; }
private:
  T *m_sptr;
  T m_flag;
};

/* RAII class to allocate an edge flag for temporary use.  You have
   to clear the flag from all edges when you are finished using it.  */

class auto_edge_flag : public auto_flag<int>
{
public:
  auto_edge_flag (function *fun)
    : auto_flag<int> (&fun->cfg->edge_flags_allocated) {}
};

/* RAII class to allocate a bb flag for temporary use.  You have
   to clear the flag from all edges when you are finished using it.  */
class auto_bb_flag : public auto_flag<int>
{
public:
  auto_bb_flag (function *fun)
    : auto_flag<int> (&fun->cfg->bb_flags_allocated) {}
};

#endif /* GCC_CFG_H */

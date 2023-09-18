/* Header file for SSA jump threading.
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

#ifndef GCC_TREE_SSA_THREADEDGE_H
#define GCC_TREE_SSA_THREADEDGE_H

// Class used to maintain path state in the jump threader and pass it
// to the jump threader simplifier.

class jt_state
{
public:
  virtual ~jt_state () { }
  virtual void push (edge);
  virtual void pop ();
  virtual void register_equiv (tree dest, tree src, bool update_range);
  virtual void register_equivs_edge (edge e);
  virtual void register_equivs_stmt (gimple *, basic_block,
				     class jt_simplifier *);
  virtual void record_ranges_from_stmt (gimple *stmt, bool temporary);
  void get_path (vec<basic_block> &);
  void append_path (basic_block);
  void dump (FILE *);
  void debug ();

private:
  auto_vec<basic_block> m_blocks;
  static const basic_block BB_MARKER;
};

// Statement simplifier callback for the jump threader.

class jt_simplifier
{
public:
  virtual ~jt_simplifier () { }
  virtual tree simplify (gimple *, gimple *, basic_block, jt_state *) = 0;
};

class hybrid_jt_state : public jt_state
{
private:
  void register_equivs_stmt (gimple *, basic_block, jt_simplifier *) override
  {
    // Ranger has no need to simplify anything.
  }
};

class hybrid_jt_simplifier : public jt_simplifier
{
public:
  hybrid_jt_simplifier (class gimple_ranger *r, class path_range_query *q);
  tree simplify (gimple *stmt, gimple *, basic_block, jt_state *) override;

private:
  void compute_ranges_from_state (gimple *stmt, jt_state *);

  gimple_ranger *m_ranger;
  path_range_query *m_query;
  auto_vec<basic_block> m_path;
};

// This is the high level threader.  The entry point is
// thread_outgoing_edges(), which calculates and registers paths to be
// threaded.  When all candidates have been registered,
// thread_through_all_blocks() is called to actually change the CFG.

class jump_threader
{
public:
  jump_threader (jt_simplifier *, class jt_state *);
  ~jump_threader ();
  void thread_outgoing_edges (basic_block);
  void remove_jump_threads_including (edge_def *);
  bool thread_through_all_blocks (bool may_peel_loop_headers);

private:
  tree simplify_control_stmt_condition (edge, gimple *);
  tree simplify_control_stmt_condition_1 (edge,
					  gimple *,
					  tree op0,
					  tree_code cond_code,
					  tree op1,
					  unsigned limit);

  bool thread_around_empty_blocks (vec<class jump_thread_edge *> *path,
				   edge, bitmap visited);
  int thread_through_normal_block (vec<jump_thread_edge *> *path,
				   edge, bitmap visited);
  void thread_across_edge (edge);
  bool record_temporary_equivalences_from_phis (edge);
  gimple *record_temporary_equivalences_from_stmts_at_dest (edge);

  // Dummy condition to avoid creating lots of throw away statements.
  gcond *dummy_cond;

  class fwd_jt_path_registry *m_registry;
  jt_simplifier *m_simplifier;
  jt_state *m_state;
};

extern void propagate_threaded_block_debug_into (basic_block, basic_block);
extern bool single_succ_to_potentially_threadable_block (basic_block);

// ?? All this ssa_name_values stuff is the store of values for
// avail_exprs_stack and const_and_copies, so it really belongs in the
// jump_threader class.  However, it's probably not worth touching
// this, since all this windable state is slated to go with the
// ranger.
extern vec<tree> ssa_name_values;
#define SSA_NAME_VALUE(x) \
    (SSA_NAME_VERSION (x) < ssa_name_values.length () \
     ? ssa_name_values[SSA_NAME_VERSION (x)] \
     : NULL_TREE)
extern void set_ssa_name_value (tree, tree);

#endif /* GCC_TREE_SSA_THREADEDGE_H */

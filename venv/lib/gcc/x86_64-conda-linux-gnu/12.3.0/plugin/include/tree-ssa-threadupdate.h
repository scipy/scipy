/* Communication between registering jump thread requests and
   updating the SSA/CFG for jump threading.
   Copyright (C) 2013-2022 Free Software Foundation, Inc.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3, or (at your option)
any later version.

GCC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

#ifndef _TREE_SSA_THREADUPDATE_H
#define _TREE_SSA_THREADUPDATE_H 1

enum jump_thread_edge_type
{
  EDGE_START_JUMP_THREAD,
  EDGE_COPY_SRC_BLOCK,
  EDGE_COPY_SRC_JOINER_BLOCK,
  EDGE_NO_COPY_SRC_BLOCK
};

// We keep the registered jump threading opportunities in this
// vector as edge pairs (original_edge, target_edge).

class jump_thread_edge
{
public:
  jump_thread_edge (edge e, jump_thread_edge_type t) : e (e), type (t) {}

  edge e;
  jump_thread_edge_type type;
};

class jump_thread_path_allocator
{
public:
  jump_thread_path_allocator ();
  ~jump_thread_path_allocator ();
  jump_thread_edge *allocate_thread_edge (edge, jump_thread_edge_type);
  vec<jump_thread_edge *> *allocate_thread_path ();
private:
  DISABLE_COPY_AND_ASSIGN (jump_thread_path_allocator);
  obstack m_obstack;
};

// Abstract class for the jump thread registry.
//
// When all candidates have been registered with
// register_jump_thread(), thread_through_all_blocks() is called to
// update the CFG.

class jt_path_registry
{
public:
  jt_path_registry (bool backedge_threads);
  virtual ~jt_path_registry ();
  bool register_jump_thread (vec<jump_thread_edge *> *);
  bool thread_through_all_blocks (bool peel_loop_headers);
  void push_edge (vec<jump_thread_edge *> *path, edge, jump_thread_edge_type);
  vec<jump_thread_edge *> *allocate_thread_path ();
  void debug ();
protected:
  void debug_path (FILE *, int pathno);
  vec<vec<jump_thread_edge *> *> m_paths;
  unsigned long m_num_threaded_edges;
private:
  virtual bool update_cfg (bool peel_loop_headers) = 0;
  bool cancel_invalid_paths (vec<jump_thread_edge *> &path);
  jump_thread_path_allocator m_allocator;
  // True if threading through back edges is allowed.  This is only
  // allowed in the generic copier in the backward threader.
  bool m_backedge_threads;
  DISABLE_COPY_AND_ASSIGN (jt_path_registry);
};

// Forward threader path registry using a custom BB copier.

class fwd_jt_path_registry : public jt_path_registry
{
public:
  fwd_jt_path_registry ();
  ~fwd_jt_path_registry ();
  void remove_jump_threads_including (edge);
private:
  bool update_cfg (bool peel_loop_headers) override;
  void mark_threaded_blocks (bitmap threaded_blocks);
  bool thread_block_1 (basic_block, bool noloop_only, bool joiners);
  bool thread_block (basic_block, bool noloop_only);
  bool thread_through_loop_header (class loop *loop,
				   bool may_peel_loop_headers);
  class redirection_data *lookup_redirection_data (edge e, enum insert_option);

  hash_table<struct removed_edges> *m_removed_edges;

  // Main data structure to hold information for duplicates of BB.
  hash_table<redirection_data> *m_redirection_data;
};

// Backward threader path registry using a generic BB copier.

class back_jt_path_registry : public jt_path_registry
{
public:
  back_jt_path_registry ();
private:
  bool update_cfg (bool peel_loop_headers) override;
  void adjust_paths_after_duplication (unsigned curr_path_num);
  bool duplicate_thread_path (edge entry, edge exit, basic_block *region,
			      unsigned n_region, unsigned current_path_no);
  bool rewire_first_differing_edge (unsigned path_num, unsigned edge_num);
};

// Rather than search all the edges in jump thread paths each time DOM
// is able to simply if control statement, we build a hash table with
// the deleted edges.  We only care about the address of the edge, not
// its contents.
struct removed_edges : nofree_ptr_hash<edge_def>
{
  static hashval_t hash (edge e) { return htab_hash_pointer (e); }
  static bool equal (edge e1, edge e2) { return e1 == e2; }
};

extern unsigned int estimate_threading_killed_stmts (basic_block);

enum bb_dom_status
{
  /* BB does not dominate latch of the LOOP.  */
  DOMST_NONDOMINATING,
  /* The LOOP is broken (there is no path from the header to its latch.  */
  DOMST_LOOP_BROKEN,
  /* BB dominates the latch of the LOOP.  */
  DOMST_DOMINATING
};

enum bb_dom_status determine_bb_domination_status (class loop *, basic_block);

// In tree-ssa-dom.cc.
extern void free_dom_edge_info (edge);

#endif

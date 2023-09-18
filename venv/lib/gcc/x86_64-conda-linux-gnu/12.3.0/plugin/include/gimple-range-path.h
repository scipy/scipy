/* Header file for jump threading path solver.
   Copyright (C) 2021-2022 Free Software Foundation, Inc.
   Contributed by Aldy Hernandez <aldyh@redhat.com>.

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

#ifndef GCC_TREE_SSA_THREADSOLVER_H
#define GCC_TREE_SSA_THREADSOLVER_H

// This class is a basic block path solver.  Given a set of BBs
// indicating a path through the CFG, range_of_expr and range_of_stmt
// will calculate the range of an SSA or STMT as if the BBs in the
// path would have been executed in order.
//
// Note that the blocks are in reverse order, thus the exit block is
// path[0].

class path_range_query : public range_query
{
public:
  path_range_query (bool resolve = true, class gimple_ranger *ranger = NULL);
  virtual ~path_range_query ();
  void compute_ranges (const vec<basic_block> &,
		       const bitmap_head *imports = NULL);
  void compute_ranges (edge e);
  void compute_imports (bitmap imports, basic_block exit);
  bool range_of_expr (irange &r, tree name, gimple * = NULL) override;
  bool range_of_stmt (irange &r, gimple *, tree name = NULL) override;
  bool unreachable_path_p ();
  void dump (FILE *) override;
  void debug ();

private:
  bool internal_range_of_expr (irange &r, tree name, gimple *);
  bool defined_outside_path (tree name);
  void range_on_path_entry (irange &r, tree name);
  path_oracle *get_path_oracle () { return (path_oracle *)m_oracle; }

  // Cache manipulation.
  void set_cache (const irange &r, tree name);
  bool get_cache (irange &r, tree name);
  void clear_cache (tree name);

  // Methods to compute ranges for the given path.
  bool range_defined_in_block (irange &, tree name, basic_block bb);
  void compute_ranges_in_block (basic_block bb);
  void compute_ranges_in_phis (basic_block bb);
  void adjust_for_non_null_uses (basic_block bb);
  void ssa_range_in_phi (irange &r, gphi *phi);
  void compute_outgoing_relations (basic_block bb, basic_block next);
  void compute_phi_relations (basic_block bb, basic_block prev);
  void maybe_register_phi_relation (gphi *, edge e);
  bool add_to_imports (tree name, bitmap imports);
  bool import_p (tree name);
  bool ssa_defined_in_bb (tree name, basic_block bb);
  bool relations_may_be_invalidated (edge);

  // Path navigation.
  void set_path (const vec<basic_block> &);
  basic_block entry_bb () { return m_path[m_path.length () - 1]; }
  basic_block exit_bb ()  { return m_path[0]; }
  basic_block curr_bb ()  { return m_path[m_pos]; }
  basic_block prev_bb ()  { return m_path[m_pos + 1]; }
  basic_block next_bb ()  { return m_path[m_pos - 1]; }
  bool at_entry ()	  { return m_pos == m_path.length () - 1; }
  bool at_exit ()	  { return m_pos == 0; }
  void move_next ()	  { --m_pos; }

  // Range cache for SSA names.
  ssa_global_cache *m_cache;

  // Set for each SSA that has an active entry in the cache.
  bitmap m_has_cache_entry;

  // Path being analyzed.
  auto_vec<basic_block> m_path;

  auto_bitmap m_imports;
  gimple_ranger *m_ranger;
  non_null_ref m_non_null;

  // Current path position.
  unsigned m_pos;

  // Use ranger to resolve anything not known on entry.
  bool m_resolve;

  // Set if there were any undefined expressions while pre-calculating path.
  bool m_undefined_path;

  // True if m_ranger was allocated in this class and must be freed at
  // destruction.
  bool m_alloced_ranger;
};

#endif // GCC_TREE_SSA_THREADSOLVER_H

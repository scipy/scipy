/* Header file for gimple ranger SSA cache.
   Copyright (C) 2017-2022 Free Software Foundation, Inc.
   Contributed by Andrew MacLeod <amacleod@redhat.com>.

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

#ifndef GCC_SSA_RANGE_CACHE_H
#define GCC_SSA_RANGE_CACHE_H

#include "gimple-range-gori.h" 

// Class used to track non-null references of an SSA name.  A vector
// of bitmaps indexed by SSA name is maintained.  When indexed by
// basic block, an on-bit indicates there is a non-null dereference
// for that SSA in that block.

class non_null_ref
{
public:
  non_null_ref ();
  ~non_null_ref ();
  bool non_null_deref_p (tree name, basic_block bb, bool search_dom = true);
  bool adjust_range (irange &r, tree name, basic_block bb,
		     bool search_dom = true);
  bool set_nonnull (basic_block bb, tree name);
private:
  vec <bitmap> m_nn;
  void process_name (tree name);
  bitmap_obstack m_bitmaps;
};

// If NAME has a non-null dereference in block BB, adjust R with the
// non-zero information from non_null_deref_p, and return TRUE.  If
// SEARCH_DOM is true, non_null_deref_p should search the dominator tree.

inline bool
non_null_ref::adjust_range (irange &r, tree name, basic_block bb,
			    bool search_dom)
{
  // Non-call exceptions mean we could throw in the middle of the
  // block, so just punt on those for now.
  if (cfun->can_throw_non_call_exceptions)
    return false;
  // We only care about the null / non-null property of pointers.
  if (!POINTER_TYPE_P (TREE_TYPE (name)))
    return false;
  if (r.undefined_p () || r.lower_bound () != 0 || r.upper_bound () == 0)
    return false;
  // Check if pointers have any non-null dereferences.
  if (non_null_deref_p (name, bb, search_dom))
    {
      // Remove zero from the range.
      unsigned prec = TYPE_PRECISION (TREE_TYPE (name));
      r.intersect (wi::one (prec), wi::max_value (prec, UNSIGNED));
      return true;
    }
  return false;
}

// This class manages a vector of pointers to ssa_block ranges.  It
// provides the basis for the "range on entry" cache for all
// SSA names.

class block_range_cache
{
public:
  block_range_cache ();
  ~block_range_cache ();

  bool set_bb_range (tree name, const_basic_block bb, const irange &r);
  bool get_bb_range (irange &r, tree name, const_basic_block bb);
  bool bb_range_p (tree name, const_basic_block bb);

  void dump (FILE *f);
  void dump (FILE *f, basic_block bb, bool print_varying = true);
private:
  vec<class ssa_block_ranges *> m_ssa_ranges;
  ssa_block_ranges &get_block_ranges (tree name);
  ssa_block_ranges *query_block_ranges (tree name);
  irange_allocator *m_irange_allocator;
  bitmap_obstack m_bitmaps;
};

// This global cache is used with the range engine as markers for what
// has been visited during this incarnation.  Once the ranger evaluates
// a name, it is typically not re-evaluated again.

class ssa_global_cache
{
public:
  ssa_global_cache ();
  ~ssa_global_cache ();
  bool get_global_range (irange &r, tree name) const;
  bool set_global_range (tree name, const irange &r);
  void clear_global_range (tree name);
  void clear ();
  void dump (FILE *f = stderr);
private:
  vec<irange *> m_tab;
  class irange_allocator *m_irange_allocator;
};

// This class provides all the caches a global ranger may need, and makes 
// them available for gori-computes to query so outgoing edges can be
// properly calculated.

class ranger_cache : public range_query
{
public:
  ranger_cache (int not_executable_flag);
  ~ranger_cache ();

  virtual bool range_of_expr (irange &r, tree name, gimple *stmt);
  virtual bool range_on_edge (irange &r, edge e, tree expr);
  bool block_range (irange &r, basic_block bb, tree name, bool calc = true);
  bool range_from_dom (irange &r, tree name, basic_block bb);

  bool get_global_range (irange &r, tree name) const;
  bool get_global_range (irange &r, tree name, bool &current_p);
  void set_global_range (tree name, const irange &r);

  void propagate_updated_value (tree name, basic_block bb);

  void block_apply_nonnull (gimple *s);
  void update_to_nonnull (basic_block bb, tree name);
  non_null_ref m_non_null;
  gori_compute m_gori;

  void dump_bb (FILE *f, basic_block bb);
  virtual void dump (FILE *f) OVERRIDE;
private:
  ssa_global_cache m_globals;
  block_range_cache m_on_entry;
  class temporal_cache *m_temporal;
  void fill_block_cache (tree name, basic_block bb, basic_block def_bb);
  void propagate_cache (tree name);

  void range_of_def (irange &r, tree name, basic_block bb = NULL);
  void entry_range (irange &r, tree expr, basic_block bb);
  void exit_range (irange &r, tree expr, basic_block bb);

  vec<basic_block> m_workback;
  class update_list *m_update;
};

#endif // GCC_SSA_RANGE_CACHE_H

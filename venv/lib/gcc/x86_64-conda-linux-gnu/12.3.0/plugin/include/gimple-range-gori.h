/* Header file for gimple range GORI structures.
   Copyright (C) 2017-2022 Free Software Foundation, Inc.
   Contributed by Andrew MacLeod <amacleod@redhat.com>
   and Aldy Hernandez <aldyh@redhat.com>.

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

#ifndef GCC_GIMPLE_RANGE_GORI_H
#define GCC_GIMPLE_RANGE_GORI_H

// RANGE_DEF_CHAIN is used to determine which SSA names in a block can
// have range information calculated for them, and what the
// dependencies on each other are.

class range_def_chain
{
public:
  range_def_chain ();
  ~range_def_chain ();
  tree depend1 (tree name) const;
  tree depend2 (tree name) const;
  bool in_chain_p (tree name, tree def);
  bool chain_import_p (tree name, tree import);
  void register_dependency (tree name, tree ssa1, basic_block bb = NULL);
  void dump (FILE *f, basic_block bb, const char *prefix = NULL);
protected:
  bool has_def_chain (tree name);
  bool def_chain_in_bitmap_p (tree name, bitmap b);
  void add_def_chain_to_bitmap (bitmap b, tree name);
  bitmap get_def_chain (tree name);
  bitmap get_imports (tree name);
  bitmap_obstack m_bitmaps;
private:
  struct rdc {
   tree ssa1;		// First direct dependency
   tree ssa2;		// Second direct dependency
   bitmap bm;		// All dependencies
   bitmap m_import;
  };
  vec<rdc> m_def_chain;	// SSA_NAME : def chain components.
  void set_import (struct rdc &data, tree imp, bitmap b);
  int m_logical_depth;
};

// Return the first direct dependency for NAME, if there is one.
// Direct dependencies are those which occur on the defintion statement.
// Only the first 2 such names are cached.

inline tree
range_def_chain::depend1 (tree name) const
{
  unsigned v = SSA_NAME_VERSION (name);
  if (v >= m_def_chain.length ())
    return NULL_TREE;
  return m_def_chain[v].ssa1;
}

// Return the second direct dependency for NAME, if there is one.

inline tree
range_def_chain::depend2 (tree name) const
{
  unsigned v = SSA_NAME_VERSION (name);
  if (v >= m_def_chain.length ())
    return NULL_TREE;
  return m_def_chain[v].ssa2;
}

// GORI_MAP is used to accumulate what SSA names in a block can
// generate range information, and provides tools for the block ranger
// to enable it to efficiently calculate these ranges.

class gori_map : public range_def_chain
{
public:
  gori_map ();
  ~gori_map ();

  bool is_export_p (tree name, basic_block bb = NULL);
  bool is_import_p (tree name, basic_block bb);
  bitmap exports (basic_block bb);
  bitmap imports (basic_block bb);
  void set_range_invariant (tree name);

  void dump (FILE *f);
  void dump (FILE *f, basic_block bb, bool verbose = true);
private:
  vec<bitmap> m_outgoing;	// BB: Outgoing ranges calculatable on edges
  vec<bitmap> m_incoming;	// BB: Incoming ranges which can affect exports.
  bitmap m_maybe_variant;	// Names which might have outgoing ranges.
  void maybe_add_gori (tree name, basic_block bb);
  void calculate_gori (basic_block bb);
};


// This class is used to determine which SSA_NAMES can have ranges
// calculated for them on outgoing edges from basic blocks.  This represents
// ONLY the effect of the basic block edge->src on a range.
//
// There are 2 primary entry points:
//
// has_edge_range_p (tree name, edge e)
//   returns true if the outgoing edge *may* be able to produce range
//   information for ssa_name NAME on edge E.
//   FALSE is returned if this edge does not affect the range of NAME.
//   if no edge is specified, return TRUE if name may have a value calculated
//   on *ANY* edge that has been seen.  FALSE indicates that the global value
//   is applicable everywhere that has been processed.
//
// outgoing_edge_range_p (irange &range, edge e, tree name)
//   Actually does the calculation of RANGE for name on E
//   This represents application of whatever static range effect edge E
//   may have on NAME, not any cumulative effect.

// There are also some internal APIs
//
// ssa_range_in_bb ()  is an internal routine which is used to start any
// calculation chain using SSA_NAMES which come from outside the block. ie
//      a_2 = b_4 - 8
//      if (a_2 < 30)
// on the true edge, a_2 is known to be [0, 29]
// b_4 can be calculated as [8, 37]
// during this calculation, b_4 is considered an "import" and ssa_range_in_bb
// is queried for a starting range which is used in the calculation.
// A default value of VARYING provides the raw static info for the edge.
//
// If there is any known range for b_4 coming into this block, it can refine
// the results.  This allows for cascading results to be propogated.
// if b_4 is [100, 200] on entry to the block, feeds into the calculation
// of a_2 = [92, 192], and finally on the true edge the range would be 
// an empty range [] because it is not possible for the true edge to be taken.
//
// expr_range_in_bb is simply a wrapper which calls ssa_range_in_bb for 
// SSA_NAMES and otherwise simply calculates the range of the expression.
//
// The constructor takes a flag value to use on edges to check for the
// NON_EXECUTABLE_EDGE property.  The zero default means no flag is checked.
// All value requests from NON_EXECUTABLE_EDGE edges are returned UNDEFINED.
//
// The remaining routines are internal use only.

class gori_compute : public gori_map
{
public:
  gori_compute (int not_executable_flag = 0);
  bool outgoing_edge_range_p (irange &r, edge e, tree name, range_query &q);
  bool condexpr_adjust (irange &r1, irange &r2, gimple *s, tree cond, tree op1,
			tree op2, fur_source &src);
  bool has_edge_range_p (tree name, basic_block bb = NULL);
  bool has_edge_range_p (tree name, edge e);
  void dump (FILE *f);
private:
  bool may_recompute_p (tree name, edge e);
  bool may_recompute_p (tree name, basic_block bb = NULL);
  bool compute_operand_range (irange &r, gimple *stmt, const irange &lhs,
			      tree name, class fur_source &src);
  bool compute_operand_range_switch (irange &r, gswitch *s, const irange &lhs,
				     tree name, fur_source &src);
  bool compute_operand1_range (irange &r, gimple *stmt, const irange &lhs,
			       tree name, fur_source &src);
  bool compute_operand2_range (irange &r, gimple *stmt, const irange &lhs,
			       tree name, fur_source &src);
  bool compute_operand1_and_operand2_range (irange &r, gimple *stmt,
					    const irange &lhs, tree name,
					    fur_source &src);
  void compute_logical_operands (irange &true_range, irange &false_range,
				 gimple *stmt, const irange &lhs,
				 tree name, fur_source &src, tree op,
				 bool op_in_chain);
  bool logical_combine (irange &r, enum tree_code code, const irange &lhs,
			const irange &op1_true, const irange &op1_false,
			const irange &op2_true, const irange &op2_false);
  int_range<2> m_bool_zero;	// Boolean false cached.
  int_range<2> m_bool_one;	// Boolean true cached.

  gimple_outgoing_range outgoing;	// Edge values for COND_EXPR & SWITCH_EXPR.
  range_tracer tracer;
  int m_not_executable_flag;
};

// These routines provide a GIMPLE interface to the range-ops code.
extern bool gimple_range_calc_op1 (irange &r, const gimple *s,
				   const irange &lhs_range);
extern bool gimple_range_calc_op1 (irange &r, const gimple *s,
				   const irange &lhs_range,
				   const irange &op2_range);
extern bool gimple_range_calc_op2 (irange &r, const gimple *s,
				   const irange &lhs_range,
				   const irange &op1_range);

// For each name that is an import into BB's exports..
#define FOR_EACH_GORI_IMPORT_NAME(gori, bb, name)			\
  for (gori_export_iterator iter ((gori).imports ((bb)));	\
       ((name) = iter.get_name ());				\
       iter.next ())

// For each name possibly exported from block BB.
#define FOR_EACH_GORI_EXPORT_NAME(gori, bb, name)		\
  for (gori_export_iterator iter ((gori).exports ((bb)));	\
       ((name) = iter.get_name ());				\
       iter.next ())

// Used to assist with iterating over the GORI export list in various ways
class gori_export_iterator {
public:
  gori_export_iterator (bitmap b);
  void next ();
  tree get_name ();
protected:
  bitmap bm;
  bitmap_iterator bi;
  unsigned y;
};

#endif // GCC_GIMPLE_RANGE_GORI_H

/* Header file for the value range relational processing.
   Copyright (C) 2020-2022 Free Software Foundation, Inc.
   Contributed by Andrew MacLeod <amacleod@redhat.com>

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

#ifndef GCC_VALUE_RELATION_H
#define GCC_VALUE_RELATION_H


// This file provides access to a relation oracle which can be used to
// maintain and query relations and equivalences between SSA_NAMES.
//
// The general range_query object provided in value-query.h provides
// access to an oracle, if one is available, via the oracle() method.
// Thre are also a couple of access routines provided, which even if there is
// no oracle, will return the default VREL_NONE no relation.
//
// Typically, when a ranger object is active, there will be an oracle, and
// any information available can be directly queried.  Ranger also sets and
// utilizes the relation information to enhance it's range calculations, this
// is totally transparent to the client, and they are free to make queries.
//
//
// relation_kind is a typedef of enum tree_code, but has restricted range
// and a couple of extra values.
//
// A query is made requesting the relation between SSA1 and SSA@ in a basic
// block, or on an edge, the possible return values are:
//
//  EQ_EXPR, NE_EXPR, LT_EXPR, LE_EXPR, GT_EXPR, and GE_EXPR mean the same.
//  VREL_NONE : No relation between the 2 names.
//  VREL_EMPTY : Impossible relation (ie, A < B && A > B produces VREL_EMPTY.
//
// The oracle maintains EQ_EXPR relations with equivalency sets, so if a
// relation comes back EQ_EXPR, it is also possible to query the set of
// equivlaencies.  These are basically bitmaps over ssa_names.
//
// Relations are maintained via the dominace trees and are optimized assuming
// they are registered in dominance order.   When a new relation is added, it
// is intersected with whatever existing relation exists in the dominance tree
// and registered at the specified block.


// Rather than introduce a new enumerated type for relations, we can use the
// existing tree_codes for relations, plus add a couple of #defines for
// the other cases.  These codes are arranged such that VREL_NONE is the first
// code, and all the rest are contiguous.

typedef enum tree_code relation_kind;

#define VREL_NONE		TRUTH_NOT_EXPR
#define VREL_EMPTY		LTGT_EXPR

// General relation kind transformations.
relation_kind relation_union (relation_kind r1, relation_kind r2);
relation_kind relation_intersect (relation_kind r1, relation_kind r2);
relation_kind relation_negate (relation_kind r);
relation_kind relation_swap (relation_kind r);
void print_relation (FILE *f, relation_kind rel);


class relation_oracle
{
public:
  virtual ~relation_oracle () { }
  // register a relation between 2 ssa names at a stmt.
  void register_stmt (gimple *, relation_kind, tree, tree);
  // register a relation between 2 ssa names on an edge.
  void register_edge (edge, relation_kind, tree, tree);

  // Return equivalency set for an SSA name in a basic block.
  virtual const_bitmap equiv_set (tree, basic_block) = 0;
  // register a relation between 2 ssa names in a basic block.
  virtual void register_relation (basic_block, relation_kind, tree, tree) = 0;
  // Query for a relation between two ssa names in a basic block.
  virtual relation_kind query_relation (basic_block, tree, tree) = 0;
  // Query for a relation between two equivalency stes in a basic block.
  virtual relation_kind query_relation (basic_block, const_bitmap,
					const_bitmap) = 0;

  virtual void dump (FILE *, basic_block) const = 0;
  virtual void dump (FILE *) const = 0;
  void debug () const;
protected:
  void valid_equivs (bitmap b, const_bitmap equivs, basic_block bb);
};

// This class represents an equivalency set, and contains a link to the next
// one in the list to be searched.

class equiv_chain
{
public:
  bitmap m_names;		// ssa-names in equiv set.
  basic_block m_bb;		// Block this belongs to
  equiv_chain *m_next;		// Next in block list.
  void dump (FILE *f) const;	// Show names in this list.
  equiv_chain *find (unsigned ssa);
};

// The equivalency oracle maintains equivalencies using the dominator tree.
// Equivalencies apply to an entire basic block.  Equivalencies on edges
// can be represented only on edges whose destination is a single-pred block,
// and the equivalence is simply applied to that succesor block.

class equiv_oracle : public relation_oracle
{
public:
  equiv_oracle ();
  ~equiv_oracle ();

  const_bitmap equiv_set (tree ssa, basic_block bb);
  void register_relation (basic_block bb, relation_kind k, tree ssa1,
			  tree ssa2);

  relation_kind query_relation (basic_block, tree, tree);
  relation_kind query_relation (basic_block, const_bitmap, const_bitmap);
  void dump (FILE *f, basic_block bb) const;
  void dump (FILE *f) const;

protected:
  bitmap_obstack m_bitmaps;
  struct obstack m_chain_obstack;
private:
  bitmap m_equiv_set;	// Index by ssa-name. true if an equivalence exists.
  vec <equiv_chain *> m_equiv;	// Index by BB.  list of equivalences.
  vec <bitmap> m_self_equiv;  // Index by ssa-name, self equivalency set.

  void limit_check (basic_block bb = NULL);
  equiv_chain *find_equiv_block (unsigned ssa, int bb) const;
  equiv_chain *find_equiv_dom (tree name, basic_block bb) const;

  bitmap register_equiv (basic_block bb, unsigned v, equiv_chain *equiv_1);
  bitmap register_equiv (basic_block bb, equiv_chain *equiv_1,
			 equiv_chain *equiv_2);
  void register_initial_def (tree ssa);
  void add_equiv_to_block (basic_block bb, bitmap equiv);
};

// Summary block header for relations.

class relation_chain_head
{
public:
  bitmap m_names;		// ssa_names with relations in this block.
  class relation_chain *m_head; // List of relations in block.
  int m_num_relations;		// Number of relations in block.
  relation_kind find_relation (const_bitmap b1, const_bitmap b2) const;
};

// A relation oracle maintains a set of relations between ssa_names using the
// dominator tree structures.  Equivalencies are considered a subset of
// a general relation and maintained by an equivalence oracle by transparently
// passing any EQ_EXPR relations to it.
// Relations are handled at the basic block level.  All relations apply to
// an entire block, and are thus kept in a summary index by block.
// Similar to the equivalence oracle, edges are handled by applying the
// relation to the destination block of the edge, but ONLY if that block
// has a single successor.  For now.

class dom_oracle : public equiv_oracle
{
public:
  dom_oracle ();
  ~dom_oracle ();

  void register_relation (basic_block bb, relation_kind k, tree op1, tree op2);

  relation_kind query_relation (basic_block bb, tree ssa1, tree ssa2);
  relation_kind query_relation (basic_block bb, const_bitmap b1,
				   const_bitmap b2);

  void dump (FILE *f, basic_block bb) const;
  void dump (FILE *f) const;
private:
  bitmap m_tmp, m_tmp2;
  bitmap m_relation_set;  // Index by ssa-name. True if a relation exists
  vec <relation_chain_head> m_relations;  // Index by BB, list of relations.
  relation_kind find_relation_block (unsigned bb, const_bitmap b1,
				     const_bitmap b2) const;
  relation_kind find_relation_block (int bb, unsigned v1, unsigned v2,
				     relation_chain **obj = NULL) const;
  relation_kind find_relation_dom (basic_block bb, unsigned v1, unsigned v2) const;
  relation_chain *set_one_relation (basic_block bb, relation_kind k, tree op1,
				    tree op2);
  void register_transitives (basic_block, const class value_relation &);

};

// A path_oracle implements relations in a list.  The only sense of ordering
// is the latest registered relation is the first found during a search.
// It can be constructed with an optional "root" oracle which will be used
// to look up any relations not found in the list.
// This allows the client to walk paths starting at some block and register
// and query relations along that path, ignoring other edges.
//
// For registering a relation, a query if made of the root oracle if there is
// any known relationship at block BB, and it is combined with this new
// relation and entered in the list.
//
// Queries are resolved by looking first in the list, and only if nothing is
// found is the root oracle queried at block BB.
//
// reset_path is used to clear all locally registered paths to initial state.

class path_oracle : public relation_oracle
{
public:
  path_oracle (relation_oracle *oracle = NULL);
  ~path_oracle ();
  const_bitmap equiv_set (tree, basic_block);
  void register_relation (basic_block, relation_kind, tree, tree);
  void killing_def (tree);
  relation_kind query_relation (basic_block, tree, tree);
  relation_kind query_relation (basic_block, const_bitmap, const_bitmap);
  void reset_path ();
  void set_root_oracle (relation_oracle *oracle) { m_root = oracle; }
  void dump (FILE *, basic_block) const;
  void dump (FILE *) const;
private:
  void register_equiv (basic_block bb, tree ssa1, tree ssa2);
  equiv_chain m_equiv;
  relation_chain_head m_relations;
  relation_oracle *m_root;
  bitmap m_killed_defs;

  bitmap_obstack m_bitmaps;
  struct obstack m_chain_obstack;
};
#endif  /* GCC_VALUE_RELATION_H */

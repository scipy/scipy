/* Data structures and function declarations for the SSA value propagation
   engine.
   Copyright (C) 2004-2022 Free Software Foundation, Inc.
   Contributed by Diego Novillo <dnovillo@redhat.com>

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

#ifndef _TREE_SSA_PROPAGATE_H
#define _TREE_SSA_PROPAGATE_H 1

#include "value-query.h"

/* If SIM_P is true, statement S will be simulated again.  */

static inline void
prop_set_simulate_again (gimple *s, bool visit_p)
{
  gimple_set_visited (s, visit_p);
}

/* Return true if statement T should be simulated again.  */

static inline bool
prop_simulate_again_p (gimple *s)
{
  return gimple_visited_p (s);
}

/* Lattice values used for propagation purposes.  Specific instances
   of a propagation engine must return these values from the statement
   and PHI visit functions to direct the engine.  */
enum ssa_prop_result {
    /* The statement produces nothing of interest.  No edges will be
       added to the work lists.  */
    SSA_PROP_NOT_INTERESTING,

    /* The statement produces an interesting value.  The set SSA_NAMEs
       returned by SSA_PROP_VISIT_STMT should be added to
       INTERESTING_SSA_EDGES.  If the statement being visited is a
       conditional jump, SSA_PROP_VISIT_STMT should indicate which edge
       out of the basic block should be marked executable.  */
    SSA_PROP_INTERESTING,

    /* The statement produces a varying (i.e., useless) value and
       should not be simulated again.  If the statement being visited
       is a conditional jump, all the edges coming out of the block
       will be considered executable.  */
    SSA_PROP_VARYING
};


extern void move_ssa_defining_stmt_for_defs (gimple *, gimple *);
extern bool stmt_makes_single_store (gimple *);
extern bool may_propagate_copy (tree, tree, bool = false);
extern bool may_propagate_copy_into_stmt (gimple *, tree);
extern bool may_propagate_copy_into_asm (tree);
extern void propagate_value (use_operand_p, tree);
extern void replace_exp (use_operand_p, tree);
extern void propagate_tree_value (tree *, tree);
extern void propagate_tree_value_into_stmt (gimple_stmt_iterator *, tree);

/* Public interface into the SSA propagation engine.  Clients should inherit
   from this class and provide their own visitors.  */

class ssa_propagation_engine
{
 public:

  virtual ~ssa_propagation_engine (void) { }

  /* Virtual functions the clients must provide to visit statements
     and phi nodes respectively.  */
  virtual enum ssa_prop_result visit_stmt (gimple *, edge *, tree *) = 0;
  virtual enum ssa_prop_result visit_phi (gphi *) = 0;

  /* Main interface into the propagation engine.  */
  void ssa_propagate (void);

 private:
  /* Internal implementation details.  */
  void simulate_stmt (gimple *stmt);
  void simulate_block (basic_block);
};

class substitute_and_fold_engine : public value_query
{
 public:
  substitute_and_fold_engine (bool fold_all_stmts = false)
    : fold_all_stmts (fold_all_stmts) { }
  virtual ~substitute_and_fold_engine (void) { }
  virtual bool fold_stmt (gimple_stmt_iterator *) { return false; }

  bool substitute_and_fold (basic_block = NULL);
  bool replace_uses_in (gimple *);
  bool replace_phi_args_in (gphi *);

  virtual void pre_fold_bb (basic_block) { }
  virtual void post_fold_bb (basic_block) { }
  virtual void pre_fold_stmt (gimple *) { }
  virtual void post_new_stmt (gimple *) { }

  bool propagate_into_phi_args (basic_block);

  /* Users like VRP can set this when they want to perform
     folding for every propagation.  */
  bool fold_all_stmts;
};

#endif /* _TREE_SSA_PROPAGATE_H  */

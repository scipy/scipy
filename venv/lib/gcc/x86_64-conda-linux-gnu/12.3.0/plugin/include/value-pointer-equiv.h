/* Header file for the context-aware pointer equivalence tracker.
   Copyright (C) 2020-2022 Free Software Foundation, Inc.
   Contributed by Aldy Hernandez <aldyh@redhat.com>

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

#ifndef GCC_VALUE_POINTER_EQUIV_H
#define GCC_VALUE_POINTER_EQUIV_H

// Simple context-aware pointer equivalency analyzer that returns what
// a pointer SSA name is equivalent to at a given point during a walk
// of the IL.
//
// Note that global equivalency take priority over conditional
// equivalency.  That is, p = &q takes priority over a later p == &t.
//
// This class is meant to be called during a DOM walk.

class pointer_equiv_analyzer
{
public:
  pointer_equiv_analyzer (gimple_ranger *r);
  ~pointer_equiv_analyzer ();
  void enter (basic_block);
  void leave (basic_block);
  void visit_stmt (gimple *stmt);
  tree get_equiv (tree ssa);

private:
  void visit_edge (edge e);
  tree get_equiv_expr (tree_code code, tree expr);
  void set_global_equiv (tree ssa, tree pointee);
  void set_cond_equiv (tree ssa, tree pointee);

  gimple_ranger *m_ranger;
  // Global pointer equivalency indexed by SSA_NAME_VERSION.
  auto_vec<tree> m_global_points;
  // Conditional pointer equivalency.
  class ssa_equiv_stack *m_cond_points;
};

inline bool
supported_pointer_equiv_p (tree expr)
{
  return TREE_CODE (expr) == SSA_NAME && POINTER_TYPE_P (TREE_TYPE (expr));
}

#endif // GCC_VALUE_POINTER_EQUIV_H

/* Array bounds checking.
   Copyright (C) 2020-2022 Free Software Foundation, Inc.

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

#ifndef GCC_GIMPLE_ARRAY_BOUNDS_H
#define GCC_GIMPLE_ARRAY_BOUNDS_H

#include "pointer-query.h"

class array_bounds_checker
{
  friend class check_array_bounds_dom_walker;

public:
  array_bounds_checker (struct function *, range_query *);
  void check ();

private:
  static tree check_array_bounds (tree *tp, int *walk_subtree, void *data);
  bool check_array_ref (location_t, tree, gimple *, bool ignore_off_by_one);
  bool check_mem_ref (location_t, tree, bool ignore_off_by_one);
  void check_addr_expr (location_t, tree, gimple *);
  const value_range *get_value_range (const_tree op, gimple *);

  /* Current function.  */
  struct function *fun;
  /* A pointer_query object to store information about pointers and
     their targets in.  */
  pointer_query m_ptr_qry;
  /* Current statement.  */
  gimple *m_stmt;
};

#endif // GCC_GIMPLE_ARRAY_BOUNDS_H

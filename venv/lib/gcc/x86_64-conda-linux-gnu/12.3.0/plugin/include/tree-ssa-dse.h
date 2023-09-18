/* Support routines for dead store elimination. 
   Copyright (C) 2019-2022 Free Software Foundation, Inc.

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

#ifndef GCC_TREE_SSA_DSE_H
#define GCC_TREE_SSA_DSE_H

/* Return value from dse_classify_store */
enum dse_store_status
{
  DSE_STORE_LIVE,
  DSE_STORE_MAYBE_PARTIAL_DEAD,
  DSE_STORE_DEAD
};

dse_store_status dse_classify_store (ao_ref *, gimple *, bool, sbitmap,
				     bool * = NULL, tree = NULL);

void delete_dead_or_redundant_assignment (gimple_stmt_iterator *, const char *,
					  bitmap = NULL, bitmap = NULL);

#endif   /* GCC_TREE_SSA_DSE_H  */

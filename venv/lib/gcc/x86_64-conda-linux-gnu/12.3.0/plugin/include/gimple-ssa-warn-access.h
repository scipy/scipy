/* Pass to detect and issue warnings for invalid accesses, including
   invalid or mismatched allocation/deallocation calls.

   Copyright (C) 2020-2022 Free Software Foundation, Inc.
   Contributed by Martin Sebor <msebor@redhat.com>.

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

#ifndef GCC_GIMPLE_SSA_WARN_ACCESS_H
#define GCC_GIMPLE_SSA_WARN_ACCESS_H

extern bool check_nul_terminated_array (tree, tree, tree = NULL_TREE);
extern void warn_string_no_nul (location_t, gimple *, const char *, tree,
				tree, tree = NULL_TREE, bool = false,
				const wide_int[2] = NULL);
extern void warn_string_no_nul (location_t, tree, const char *, tree,
				tree, tree = NULL_TREE, bool = false,
				const wide_int[2] = NULL);
extern tree unterminated_array (tree, tree * = NULL, bool * = NULL);

extern bool maybe_warn_nonstring_arg (tree, gimple *);
extern bool maybe_warn_nonstring_arg (tree, tree);

class access_data;
extern bool maybe_warn_for_bound (opt_code, location_t, gimple *, tree,
				  tree[2], tree, const access_data * = NULL);
extern bool maybe_warn_for_bound (opt_code, location_t, tree, tree,
				  tree[2], tree, const access_data * = NULL);

class access_data;
extern bool check_access (tree, tree, tree, tree, tree, access_mode,
			  const access_data * = NULL);

#endif   // GCC_GIMPLE_SSA_WARN_ACCESS_H

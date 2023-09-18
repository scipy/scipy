/* Traits for hashing trees.
   Copyright (C) 2014-2022 Free Software Foundation, Inc.

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

#ifndef tree_hash_traits_h
#define tree_hash_traits_h

/* Hash for trees based on operand_equal_p.  */
struct tree_operand_hash : ggc_ptr_hash <tree_node>
{
  static inline hashval_t hash (const value_type &);
  static inline bool equal (const value_type &,
			    const compare_type &);
};

inline hashval_t
tree_operand_hash::hash (const value_type &t)
{
  return iterative_hash_expr (t, 0);
}

inline bool
tree_operand_hash::equal (const value_type &t1,
			  const compare_type &t2)
{
  return operand_equal_p (t1, t2, 0);
}

#endif

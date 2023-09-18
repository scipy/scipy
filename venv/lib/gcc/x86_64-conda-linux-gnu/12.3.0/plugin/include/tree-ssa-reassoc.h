/* Reassociation for trees.
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

#ifndef GCC_SSA_REASSOC_H
#define GCC_SSA_REASSOC_H

/* Operator, rank pair.  */
struct operand_entry
{
  unsigned int rank;
  unsigned int id;
  tree op;
  unsigned int count;
  gimple *stmt_to_insert;
};

struct range_entry
{
  tree exp;
  tree low;
  tree high;
  bool in_p;
  bool strict_overflow_p;
  unsigned int idx, next;
};

void dump_range_entry (FILE *file, struct range_entry *r);
void debug_range_entry (struct range_entry *r);
void init_range_entry (struct range_entry *r, tree exp, gimple *stmt);
bool no_side_effect_bb (basic_block bb);

#endif  /* GCC_SSA_REASSOC_H  */

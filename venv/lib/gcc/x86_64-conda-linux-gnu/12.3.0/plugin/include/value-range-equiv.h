/* Support routines for value ranges with equivalences.
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

#ifndef GCC_VALUE_RANGE_EQUIV_H
#define GCC_VALUE_RANGE_EQUIV_H

#include "value-range.h"

/* Note value_range_equiv cannot currently be used with GC memory,
   only value_range is fully set up for this.  */
class GTY((user)) value_range_equiv : public value_range
{
 public:
  value_range_equiv () : value_range () { m_equiv = NULL; }
  value_range_equiv (const value_range &);
  /* Deep-copies equiv bitmap argument.  */
  value_range_equiv (tree, tree, bitmap = NULL, value_range_kind = VR_RANGE);

  /* Shallow-copies equiv bitmap.  */
  value_range_equiv (const value_range_equiv &) /* = delete */;
  /* Shallow-copies equiv bitmap.  */
  value_range_equiv& operator=(const value_range_equiv &) /* = delete */;

  /* Move equiv bitmap from source range.  */
  void move (value_range_equiv *);

  /* Leaves equiv bitmap alone.  */
  void update (tree, tree, value_range_kind = VR_RANGE);
  /* Deep-copies equiv bitmap argument.  */
  void set (tree, tree, bitmap = NULL, value_range_kind = VR_RANGE);
  void set (tree);

  bool operator== (const value_range_equiv &) const /* = delete */;
  bool operator!= (const value_range_equiv &) const /* = delete */;
  void intersect (const value_range_equiv *);
  void union_ (const value_range_equiv *);
  bool equal_p (const value_range_equiv &, bool ignore_equivs) const;

  /* Types of value ranges.  */
  void set_undefined ();
  void set_varying (tree);

  /* Equivalence bitmap methods.  */
  bitmap equiv () const { return m_equiv; }
  void equiv_clear ();
  void equiv_add (const_tree, const value_range_equiv *,
		  bitmap_obstack * = NULL);

  /* Misc methods.  */
  void deep_copy (const value_range_equiv *);
  void dump (FILE *) const;
  void dump () const;

 private:
  /* Deep-copies bitmap argument.  */
  void set_equiv (bitmap);
  void check ();

  /* Set of SSA names whose value ranges are equivalent to this one.
     This set is only valid when TYPE is VR_RANGE or VR_ANTI_RANGE.  */
  bitmap m_equiv;
};

extern void dump_value_range (FILE *, const value_range_equiv *);

#endif // GCC_VALUE_RANGE_EQUIV_H

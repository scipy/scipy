/* Representation of thunks inside symbol table.
   Copyright (C) 2003-2022 Free Software Foundation, Inc.
   Contributed by Jan Hubicka

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

#ifndef GCC_SYMTAB_THUNKS_H
#define GCC_SYMTAB_THUNKS_H

/* This symbol annotation holds information about thunk.

   Thunks are basically wrappers around methods which are introduced in case
   of multiple inheritance in order to adjust the value of the "this" pointer
   or of the returned value.

   In the case of this-adjusting thunks, each back-end can override the
   can_output_mi_thunk/output_mi_thunk target hooks to generate a minimal thunk
   (with a tail call for instance) directly as assembly.  For the default hook
   or for the case where the can_output_mi_thunk hooks return false, the thunk
   is gimplified and lowered using the regular machinery.  */

struct GTY(()) thunk_info {
  /* Constructor.  */
  thunk_info ()
    : fixed_offset (0),
      virtual_value (0),
      indirect_offset (0),
      alias (NULL),
      this_adjusting (false),
      virtual_offset_p (false)
  {
  }
  /* Copy constructor.  */
  thunk_info (const thunk_info &t)
    : fixed_offset (t.fixed_offset),
      virtual_value (t.virtual_value),
      indirect_offset (t.indirect_offset),
      alias (t.alias),
      this_adjusting (t.this_adjusting),
      virtual_offset_p (t.virtual_offset_p)
  {
  }

  /* Compare for equiality.  */
  bool
  operator==(const thunk_info &other) const
  {
    return fixed_offset == other.fixed_offset
	   && virtual_value == other.virtual_value
	   && indirect_offset == other.indirect_offset
	   && this_adjusting == other.this_adjusting
	   && virtual_offset_p == other.virtual_offset_p;
  }
  bool
  operator!=(const thunk_info &other) const
  {
    return !(*this == other);
  }
  /* Copy operator.  */
  thunk_info &
  operator=(const thunk_info &other)
  {
    fixed_offset = other.fixed_offset;
    virtual_value = other.virtual_value;
    indirect_offset = other.indirect_offset;
    alias = other.alias;
    this_adjusting = other.this_adjusting;
    virtual_offset_p = other.virtual_offset_p;
    return *this;
  }

  /* Offset used to adjust "this".  */
  HOST_WIDE_INT fixed_offset;

  /* Offset in the virtual table to get the offset to adjust "this".  Valid iff
     VIRTUAL_OFFSET_P is true.  */
  HOST_WIDE_INT virtual_value;

  /* Offset from "this" to get the offset to adjust "this".  Zero means: this
     offset is to be ignored.  */
  HOST_WIDE_INT indirect_offset;

  /* Thunk target, i.e. the method that this thunk wraps.  Depending on the
     TARGET_USE_LOCAL_THUNK_ALIAS_P macro, this may have to be a new alias.  */
  tree alias;

  /* Nonzero for a "this" adjusting thunk and zero for a result adjusting
     thunk.  */
  bool this_adjusting;

  /* If true, this thunk is what we call a virtual thunk.  In this case:
     * for this-adjusting thunks, after the FIXED_OFFSET based adjustment is
       done, add to the result the offset found in the vtable at:
	 vptr + VIRTUAL_VALUE
     * for result-adjusting thunks, the FIXED_OFFSET adjustment is done after
       the virtual one.  */
  bool virtual_offset_p;



  /* Dump thunk_info.  */
  void dump (FILE *);

  /* Stream out thunk_info.  */
  void stream_out (class lto_simple_output_block *);

  /* Stream in trunk_info.  */
  void stream_in (class lto_input_block *);

  hashval_t hash ();



  /* Return thunk_info, if available.  */
  static thunk_info *get (cgraph_node *node);

  /* Return thunk_info possibly creating new one.  */
  static thunk_info *get_create (cgraph_node *node);

  /* Remove thunk_info.  */
  static void remove (cgraph_node *node);

  /* Add unprocessed thunk.  */
  void register_early (cgraph_node *node);

  /* Attach recorded thunks to cgraph_nodes.  */
  static void process_early_thunks ();

  /* Release all thunk_infos.  */
  static void release (void);
};

bool expand_thunk (cgraph_node *, bool, bool);

/* Return thunk_info, if available.  */
inline thunk_info *
thunk_info::get (cgraph_node *node)
{
  if (!symtab->m_thunks)
    return NULL;
  return symtab->m_thunks->get (node);
}

/* Remove thunk_info association for NODE.  */
inline void
thunk_info::remove (cgraph_node *node)
{
  symtab->m_thunks->remove (node);
}

/* Free thunk info summaries.  */
inline void
thunk_info::release ()
{
  if (symtab->m_thunks)
    ggc_delete (symtab->m_thunks);
  symtab->m_thunks = NULL;
}
#endif  /* GCC_SYMTAB_THUNKS_H  */

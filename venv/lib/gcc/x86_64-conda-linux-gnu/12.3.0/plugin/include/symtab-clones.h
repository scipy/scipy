/* Representation of adjustment made to virtual clones in the symbol table.
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

#ifndef GCC_SYMTAB_CLONES_H
#define GCC_SYMTAB_CLONES_H

struct GTY(()) clone_info
{
  /* Constructor.  */
  clone_info ()
    : tree_map (NULL),
      param_adjustments (NULL)
  {
  }
  /* Constants discovered by IPA-CP, i.e. which parameter should be replaced
     with what.  */
  vec<ipa_replace_map *, va_gc> *tree_map;
  /* Parameter modification that IPA-SRA decided to perform.  */
  ipa_param_adjustments *param_adjustments;

  /* Return clone_info, if available.  */
  static clone_info *get (cgraph_node *node);

  /* Return clone_info possibly creating new one.  */
  static clone_info *get_create (cgraph_node *node);

  /* Remove clone_info.  */
  static void remove (cgraph_node *node);

  /* Release all clone_infos.  */
  static void release (void);
};

/* Return clone_info, if available.  */
inline clone_info *
clone_info::get (cgraph_node *node)
{
  if (!symtab->m_clones)
    return NULL;
  return symtab->m_clones->get (node);
}


/* Remove clone_info association for NODE.  */
inline void
clone_info::remove (cgraph_node *node)
{
  symtab->m_clones->remove (node);
}

/* Free clone info summaries.  */
inline void
clone_info::release ()
{
  if (symtab->m_clones)
    ggc_delete (symtab->m_clones);
  symtab->m_clones = NULL;
}

#endif  /* GCC_SYMTAB_CLONES_H  */

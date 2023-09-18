/* Header file for Nested function decomposition for GIMPLE.
   Copyright (C) 2013-2022 Free Software Foundation, Inc.

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

#ifndef GCC_TREE_NESTED_H
#define GCC_TREE_NESTED_H

extern tree build_addr (tree);
extern void insert_field_into_struct (tree, tree);
extern void lower_nested_functions (tree);

class nested_function_info
{
public:
  /* Constructor.  */
  nested_function_info ()
    : origin (NULL),
      nested (NULL),
      next_nested (NULL)
  {
  }
  /* Copy constructor.  We can not simply copy the structure,
     because the linked lists would go wrong.  However we should never
     need that.  */
  nested_function_info (const nested_function_info &)
  {
     gcc_unreachable ();
  }
  ~nested_function_info ();

  /* Return nested_function_info, if available.  */
  static nested_function_info *get (cgraph_node *node);

  /* Return nested_function_info possibly creating new one.  */
  static nested_function_info *get_create (cgraph_node *node);

  /* Release all nested_function_infos.  */
  static void release (void);

  /* For nested functions points to function the node is nested in.  */
  cgraph_node *origin;
  /* Points to first nested function, if any.  */
  cgraph_node *nested;
  /* Pointer to the next function with same origin, if any.  */
  cgraph_node *next_nested;
};

extern void maybe_record_nested_function (cgraph_node *node);
extern void unnest_function (cgraph_node *node);

/* If there are functions nested in NODE, return first one.  */
inline cgraph_node *
first_nested_function (cgraph_node *node)
{
  nested_function_info *info = nested_function_info::get (node);
  return info ? info->nested : NULL;
}

/* Return next nested function (used to iterate from first_nested_function).  */
inline cgraph_node *
next_nested_function (cgraph_node *node)
{
  return nested_function_info::get (node)->next_nested;
}

/* Return origin of nested function (and NULL otherwise).  */
inline cgraph_node *
nested_function_origin (cgraph_node *node)
{
  nested_function_info *info = nested_function_info::get (node);
  return info ? info->origin : NULL;
}

#endif /* GCC_TREE_NESTED_H */

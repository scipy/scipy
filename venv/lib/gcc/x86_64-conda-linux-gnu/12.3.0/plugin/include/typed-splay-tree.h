/* A typesafe wrapper around libiberty's splay-tree.h.
   Copyright (C) 2015-2022 Free Software Foundation, Inc.

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

#ifndef GCC_TYPED_SPLAY_TREE_H
#define GCC_TYPED_SPLAY_TREE_H

/* Typesafe wrapper around libiberty's splay-tree.h.  */
template <typename KEY_TYPE, typename VALUE_TYPE>
class typed_splay_tree
{
 public:
  typedef KEY_TYPE key_type;
  typedef VALUE_TYPE value_type;

  typedef int (*compare_fn) (key_type, key_type);
  typedef void (*delete_key_fn) (key_type);
  typedef void (*delete_value_fn) (value_type);
  typedef int (*foreach_fn) (key_type, value_type, void *);

  typed_splay_tree (compare_fn,
		    delete_key_fn,
		    delete_value_fn);
  ~typed_splay_tree ();

  value_type lookup (key_type k);
  value_type predecessor (key_type k);
  value_type successor (key_type k);
  void insert (key_type k, value_type v);
  void remove (key_type k);
  value_type max ();
  value_type min ();
  int foreach (foreach_fn, void *);

 private:
  /* Copy and assignment ops are not supported.  */
  typed_splay_tree (const typed_splay_tree &);
  typed_splay_tree & operator = (const typed_splay_tree &);

  typedef key_type splay_tree_key;
  typedef value_type splay_tree_value;

  /* The nodes in the splay tree.  */
  struct splay_tree_node_s {
    /* The key.  */
    splay_tree_key key;

    /* The value.  */
    splay_tree_value value;

    /* The left and right children, respectively.  */
    splay_tree_node_s *left, *right;

    /* Used as temporary value for tree traversals.  */
    splay_tree_node_s *back;
  };
  typedef splay_tree_node_s *splay_tree_node;

  inline void KDEL (splay_tree_key);
  inline void VDEL (splay_tree_value);
  void splay_tree_delete_helper (splay_tree_node);
  static inline void rotate_left (splay_tree_node *,
				  splay_tree_node, splay_tree_node);
  static inline void rotate_right (splay_tree_node *,
				   splay_tree_node, splay_tree_node);
  void splay_tree_splay (splay_tree_key);
  static int splay_tree_foreach_helper (splay_tree_node,
					foreach_fn, void*);
  splay_tree_node splay_tree_insert (splay_tree_key, splay_tree_value);
  void splay_tree_remove (splay_tree_key key);
  splay_tree_node splay_tree_lookup (splay_tree_key key);
  splay_tree_node splay_tree_predecessor (splay_tree_key);
  splay_tree_node splay_tree_successor (splay_tree_key);
  splay_tree_node splay_tree_max ();
  splay_tree_node splay_tree_min ();

  static value_type node_to_value (splay_tree_node node);

  /* The root of the tree.  */
  splay_tree_node root;

  /* The comparision function.  */
  compare_fn comp;

  /* The deallocate-key function.  NULL if no cleanup is necessary.  */
  delete_key_fn delete_key;

  /* The deallocate-value function.  NULL if no cleanup is necessary.  */
  delete_value_fn delete_value;
};

/* Constructor for typed_splay_tree <K, V>.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline typed_splay_tree<KEY_TYPE, VALUE_TYPE>::
  typed_splay_tree (compare_fn compare_fn,
		    delete_key_fn delete_key_fn,
		    delete_value_fn delete_value_fn)
{
  root = NULL;
  comp = compare_fn;
  delete_key = delete_key_fn;
  delete_value = delete_value_fn;
}

/* Destructor for typed_splay_tree <K, V>.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline typed_splay_tree<KEY_TYPE, VALUE_TYPE>::
  ~typed_splay_tree ()
{
  splay_tree_delete_helper (root);
}

/* Lookup KEY, returning a value if present, and NULL
   otherwise.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline VALUE_TYPE
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::lookup (key_type key)
{
  splay_tree_node node = splay_tree_lookup (key);
  return node_to_value (node);
}

/* Return the immediate predecessor of KEY, or NULL if there is no
   predecessor.  KEY need not be present in the tree.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline VALUE_TYPE
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::predecessor (key_type key)
{
  splay_tree_node node = splay_tree_predecessor (key);
  return node_to_value (node);
}

/* Return the immediate successor of KEY, or NULL if there is no
   successor.  KEY need not be present in the tree.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline VALUE_TYPE
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::successor (key_type key)
{
  splay_tree_node node = splay_tree_successor (key);
  return node_to_value (node);
}

/* Insert a new node (associating KEY with VALUE).  If a
   previous node with the indicated KEY exists, its data is replaced
   with the new value.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline void
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::insert (key_type key,
						value_type value)
{
  splay_tree_insert (key, value);
}

/* Remove a node (associating KEY with VALUE).  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline void
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::remove (key_type key)
{
  splay_tree_remove (key);
}

/* Get the value with maximal key.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline VALUE_TYPE
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::max ()
{
  return node_to_value (splay_tree_max ());
}

/* Get the value with minimal key.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline VALUE_TYPE
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::min ()
{
  return node_to_value (splay_tree_min ());
}

/* Call OUTER_CB, passing it the OUTER_USER_DATA, for every node,
   following an in-order traversal.  If OUTER_CB ever returns a non-zero
   value, the iteration ceases immediately, and the value is returned.
   Otherwise, this function returns 0.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline int
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::foreach (foreach_fn foreach_fn,
						 void *user_data)
{
  return splay_tree_foreach_helper (root, foreach_fn, user_data);
}

/* Internal function for converting from splay_tree_node to
   VALUE_TYPE.  */
template <typename KEY_TYPE, typename VALUE_TYPE>
inline VALUE_TYPE
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::node_to_value (splay_tree_node node)
{
  if (node)
    return node->value;
  else
    return 0;
}

template <typename KEY_TYPE, typename VALUE_TYPE>
inline void
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::KDEL(splay_tree_key x)
{
  if (delete_key)
    (*delete_key)(x);
}

template <typename KEY_TYPE, typename VALUE_TYPE>
inline void
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::VDEL(splay_tree_value x)
{
  if (delete_value)
    (*delete_value)(x);
}

/* Deallocate NODE (a member of SP), and all its sub-trees.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
void
typed_splay_tree<KEY_TYPE,
		 VALUE_TYPE>::splay_tree_delete_helper (splay_tree_node node)
{
  splay_tree_node pending = NULL;
  splay_tree_node active = NULL;

  if (!node)
    return;

  KDEL (node->key);
  VDEL (node->value);

  /* We use the "back" field to hold the "next" pointer.  */
  node->back = pending;
  pending = node;

  /* Now, keep processing the pending list until there aren't any
     more.  This is a little more complicated than just recursing, but
     it doesn't toast the stack for large trees.  */

  while (pending)
    {
      active = pending;
      pending = NULL;
      while (active)
	{
	  splay_tree_node temp;

	  /* active points to a node which has its key and value
	     deallocated, we just need to process left and right.  */

	  if (active->left)
	    {
	      KDEL (active->left->key);
	      VDEL (active->left->value);
	      active->left->back = pending;
	      pending = active->left;
	    }
	  if (active->right)
	    {
	      KDEL (active->right->key);
	      VDEL (active->right->value);
	      active->right->back = pending;
	      pending = active->right;
	    }

	  temp = active;
	  active = temp->back;
	  delete temp;
	}
    }
}

/* Rotate the edge joining the left child N with its parent P.  PP is the
   grandparents' pointer to P.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline void
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::rotate_left (splay_tree_node *pp,
						     splay_tree_node p,
						     splay_tree_node n)
{
  splay_tree_node tmp;
  tmp = n->right;
  n->right = p;
  p->left = tmp;
  *pp = n;
}

/* Rotate the edge joining the right child N with its parent P.  PP is the
   grandparents' pointer to P.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
inline void
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::rotate_right (splay_tree_node *pp,
						      splay_tree_node p,
						      splay_tree_node n)
{
  splay_tree_node tmp;
  tmp = n->left;
  n->left = p;
  p->right = tmp;
  *pp = n;
}

/* Bottom up splay of key.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
void
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_splay (splay_tree_key key)
{
  if (root == NULL)
    return;

  do {
    int cmp1, cmp2;
    splay_tree_node n, c;

    n = root;
    cmp1 = (*comp) (key, n->key);

    /* Found.  */
    if (cmp1 == 0)
      return;

    /* Left or right?  If no child, then we're done.  */
    if (cmp1 < 0)
      c = n->left;
    else
      c = n->right;
    if (!c)
      return;

    /* Next one left or right?  If found or no child, we're done
       after one rotation.  */
    cmp2 = (*comp) (key, c->key);
    if (cmp2 == 0
	|| (cmp2 < 0 && !c->left)
	|| (cmp2 > 0 && !c->right))
      {
	if (cmp1 < 0)
	  rotate_left (&root, n, c);
	else
	  rotate_right (&root, n, c);
	return;
      }

    /* Now we have the four cases of double-rotation.  */
    if (cmp1 < 0 && cmp2 < 0)
      {
	rotate_left (&n->left, c, c->left);
	rotate_left (&root, n, n->left);
      }
    else if (cmp1 > 0 && cmp2 > 0)
      {
	rotate_right (&n->right, c, c->right);
	rotate_right (&root, n, n->right);
      }
    else if (cmp1 < 0 && cmp2 > 0)
      {
	rotate_right (&n->left, c, c->right);
	rotate_left (&root, n, n->left);
      }
    else if (cmp1 > 0 && cmp2 < 0)
      {
	rotate_left (&n->right, c, c->left);
	rotate_right (&root, n, n->right);
      }
  } while (1);
}

/* Call FN, passing it the DATA, for every node below NODE, all of
   which are from SP, following an in-order traversal.  If FN every
   returns a non-zero value, the iteration ceases immediately, and the
   value is returned.  Otherwise, this function returns 0.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
int
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_foreach_helper (
						splay_tree_node node,
						foreach_fn fn, void *data)
{
  int val;
  splay_tree_node stack;

  /* A non-recursive implementation is used to avoid filling the stack
     for large trees.  Splay trees are worst case O(n) in the depth of
     the tree.  */

  stack = NULL;
  val = 0;

  for (;;)
    {
      while (node != NULL)
	{
	  node->back = stack;
	  stack = node;
	  node = node->left;
	}

      if (stack == NULL)
	break;

      node = stack;
      stack = stack->back;

      val = (*fn) (node->key, node->value, data);
      if (val)
	break;

      node = node->right;
    }

  return val;
}

/* Insert a new node (associating KEY with DATA) into SP.  If a
   previous node with the indicated KEY exists, its data is replaced
   with the new value.  Returns the new node.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
typename typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_node
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_insert (
						splay_tree_key key,
						splay_tree_value value)
{
  int comparison = 0;

  splay_tree_splay (key);

  if (root)
    comparison = (*comp)(root->key, key);

  if (root && comparison == 0)
    {
      /* If the root of the tree already has the indicated KEY, just
	 replace the value with VALUE.  */
      VDEL(root->value);
      root->value = value;
    }
  else
    {
      /* Create a new node, and insert it at the root.  */
      splay_tree_node node;

      node = new splay_tree_node_s;
      node->key = key;
      node->value = value;

      if (!root)
	node->left = node->right = 0;
      else if (comparison < 0)
	{
	  node->left = root;
	  node->right = node->left->right;
	  node->left->right = 0;
	}
      else
	{
	  node->right = root;
	  node->left = node->right->left;
	  node->right->left = 0;
	}

      root = node;
    }

  return root;
}

/* Remove KEY from SP.  It is not an error if it did not exist.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
void
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_remove (splay_tree_key key)
{
  splay_tree_splay (key);

  if (root && (*comp) (root->key, key) == 0)
    {
      splay_tree_node left, right;

      left = root->left;
      right = root->right;

      /* Delete the root node itself.  */
      VDEL (root->value);
      delete root;

      /* One of the children is now the root.  Doesn't matter much
	 which, so long as we preserve the properties of the tree.  */
      if (left)
	{
	  root = left;

	  /* If there was a right child as well, hang it off the
	     right-most leaf of the left child.  */
	  if (right)
	    {
	      while (left->right)
		left = left->right;
	      left->right = right;
	    }
	}
      else
	root = right;
    }
}

/* Lookup KEY in SP, returning VALUE if present, and NULL
   otherwise.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
typename typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_node
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_lookup (splay_tree_key key)
{
  splay_tree_splay (key);

  if (root && (*comp)(root->key, key) == 0)
    return root;
  else
    return 0;
}

/* Return the node in SP with the greatest key.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
typename typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_node
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_max ()
{
  splay_tree_node n = root;

  if (!n)
    return NULL;

  while (n->right)
    n = n->right;

  return n;
}

/* Return the node in SP with the smallest key.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
typename typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_node
typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_min ()
{
  splay_tree_node n = root;

  if (!n)
    return NULL;

  while (n->left)
    n = n->left;

  return n;
}

/* Return the immediate predecessor KEY, or NULL if there is no
   predecessor.  KEY need not be present in the tree.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
typename typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_node
typed_splay_tree<KEY_TYPE,
		 VALUE_TYPE>::splay_tree_predecessor (splay_tree_key key)
{
  int comparison;
  splay_tree_node node;

  /* If the tree is empty, there is certainly no predecessor.  */
  if (!root)
    return NULL;

  /* Splay the tree around KEY.  That will leave either the KEY
     itself, its predecessor, or its successor at the root.  */
  splay_tree_splay (key);
  comparison = (*comp)(root->key, key);

  /* If the predecessor is at the root, just return it.  */
  if (comparison < 0)
    return root;

  /* Otherwise, find the rightmost element of the left subtree.  */
  node = root->left;
  if (node)
    while (node->right)
      node = node->right;

  return node;
}

/* Return the immediate successor KEY, or NULL if there is no
   successor.  KEY need not be present in the tree.  */

template <typename KEY_TYPE, typename VALUE_TYPE>
typename typed_splay_tree<KEY_TYPE, VALUE_TYPE>::splay_tree_node
typed_splay_tree<KEY_TYPE,
		 VALUE_TYPE>::splay_tree_successor (splay_tree_key key)
{
  int comparison;
  splay_tree_node node;

  /* If the tree is empty, there is certainly no successor.  */
  if (!root)
    return NULL;

  /* Splay the tree around KEY.  That will leave either the KEY
     itself, its predecessor, or its successor at the root.  */
  splay_tree_splay (key);
  comparison = (*comp)(root->key, key);

  /* If the successor is at the root, just return it.  */
  if (comparison > 0)
    return root;

  /* Otherwise, find the leftmost element of the right subtree.  */
  node = root->right;
  if (node)
    while (node->left)
      node = node->left;

  return node;
}

#endif  /* GCC_TYPED_SPLAY_TREE_H  */

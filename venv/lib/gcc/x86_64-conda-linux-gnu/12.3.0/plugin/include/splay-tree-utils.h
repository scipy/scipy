// Splay tree utilities                                             -*- C++ -*-
// Copyright (C) 2020-2022 Free Software Foundation, Inc.
//
// This file is part of GCC.
//
// GCC is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 3, or (at your option) any later
// version.
//
// GCC is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// along with GCC; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

// Implement splay tree node accessors for a class that stores its
// two child nodes in a member variable of the form:
//
//    Node m_children[2];
template<typename Node>
class default_splay_tree_accessors
{
public:
  using node_type = Node;

  static auto
  child (node_type node, unsigned int index)
    -> decltype (node->m_children[index]) &
  {
    return node->m_children[index];
  }
};

// Implement splay tree node accessors for a class that stores its
// two child nodes in a member variable of the form:
//
//    Node m_children[2];
//
// and also stores its parent node in a member variable of the form:
//
//    Node m_parent;
template<typename Node>
class default_splay_tree_accessors_with_parent
  : public default_splay_tree_accessors<Node>
{
public:
  using node_type = Node;

  static auto
  parent (node_type node) -> decltype (node->m_parent) &
  {
    return node->m_parent;
  }
};

// Base is a splay tree accessor class for nodes that have no parent field.
// Base therefore provides a Base::child method but does not provide a
// Base::parent method.  Extend Base with dummy routines for setting the
// parent, which is a no-op when the parent is not stored.
template<typename Base>
class splay_tree_accessors_without_parent : public Base
{
public:
  using typename Base::node_type;

  static void set_parent (node_type, node_type) {}
};

// Base is splay tree accessor class for nodes that have a parent field.
// Base therefore provides both Base::child and Base::parent methods.
// Extend Base with routines for setting the parent.
template<typename Base>
class splay_tree_accessors_with_parent : public Base
{
public:
  using typename Base::node_type;

  // Record that NODE's parent is now NEW_PARENT.
  static void
  set_parent (node_type node, node_type new_parent)
  {
    Base::parent (node) = new_parent;
  }
};

// A base class that provides some splay tree operations that are common
// to both rooted_splay_tree and rootless_splay_tree.
//
// Nodes in the splay tree have type Accessors::node_type; this is
// usually a pointer type.  The Accessors class provides the following
// static member functions for accessing nodes:
//
// - Accessors::child (NODE, INDEX)
//     INDEX is guaranteed to be 0 or 1.  If INDEX is 0, return a reference
//     to where NODE's left child is stored, otherwise return a reference
//     to where NODE's right child is stored.
//
// - Accessors::set_parent (NODE, PARENT)
//     Record that NODE's parent node is now PARENT.
template<typename Accessors>
class base_splay_tree : protected Accessors
{
public:
  using typename Accessors::node_type;

  // INDEX is either 0 or 1.  If INDEX is 0, insert CHILD immediately
  // before NODE, otherwise insert CHILD immediately after NODE.
  //
  // Complexity: O(1).
  static void insert_child (node_type node, unsigned int index,
			    node_type child);

  // Print NODE and its child nodes to PP for debugging purposes,
  // using PRINTER (PP, N) to print the data for node N.
  template<typename Printer>
  static void print (pretty_printer *pp, node_type node, Printer printer);

protected:
  using Accessors::set_parent;

  static node_type get_child (node_type, unsigned int);
  static void set_child (node_type, unsigned int, node_type);
  static node_type promote_child (node_type, unsigned int);
  static void promote_child (node_type, unsigned int, node_type);

  template<unsigned int N>
  static node_type splay_limit (node_type);

  static node_type remove_node_internal (node_type);

  template<typename Printer>
  static void print (pretty_printer *pp, node_type node, Printer printer,
		     char, vec<char> &);
};

// This class provides splay tree routines for cases in which the root
// of the splay tree is known.  It works with both nodes that store
// their parent node and nodes that don't.
//
// The class is lightweight: it only contains a single root node.
template<typename Accessors>
class rooted_splay_tree : public base_splay_tree<Accessors>
{
  using parent = base_splay_tree<Accessors>;

public:
  using typename Accessors::node_type;

protected:
  // The root of the splay tree, or node_type () if the tree is empty.
  node_type m_root;

public:
  rooted_splay_tree () : m_root () {}

  // Construct a tree with the specified root node.
  rooted_splay_tree (node_type root) : m_root (root) {}

  // Return the root of the tree.
  node_type root () const { return m_root; }

  // Return true if the tree contains any nodes.
  explicit operator bool () const { return m_root; }

  // Dereference the root node.
  node_type operator-> () { return m_root; }

  // Insert NEW_NODE into the splay tree, if no equivalent node already
  // exists.  For a given node N, COMPARE (N) should return:
  //
  // - a negative value if NEW_NODE should come before N
  // - zero if NEW_NODE and N are the same
  // - a positive value if NEW_NODE should come after N
  //
  // Return true if NEW_NODE was inserted.
  //
  // On return, NEW_NODE or its equivalent is the root of the tree.
  //
  // Complexity: amortized O(C log N), worst-cast O(C N), where C is
  // the complexity of the comparison.
  template<typename Comparator>
  bool insert (node_type new_node, Comparator compare);

  // Insert NEW_NODE into the splay tree, given that NEW_NODE is the
  // maximum node of the new tree.  On return, NEW_NODE is also the
  // root of the tree.
  //
  // Complexity: O(1).
  void insert_max_node (node_type new_node);

  // Splice NEXT_TREE onto this one, given that all nodes in NEXT_TREE
  // are greater than the maximum node in this tree.  NEXT_TREE should
  // not be used afterwards.
  //
  // Complexity: O(1) if the root of the splay tree is already the maximum
  // node.  Otherwise amortized O(log N), worst-cast O(N).
  void splice_next_tree (rooted_splay_tree next_tree);

  // The root of the tree is currently the maximum node.  Replace it
  // with NEW_NODE.
  //
  // Complexity: O(1).
  void replace_max_node_at_root (node_type new_node);

  // Remove the root node of the splay tree.
  //
  // Complexity: O(1) if removing the maximum or minimum node.
  // Otherwise amortized O(log N), worst-cast O(N).
  void remove_root ();

  // Split the left child of the current root out into a separate tree
  // and return the new tree.
  rooted_splay_tree split_before_root ();

  // Split the right child of the current root out into a separate tree
  // and return the new tree.
  rooted_splay_tree split_after_root ();

  // If the root is not the minimum node of the splay tree, bring the previous
  // node to the root and return true, otherwise return false.
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  bool splay_prev_node ();

  // If the root is not the maximum node of the splay tree, bring the next
  // node to the root and return true, otherwise return false.
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  bool splay_next_node ();

  // Bring the minimum node of the splay tree to the root.
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  void splay_min_node ();

  // Bring the maximum node of the splay tree to the root.
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  void splay_max_node ();

  // Return the minimum node of the splay tree, or node_type () if the
  // tree is empty.  On return, the minimum node (if any) is also the
  // root of the tree.
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  node_type min_node ();

  // Return the maximum node of the splay tree, or node_type () if the
  // tree is empty.  On return, the maximum node (if any) is also the
  // root of the tree.
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  node_type max_node ();

  // Search the splay tree.  For a given node N, COMPARE (N) should return:
  //
  // - a negative value if N is bigger than the node being searched for
  // - zero if N is the node being searched for
  // - a positive value if N is smaller than the node being searched for
  //
  // If the node that COMPARE is looking for exists, install it as the root
  // node of the splay tree.  Otherwise, arbitrarily pick either:
  //
  // - the maximum node that is smaller than the node being searched for or
  // - the minimum node that is bigger than the node being searched for
  //
  // and install that node as the root instead.
  //
  // Return the result of COMPARE for the new root.
  //
  // This form of lookup is intended for cases in which both the following
  // are true:
  //
  // (a) The work that COMPARE needs to do to detect if a node is too big
  //     is the same as the work that COMPARE needs to do to detect if a
  //     node is too small.  (This is not true of range comparisons,
  //     for example.)
  //
  // (b) COMPARE is (or might be) relatively complex.
  //
  // This form of lookup is also useful if the items being compared naturally
  // provide a <=>-style comparison result, without the result having to be
  // forced by the equivalent of a ?: expression.
  //
  // The implementation only invokes COMPARE once per node.
  //
  // Complexity: amortized O(C log N), worst-cast O(C N), where C is
  // the complexity of the comparison.
  template<typename Comparator>
  auto lookup (Comparator compare) -> decltype (compare (m_root));

  // Search the splay tree.  For a given node N, WANT_SOMETHING_SMALLER (N)
  // is true if N is too big and WANT_SOMETHING_BIGGER (N) is true if N
  // is too small.  Both functions return false if N is the node being
  // searched for.
  //
  // If the node that is being searched for exists, install it as the root
  // node of the splay tree and return 0.  Otherwise, arbitrarily choose
  // between these two options:
  //
  // - Install the maximum node that is smaller than the node being
  //   searched for as the root of the splay tree and return 1.
  //
  // - Install the minimum node that is bigger than the node being
  //   searched for and return -1.
  //
  // This form of lookup is intended for cases in which either of the
  // following are true:
  //
  // (a) WANT_SOMETHING_SMALLER and WANT_SOMETHING_BIGGER test different
  //     parts of the node's data.  For example, when comparing ranges,
  //     WANT_SOMETHING_SMALLER would test the lower limit of the given
  //     node's range while WANT_SOMETHING_BIGGER would test the upper
  //     limit of the given node's range.
  //
  // (b) There is no significant overhead to calling both
  //     WANT_SOMETHING_SMALLER and WANT_SOMETHING_BIGGER for the same node.
  //
  // Complexity: amortized O(C log N), worst-cast O(C N), where C is
  // the complexity of the comparisons.
  template<typename LeftPredicate, typename RightPredicate>
  int lookup (LeftPredicate want_something_smaller,
	      RightPredicate want_something_bigger);

  // Keep the ability to print subtrees.
  using parent::print;

  // Print the tree to PP for debugging purposes, using PRINTER (PP, N)
  // to print the data for node N.
  template<typename Printer>
  void print (pretty_printer *pp, Printer printer) const;

protected:
  using parent::get_child;
  using parent::set_child;
  using parent::promote_child;

  using parent::set_parent;

  template<unsigned int N>
  bool splay_neighbor ();
};

// Provide splay tree routines for nodes of type Accessors::node_type,
// which doesn't have a parent field.  Use Accessors::child to access
// the children of a node.
template<typename Accessors>
using splay_tree_without_parent
  = rooted_splay_tree<splay_tree_accessors_without_parent<Accessors>>;

// A splay tree for nodes of type Node, which is usually a pointer type.
// The child nodes are stored in a member variable:
//
//    Node m_children[2];
//
// Node does not have a parent field.
template<typename Node>
using default_splay_tree
  = splay_tree_without_parent<default_splay_tree_accessors<Node>>;

// A simple splay tree node that stores a value of type T.
template<typename T>
class splay_tree_node
{
  friend class default_splay_tree_accessors<splay_tree_node *>;

public:
  splay_tree_node () = default;
  splay_tree_node (T value) : m_value (value), m_children () {}

  T &value () { return m_value; }
  const T &value () const { return m_value; }

private:
  T m_value;
  splay_tree_node *m_children[2];
};

// A splay tree whose nodes hold values of type T.
template<typename T>
using splay_tree = default_splay_tree<splay_tree_node<T> *>;

// Provide splay tree routines for cases in which the root of the tree
// is not explicitly stored.
//
// The nodes of the tree have type Accessors::node_type, which is usually
// a pointer type.  The nodes have a link back to their parent.
//
// The Accessors class provides the following static member functions:
//
// - Accessors::child (NODE, INDEX)
//     INDEX is guaranteed to be 0 or 1.  If INDEX is 0, return a reference
//     to where NODE's left child is stored, otherwise return a reference
//     to where NODE's right child is stored.
//
// - Accessors::parent (NODE)
//     Return a reference to where NODE's parent is stored.
template<typename Accessors>
class rootless_splay_tree
  : public base_splay_tree<splay_tree_accessors_with_parent<Accessors>>
{
  using full_accessors = splay_tree_accessors_with_parent<Accessors>;
  using parent = base_splay_tree<full_accessors>;

public:
  using rooted = rooted_splay_tree<full_accessors>;

  using typename Accessors::node_type;

  // Remove NODE from the splay tree.  Return the node that replaces it,
  // or null if NODE had no children.
  //
  // Complexity: O(1) if removing the maximum or minimum node.
  // Otherwise amortized O(log N), worst-cast O(N).
  static node_type remove_node (node_type node);

  // Splay NODE so that it becomes the root of the splay tree.
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  static void splay (node_type node);

  // Like splay, but take advantage of the fact that NODE is known to be
  // the minimum node in the tree.
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  static void splay_known_min_node (node_type node);

  // Like splay, but take advantage of the fact that NODE is known to be
  // the maximum node in the tree.
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  static void splay_known_max_node (node_type node);

  // Splay NODE while looking for an ancestor node N for which PREDICATE (N)
  // is true.  If such an ancestor node exists, stop the splay operation
  // early and return PREDICATE (N).  Otherwise, complete the splay operation
  // and return DEFAULT_RESULT.  In the latter case, NODE is now the root of
  // the splay tree.
  //
  // Note that this routine only examines nodes that happen to be ancestors
  // of NODE.  It does not search the full tree.
  //
  // Complexity: amortized O(P log N), worst-cast O(P N), where P is the
  // complexity of the predicate.
  template<typename DefaultResult, typename Predicate>
  static auto splay_and_search (node_type node, DefaultResult default_result,
				Predicate predicate)
    -> decltype (predicate (node, 0));

  // NODE1 and NODE2 are known to belong to the same splay tree.  Return:
  //
  // -1 if NODE1 < NODE2
  // 0 if NODE1 == NODE2
  // 1 if NODE1 > NODE2
  //
  // Complexity: amortized O(log N), worst-cast O(N).
  static int compare_nodes (node_type node1, node_type node2);

protected:
  using parent::get_child;
  using parent::set_child;
  using parent::promote_child;

  static node_type get_parent (node_type);
  using parent::set_parent;

  static unsigned int child_index (node_type, node_type);

  static int compare_nodes_one_way (node_type, node_type);

  template<unsigned int N>
  static void splay_known_limit (node_type);
};

// Provide rootless splay tree routines for nodes of type Node.
// The child nodes are stored in a member variable:
//
//    Node m_children[2];
//
// and the parent node is stored in a member variable:
//
//    Node m_parent;
template<typename Node>
using default_rootless_splay_tree
  = rootless_splay_tree<default_splay_tree_accessors_with_parent<Node>>;

#include "splay-tree-utils.tcc"

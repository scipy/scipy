/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SPLAY_H_
#define HIGHS_SPLAY_H_

#include <cassert>

#include "util/HighsType.h"

/// top down splay operation to maintain a binary search tree. The search tree
/// is assumed to be stored in an array/vector and therefore uses integers
/// instead of pointers to link nodes. The splay operation on a given key will
/// return the new root node index. The new root node will be the node with the
/// given key, if it exists in the tree, and otherwise it is the node traversed
/// last in a binary tree search for the key.
///
/// The GetLeft/GetRight lambdas receive an index of a node and must
/// return references to an integer holding the left and the right child node
/// indices. The GetKey lambda must return a type that is comparable
/// to KeyT.
template <typename KeyT, typename GetLeft, typename GetRight, typename GetKey>
HighsInt highs_splay(const KeyT& key, HighsInt root, GetLeft&& get_left,
                     GetRight&& get_right, GetKey&& get_key) {
  if (root == -1) return -1;

  HighsInt Nleft = -1;
  HighsInt Nright = -1;
  HighsInt* lright = &Nright;
  HighsInt* rleft = &Nleft;

  while (true) {
    if (key < get_key(root)) {
      HighsInt left = get_left(root);
      if (left == -1) break;
      if (key < get_key(left)) {
        HighsInt y = left;
        get_left(root) = get_right(y);
        get_right(y) = root;
        root = y;
        if (get_left(root) == -1) break;
      }

      *rleft = root;
      rleft = &get_left(root);
      root = get_left(root);
    } else if (key > get_key(root)) {
      HighsInt right = get_right(root);
      if (right == -1) break;
      if (key > get_key(right)) {
        HighsInt y = right;
        get_right(root) = get_left(y);
        get_left(y) = root;
        root = y;
        if (get_right(root) == -1) break;
      }

      *lright = root;
      lright = &get_right(root);
      root = get_right(root);
    } else
      break;
  }

  *lright = get_left(root);
  *rleft = get_right(root);
  get_left(root) = Nright;
  get_right(root) = Nleft;

  return root;
}

/// links a new node into the binary tree rooted at the given reference to the
/// root node. Lambdas must behave as described in highs_splay above.
/// Equal keys are put to the right subtree.
template <typename GetLeft, typename GetRight, typename GetKey>
void highs_splay_link(HighsInt linknode, HighsInt& root, GetLeft&& get_left,
                      GetRight&& get_right, GetKey&& get_key) {
  if (root == -1) {
    get_left(linknode) = -1;
    get_right(linknode) = -1;
    root = linknode;
    return;
  }

  root = highs_splay(get_key(linknode), root, get_left, get_right, get_key);

  if (get_key(linknode) < get_key(root)) {
    get_left(linknode) = get_left(root);
    get_right(linknode) = root;
    get_left(root) = -1;
  } else {
    assert(get_key(linknode) > get_key(root));
    get_right(linknode) = get_right(root);
    get_left(linknode) = root;
    get_right(root) = -1;
  }

  root = linknode;
}

/// unlinks a new node into the binary tree rooted at the given reference to the
/// root node. Lambdas must behave as described in highs_splay above.
template <typename GetLeft, typename GetRight, typename GetKey>
void highs_splay_unlink(HighsInt unlinknode, HighsInt& root, GetLeft&& get_left,
                        GetRight&& get_right, GetKey&& get_key) {
  assert(root != -1);
  root = highs_splay(get_key(unlinknode), root, get_left, get_right, get_key);
  assert(get_key(root) == get_key(unlinknode));

  // in case keys can be equal it might happen that we did not splay the correct
  // node to the top since equal keys are put to the right subtree, we recurse
  // into that part of the tree
  if (root != unlinknode) {
    highs_splay_unlink(unlinknode, get_right(root), get_left, get_right,
                       get_key);
    return;
  }

  assert(root == unlinknode);

  if (get_left(unlinknode) == -1) {
    root = get_right(unlinknode);
  } else {
    root = highs_splay(get_key(unlinknode), get_left(unlinknode), get_left,
                       get_right, get_key);
    get_right(root) = get_right(unlinknode);
  }
}

#endif

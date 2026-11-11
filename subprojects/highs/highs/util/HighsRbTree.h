/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_RBTREE_H_
#define HIGHS_RBTREE_H_

#include <algorithm>
#include <cassert>
#include <type_traits>

#include "util/HighsType.h"

namespace highs {

template <typename T>
struct RbTreeLinks {
  enum Direction {
    kLeft = 0,
    kRight = 1,
  };

  T child[2];

  using ParentStorageType =
      typename std::make_unsigned<typename std::conditional<
          std::is_pointer<T>::value, uintptr_t, T>::type>::type;

  static constexpr int colorBitPos() {
    return std::is_pointer<T>::value ? 0 : sizeof(ParentStorageType) * 8 - 1;
  }

  static constexpr ParentStorageType colorBitMask() {
    return ParentStorageType{1} << colorBitPos();
  }

  ParentStorageType parentAndColor;

  template <typename U = T,
            typename std::enable_if<std::is_integral<U>::value, int>::type = 0>
  static constexpr T noLink() {
    return -1;
  }

  template <typename U = T,
            typename std::enable_if<std::is_pointer<U>::value, int>::type = 0>
  static constexpr T noLink() {
    return nullptr;
  }

  bool getColor() const { return (parentAndColor >> colorBitPos()); }

  bool isBlack() const { return getColor() == 0; }

  bool isRed() const { return getColor() == 1; }

  void makeRed() { parentAndColor |= colorBitMask(); }

  void makeBlack() { parentAndColor &= ~colorBitMask(); }

  void setColor(bool color) {
    makeBlack();
    parentAndColor |= ParentStorageType(color) << colorBitPos();
  }

  void setParent(T p) {
    if (std::is_pointer<T>::value)
      parentAndColor =
          (colorBitMask() & parentAndColor) | (ParentStorageType)(p);
    else
      parentAndColor =
          (colorBitMask() & parentAndColor) | (ParentStorageType)(p + 1);
  }

  T getParent() const {
    if (std::is_pointer<T>::value)
      return (T)(parentAndColor & ~colorBitMask());
    else
      return ((T)(parentAndColor & ~colorBitMask())) - 1;
  }
};

template <typename Impl>
struct RbTreeTraits;

template <typename Impl>
class RbTree {
  enum Dir {
    kLeft = 0,
    kRight = 1,
  };

  using KeyType = typename RbTreeTraits<Impl>::KeyType;
  using LinkType = typename RbTreeTraits<Impl>::LinkType;

  LinkType& rootNode;
  // static_assert(
  //     std::is_same<RbTreeLinks&, decltype(static_cast<Impl*>(nullptr)
  //                                             ->getRbTreeLinks(0))>::value,
  //     "RbTree implementation must provide a getRbTreeLinks() function that "
  //     "returns a non-const reference of RbTreeLinks given the index of a
  //     node");

  // static_assert(std::is_same<const RbTreeLinks&,
  //                           decltype(static_cast<const Impl*>(nullptr)
  //                                        ->getRbTreeLinks(0))>::value,
  //              "RbTree implementation must provide a getRbTreeLinks() const "
  //              "function that returns a const reference of RbTreeLinks given
  //              " "the index of a node");
  // static_assert(
  //    std::is_same<bool, decltype((*static_cast<KeyType*>(nullptr)) <
  //                                (*static_cast<KeyType*>(nullptr)))>::value,
  //    "RbTree implementation must provide a getKey() function that, given the
  //    " "index of a node, returns its key which must have a type that has "
  //    "operator< defined");

  static constexpr LinkType kNoLink = RbTreeLinks<LinkType>::noLink();

  bool isRed(LinkType node) const {
    return node != kNoLink &&
           static_cast<const Impl*>(this)->getRbTreeLinks(node).isRed();
  }

  bool isBlack(LinkType node) const {
    return node == kNoLink ||
           static_cast<const Impl*>(this)->getRbTreeLinks(node).isBlack();
  }

  void makeRed(LinkType node) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).makeRed();
  }

  void makeBlack(LinkType node) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).makeBlack();
  }

  void setColor(LinkType node, HighsUInt color) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).setColor(color != 0);
  }

  HighsUInt getColor(LinkType node) const {
    return static_cast<const Impl*>(this)->getRbTreeLinks(node).getColor();
  }

  void setParent(LinkType node, LinkType parent) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).setParent(parent);
  }

  LinkType getParent(LinkType node) const {
    return static_cast<const Impl*>(this)->getRbTreeLinks(node).getParent();
  }

  KeyType getKey(LinkType node) const {
    return static_cast<const Impl*>(this)->getKey(node);
  }

  static constexpr Dir opposite(Dir dir) { return Dir(1 - dir); }

  LinkType getChild(LinkType node, Dir dir) const {
    return static_cast<const Impl*>(this)->getRbTreeLinks(node).child[dir];
  }

  void setChild(LinkType node, Dir dir, LinkType child) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).child[dir] = child;
  }

  void rotate(LinkType x, Dir dir) {
    LinkType y = getChild(x, opposite(dir));
    LinkType yDir = getChild(y, dir);
    setChild(x, opposite(dir), yDir);
    if (yDir != kNoLink) setParent(yDir, x);

    LinkType pX = getParent(x);
    setParent(y, pX);

    if (pX == kNoLink)
      rootNode = y;
    else
      setChild(pX, Dir((x != getChild(pX, dir)) ^ dir), y);

    setChild(y, dir, x);
    setParent(x, y);
  }

  void insertFixup(LinkType z) {
    LinkType pZ = getParent(z);
    while (isRed(pZ)) {
      LinkType zGrandParent = getParent(pZ);
      assert(zGrandParent != kNoLink);

      Dir dir = Dir(getChild(zGrandParent, kLeft) == pZ);

      LinkType y = getChild(zGrandParent, dir);
      if (isRed(y)) {
        makeBlack(pZ);
        makeBlack(y);
        makeRed(zGrandParent);
        z = zGrandParent;
      } else {
        if (z == getChild(pZ, dir)) {
          z = pZ;
          rotate(z, opposite(dir));
          pZ = getParent(z);
          zGrandParent = getParent(pZ);
          assert(zGrandParent != kNoLink);
        }

        makeBlack(pZ);
        makeRed(zGrandParent);
        rotate(zGrandParent, dir);
      }

      pZ = getParent(z);
    }

    makeBlack(rootNode);
  }

  void transplant(LinkType u, LinkType v, LinkType& nilParent) {
    LinkType p = getParent(u);

    if (p == kNoLink)
      rootNode = v;
    else
      setChild(p, Dir(u != getChild(p, kLeft)), v);

    if (v == kNoLink)
      nilParent = p;
    else
      setParent(v, p);
  }

  void deleteFixup(LinkType x, const LinkType nilParent) {
    while (x != rootNode && isBlack(x)) {
      Dir dir;

      LinkType p = x == kNoLink ? nilParent : getParent(x);
      dir = Dir(x == getChild(p, kLeft));
      LinkType w = getChild(p, dir);
      assert(w != kNoLink);

      if (isRed(w)) {
        makeBlack(w);
        makeRed(p);
        rotate(p, opposite(dir));
        assert((x == kNoLink && p == nilParent) ||
               (x != kNoLink && p == getParent(x)));
        w = getChild(p, dir);
        assert(w != kNoLink);
      }

      if (isBlack(getChild(w, kLeft)) && isBlack(getChild(w, kRight))) {
        makeRed(w);
        x = p;
      } else {
        if (isBlack(getChild(w, dir))) {
          makeBlack(getChild(w, opposite(dir)));
          makeRed(w);
          rotate(w, dir);
          assert((x == kNoLink && p == nilParent) ||
                 (x != kNoLink && p == getParent(x)));
          w = getChild(p, dir);
        }
        setColor(w, getColor(p));
        makeBlack(p);
        makeBlack(getChild(w, dir));
        rotate(p, opposite(dir));
        x = rootNode;
      }
    }

    if (x != kNoLink) makeBlack(x);
  }

 public:
  RbTree(LinkType& rootNode) : rootNode(rootNode) {}

  bool empty() const { return rootNode == kNoLink; }

  LinkType first(LinkType x) const {
    if (x == kNoLink) return kNoLink;

    while (true) {
      LinkType lX = getChild(x, kLeft);
      if (lX == kNoLink) return x;
      x = lX;
    }
  }

  LinkType last(LinkType x) const {
    if (x == kNoLink) return kNoLink;

    while (true) {
      LinkType rX = getChild(x, kRight);
      if (rX == kNoLink) return x;
      x = rX;
    }
  }

  LinkType first() const { return first(rootNode); }

  LinkType last() const { return last(rootNode); }

  LinkType successor(LinkType x) const {
    LinkType y = getChild(x, kRight);
    if (y != kNoLink) return first(y);

    y = getParent(x);
    while (y != kNoLink && x == getChild(y, kRight)) {
      x = y;
      y = getParent(x);
    }

    return y;
  }

  LinkType predecessor(LinkType x) const {
    LinkType y = getChild(x, kLeft);
    if (y != kNoLink) return last(y);

    y = getParent(x);
    while (y != kNoLink && x == getChild(y, kLeft)) {
      x = y;
      y = getParent(x);
    }

    return y;
  }

  std::pair<LinkType, bool> find(const KeyType& key, LinkType treeRoot) const {
    LinkType y = kNoLink;
    LinkType x = treeRoot;
    while (x != kNoLink) {
      HighsInt cmp = 1 - (getKey(x) < key) + (key < getKey(x));
      switch (cmp) {
        case 0:
          y = x;
          x = getChild(y, kRight);
          break;
        case 1:
          return std::make_pair(x, true);
        case 2:
          y = x;
          x = getChild(y, kLeft);
      }
    }

    return std::make_pair(y, false);
  }

  std::pair<LinkType, bool> find(const KeyType& key) const {
    return find(key, rootNode);
  }

  void link(LinkType z, LinkType parent) {
    setParent(z, parent);
    if (parent == kNoLink)
      rootNode = z;
    else
      setChild(parent, Dir(getKey(parent) < getKey(z)), z);

    setChild(z, kLeft, kNoLink);
    setChild(z, kRight, kNoLink);
    makeRed(z);
    insertFixup(z);
  }

  void link(LinkType z) {
    LinkType y = kNoLink;
    LinkType x = rootNode;
    while (x != kNoLink) {
      y = x;
      x = getChild(y, Dir(getKey(x) < getKey(z)));
    }

    static_cast<Impl*>(this)->link(z, y);
  }

  void unlink(LinkType z) {
    LinkType nilParent = kNoLink;
    LinkType y = z;
    bool yWasBlack = isBlack(y);
    LinkType x;

    if (getChild(z, kLeft) == kNoLink) {
      x = getChild(z, kRight);
      transplant(z, x, nilParent);
    } else if (getChild(z, kRight) == kNoLink) {
      x = getChild(z, kLeft);
      transplant(z, x, nilParent);
    } else {
      y = first(getChild(z, kRight));
      yWasBlack = isBlack(y);
      x = getChild(y, kRight);
      if (getParent(y) == z) {
        if (x == kNoLink)
          nilParent = y;
        else
          setParent(x, y);
      } else {
        transplant(y, getChild(y, kRight), nilParent);
        LinkType zRight = getChild(z, kRight);
        setChild(y, kRight, zRight);
        setParent(zRight, y);
      }
      transplant(z, y, nilParent);
      LinkType zLeft = getChild(z, kLeft);
      setChild(y, kLeft, zLeft);
      setParent(zLeft, y);
      setColor(y, getColor(z));
    }

    if (yWasBlack) deleteFixup(x, nilParent);
  }
};

template <typename Impl>
class CacheMinRbTree : public RbTree<Impl> {
  using LinkType = typename RbTreeTraits<Impl>::LinkType;
  LinkType& first_;

 public:
  CacheMinRbTree(LinkType& rootNode, LinkType& first)
      : RbTree<Impl>(rootNode), first_(first) {}

  LinkType first() const { return first_; };
  using RbTree<Impl>::first;

  void link(LinkType z, LinkType parent) {
    if (first_ == parent) {
      if (parent == RbTreeLinks<LinkType>::noLink() ||
          static_cast<const Impl*>(this)->getKey(z) <
              static_cast<const Impl*>(this)->getKey(parent))
        first_ = z;
    }

    RbTree<Impl>::link(z, parent);
  }
  using RbTree<Impl>::link;

  void unlink(LinkType z) {
    if (z == first_) first_ = this->successor(first_);
    RbTree<Impl>::unlink(z);
  }
};

}  // namespace highs

#endif

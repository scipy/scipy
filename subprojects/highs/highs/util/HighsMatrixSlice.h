/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsMatrixSlice.h
 * @brief Provides a uniform interface to iterate rows and columns in different
 * underlying matrix storage formats
 */

#ifndef UTIL_HIGHS_MATRIX_SLICE_H_
#define UTIL_HIGHS_MATRIX_SLICE_H_

#include <cstddef>
#include <iterator>
#include <vector>

#include "util/HighsType.h"

template <typename StorageFormat>
class HighsMatrixSlice;

struct HighsEmptySlice;
struct HighsCompressedSlice;
struct HighsIndexedSlice;
struct HighsTripletListSlice;
struct HighsTripletTreeSliceInOrder;
struct HighsTripletTreeSlicePreOrder;
struct HighsTripletPositionSlice;

class HighsSliceNonzero {
  template <typename>
  friend class HighsMatrixSlice;

 private:
  const HighsInt* index_;
  const double* value_;

 public:
  HighsSliceNonzero() = default;
  HighsSliceNonzero(const HighsInt* index, const double* value)
      : index_(index), value_(value) {}
  HighsInt index() const { return *index_; }
  double value() const { return *value_; }
};

template <>
class HighsMatrixSlice<HighsEmptySlice> {
 public:
  using iterator = const HighsSliceNonzero*;
  static constexpr const HighsSliceNonzero* begin() { return nullptr; }
  static constexpr const HighsSliceNonzero* end() { return nullptr; }
};

template <>
class HighsMatrixSlice<HighsCompressedSlice> {
  const HighsInt* index;
  const double* value;
  HighsInt len;

 public:
  class iterator {
    HighsSliceNonzero pos_;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(const HighsInt* index, const double* value) : pos_(index, value) {}
    iterator() = default;

    iterator operator++(int) {
      iterator prev = *this;
      ++pos_.index_;
      ++pos_.value_;
      return prev;
    }
    iterator& operator++() {
      ++pos_.index_;
      ++pos_.value_;
      return *this;
    }
    reference operator*() const { return pos_; }
    pointer operator->() const { return &pos_; }
    iterator operator+(difference_type v) const {
      iterator i = *this;
      i.pos_.index_ += v;
      i.pos_.value_ += v;
      return i;
    }

    bool operator==(const iterator& rhs) const {
      return pos_.index_ == rhs.pos_.index_;
    }
    bool operator!=(const iterator& rhs) const {
      return pos_.index_ != rhs.pos_.index_;
    }
  };

  HighsMatrixSlice(const HighsInt* index_, const double* value_, HighsInt len_)
      : index(index_), value(value_), len(len_) {}
  iterator begin() const { return iterator{index, value}; }
  iterator end() const { return iterator{index + len, nullptr}; }
};

template <>
class HighsMatrixSlice<HighsIndexedSlice> {
  const HighsInt* index;
  const double* denseValues;
  HighsInt len;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const double* denseValues;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(const HighsInt* index_, const double* denseValues_)
        : pos_(index_, denseValues_), denseValues(denseValues_) {}
    iterator() = default;

    iterator operator++(int) {
      iterator prev = *this;
      ++pos_.index_;
      return prev;
    }
    iterator& operator++() {
      ++pos_.index_;
      return *this;
    }
    reference operator*() {
      pos_.value_ = &denseValues[*pos_.index_];
      return pos_;
    }
    pointer operator->() {
      pos_.value_ = &denseValues[*pos_.index_];
      return &pos_;
    }
    iterator operator+(difference_type v) const {
      iterator i = *this;

      while (v > 0) {
        --v;
        ++i;
      }

      return i;
    }

    bool operator==(const iterator& rhs) const {
      return pos_.index_ == rhs.pos_.index_;
    }
    bool operator!=(const iterator& rhs) const {
      return pos_.index_ != rhs.pos_.index_;
    }
  };

  HighsMatrixSlice(const HighsInt* index, const double* denseValues,
                   HighsInt len)
      : index(index), denseValues(denseValues), len(len) {}
  iterator begin() const { return iterator{index, denseValues}; }
  iterator end() const { return iterator{index + len, nullptr}; }
};

template <>
class HighsMatrixSlice<HighsTripletListSlice> {
  const HighsInt* nodeIndex;
  const double* nodeValue;
  const HighsInt* nodeNext;
  HighsInt head;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const HighsInt* nodeNext;
    HighsInt currentNode;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(HighsInt node) : pos_(), nodeNext(nullptr), currentNode(node) {}
    iterator(const HighsInt* nodeIndex, const double* nodeValue,
             const HighsInt* nodeNext, HighsInt node)
        : pos_(node == -1 ? nullptr : nodeIndex + node,
               node == -1 ? nullptr : nodeValue + node),
          nodeNext(nodeNext),
          currentNode(node) {}
    iterator() = default;

    iterator& operator++() {
      pos_.index_ -= currentNode;
      pos_.value_ -= currentNode;
      currentNode = nodeNext[currentNode];
      pos_.index_ += currentNode;
      pos_.value_ += currentNode;
      return *this;
    }
    iterator operator++(int) {
      iterator prev = *this;
      ++(*this);
      return prev;
    }
    reference operator*() { return pos_; }
    pointer operator->() { return &pos_; }
    iterator operator+(difference_type v) const {
      iterator i = *this;

      while (v > 0) {
        --v;
        ++i;
      }

      return i;
    }

    HighsInt position() const { return currentNode; }

    bool operator==(const iterator& rhs) const {
      return currentNode == rhs.currentNode;
    }
    bool operator!=(const iterator& rhs) const {
      return currentNode != rhs.currentNode;
    }
  };

  HighsMatrixSlice(const HighsInt* nodeIndex, const double* nodeValue,
                   const HighsInt* nodeNext, HighsInt head)
      : nodeIndex(nodeIndex),
        nodeValue(nodeValue),
        nodeNext(nodeNext),
        head(head) {}
  iterator begin() const {
    return iterator{nodeIndex, nodeValue, nodeNext, head};
  }
  iterator end() const { return iterator{-1}; }
};

template <>
class HighsMatrixSlice<HighsTripletTreeSlicePreOrder> {
  const HighsInt* nodeIndex;
  const double* nodeValue;
  const HighsInt* nodeLeft;
  const HighsInt* nodeRight;
  HighsInt root;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const HighsInt* nodeLeft;
    const HighsInt* nodeRight;
    std::vector<HighsInt> stack;
    HighsInt currentNode;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(HighsInt node)
        : pos_(), nodeLeft(nullptr), nodeRight(nullptr), currentNode(node) {}
    iterator(const HighsInt* nodeIndex, const double* nodeValue,
             const HighsInt* nodeLeft, const HighsInt* nodeRight, HighsInt node)
        : pos_(nodeIndex + node, nodeValue + node),
          nodeLeft(nodeLeft),
          nodeRight(nodeRight),
          currentNode(node) {
      stack.reserve(16);
      stack.push_back(-1);
    }
    iterator() = default;

    iterator& operator++() {
      HighsInt offset = -currentNode;
      if (nodeLeft[currentNode] != -1) {
        if (nodeRight[currentNode] != -1)
          stack.push_back(nodeRight[currentNode]);
        currentNode = nodeLeft[currentNode];
      } else if (nodeRight[currentNode] != -1) {
        currentNode = nodeRight[currentNode];
      } else {
        currentNode = stack.back();
        stack.pop_back();
      }
      offset += currentNode;
      pos_.index_ += offset;
      pos_.value_ += offset;
      return *this;
    }

    iterator operator++(int) {
      iterator prev = *this;
      ++(*this);
      return prev;
    }
    reference operator*() { return pos_; }
    pointer operator->() { return &pos_; }
    iterator operator+(difference_type v) const {
      iterator i = *this;

      while (v > 0) {
        --v;
        ++i;
      }

      return i;
    }

    HighsInt position() const { return currentNode; }

    bool operator==(const iterator& rhs) const {
      return currentNode == rhs.currentNode;
    }
    bool operator!=(const iterator& rhs) const {
      return currentNode != rhs.currentNode;
    }
  };

  HighsMatrixSlice(const HighsInt* nodeIndex, const double* nodeValue,
                   const HighsInt* nodeLeft, const HighsInt* nodeRight,
                   HighsInt root)
      : nodeIndex(nodeIndex),
        nodeValue(nodeValue),
        nodeLeft(nodeLeft),
        nodeRight(nodeRight),
        root(root) {}

  iterator end() const { return iterator{-1}; }
  iterator begin() const {
    // avoid allocation if there are no elements
    if (root == -1) return end();
    return iterator{nodeIndex, nodeValue, nodeLeft, nodeRight, root};
  }
};

template <>
class HighsMatrixSlice<HighsTripletTreeSliceInOrder> {
  const HighsInt* nodeIndex;
  const double* nodeValue;
  const HighsInt* nodeLeft;
  const HighsInt* nodeRight;
  HighsInt root;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const HighsInt* nodeLeft;
    const HighsInt* nodeRight;
    std::vector<HighsInt> stack;
    HighsInt currentNode;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(HighsInt node)
        : pos_(), nodeLeft(nullptr), nodeRight(nullptr), currentNode(node) {}
    iterator(const HighsInt* nodeIndex, const double* nodeValue,
             const HighsInt* nodeLeft, const HighsInt* nodeRight, HighsInt node)
        : pos_(nodeIndex, nodeValue),
          nodeLeft(nodeLeft),
          nodeRight(nodeRight),
          currentNode(node) {
      stack.reserve(16);
      stack.push_back(-1);
      if (currentNode == -1) return;
      while (nodeLeft[currentNode] != -1) {
        stack.push_back(currentNode);
        currentNode = nodeLeft[currentNode];
      }

      pos_.index_ += currentNode;
      pos_.value_ += currentNode;
    }
    iterator() = default;

    iterator& operator++() {
      HighsInt offset = -currentNode;
      if (nodeRight[currentNode] != -1) {
        currentNode = nodeRight[currentNode];
        while (nodeLeft[currentNode] != -1) {
          stack.push_back(currentNode);
          currentNode = nodeLeft[currentNode];
        }
      } else {
        currentNode = stack.back();
        stack.pop_back();
      }
      offset += currentNode;
      pos_.index_ += offset;
      pos_.value_ += offset;
      return *this;
    }

    iterator operator++(int) {
      iterator prev = *this;
      ++(*this);
      return prev;
    }
    reference operator*() { return pos_; }
    pointer operator->() { return &pos_; }
    iterator operator+(difference_type v) const {
      iterator i = *this;

      while (v > 0) {
        --v;
        ++i;
      }

      return i;
    }

    HighsInt position() const { return currentNode; }

    bool operator==(const iterator& rhs) const {
      return currentNode == rhs.currentNode;
    }
    bool operator!=(const iterator& rhs) const {
      return currentNode != rhs.currentNode;
    }
  };

  HighsMatrixSlice(const HighsInt* nodeIndex, const double* nodeValue,
                   const HighsInt* nodeLeft, const HighsInt* nodeRight,
                   HighsInt root)
      : nodeIndex(nodeIndex),
        nodeValue(nodeValue),
        nodeLeft(nodeLeft),
        nodeRight(nodeRight),
        root(root) {}

  iterator end() const { return iterator{-1}; }
  iterator begin() const {
    // avoid allocation if there are no elements
    if (root == -1) return end();
    return iterator{nodeIndex, nodeValue, nodeLeft, nodeRight, root};
  }
};

template <>
class HighsMatrixSlice<HighsTripletPositionSlice> {
  const HighsInt* nodeIndex;
  const double* nodeValue;
  const HighsInt* nodePositions;
  HighsInt len;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const HighsInt* node;
    HighsInt currentNode;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(const HighsInt* node) : pos_(), node(node), currentNode(0) {}
    iterator(const HighsInt* nodeIndex, const double* nodeValue,
             const HighsInt* node)
        : pos_(nodeIndex, nodeValue), node(node), currentNode(0) {}
    iterator() = default;

    iterator& operator++() {
      ++node;
      return *this;
    }

    iterator operator++(int) {
      iterator prev = *this;
      ++(*this);
      return prev;
    }
    reference operator*() {
      HighsInt offset = -currentNode + *node;
      currentNode = *node;
      pos_.index_ += offset;
      pos_.value_ += offset;
      return pos_;
    }
    pointer operator->() {
      HighsInt offset = -currentNode + *node;
      currentNode = *node;
      pos_.index_ += offset;
      pos_.value_ += offset;
      return &pos_;
    }
    iterator operator+(difference_type v) const {
      iterator i = *this;
      i.node += v;
      return i;
    }

    HighsInt position() const { return currentNode; }

    bool operator==(const iterator& rhs) const { return node == rhs.node; }

    bool operator!=(const iterator& rhs) const { return node != rhs.node; }
  };

  HighsMatrixSlice(const HighsInt* nodeIndex, const double* nodeValue,
                   const HighsInt* nodePositions, HighsInt len)
      : nodeIndex(nodeIndex),
        nodeValue(nodeValue),
        nodePositions(nodePositions),
        len(len) {}

  iterator begin() const {
    return iterator{nodeIndex, nodeValue, nodePositions};
  }

  iterator end() const { return iterator{nodePositions + len}; }
};

struct HighsEmptySlice : public HighsMatrixSlice<HighsEmptySlice> {
  using HighsMatrixSlice<HighsEmptySlice>::HighsMatrixSlice;
};
struct HighsCompressedSlice : public HighsMatrixSlice<HighsCompressedSlice> {
  using HighsMatrixSlice<HighsCompressedSlice>::HighsMatrixSlice;
};
struct HighsIndexedSlice : public HighsMatrixSlice<HighsIndexedSlice> {
  using HighsMatrixSlice<HighsIndexedSlice>::HighsMatrixSlice;
};
struct HighsTripletListSlice : public HighsMatrixSlice<HighsTripletListSlice> {
  using HighsMatrixSlice<HighsTripletListSlice>::HighsMatrixSlice;
};
struct HighsTripletTreeSliceInOrder
    : public HighsMatrixSlice<HighsTripletTreeSliceInOrder> {
  using HighsMatrixSlice<HighsTripletTreeSliceInOrder>::HighsMatrixSlice;
};
struct HighsTripletTreeSlicePreOrder
    : public HighsMatrixSlice<HighsTripletTreeSlicePreOrder> {
  using HighsMatrixSlice<HighsTripletTreeSlicePreOrder>::HighsMatrixSlice;
};
struct HighsTripletPositionSlice
    : public HighsMatrixSlice<HighsTripletPositionSlice> {
  using HighsMatrixSlice<HighsTripletPositionSlice>::HighsMatrixSlice;
};

#endif

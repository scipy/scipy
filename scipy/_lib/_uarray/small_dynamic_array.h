#pragma once
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <new>
#include <type_traits>


/** Fixed size dynamic array with small buffer optimisation */
template <typename T, ptrdiff_t SmallCapacity = 1>
class SmallDynamicArray {
  ptrdiff_t size_;
  union {
    T elements[SmallCapacity];
    T * array;
  } storage_;

  bool is_small() const { return size_ <= SmallCapacity; }

  void destroy_buffer(T * first, T * last) noexcept {
    for (; first < last; ++first) {
      first->~T();
    }
  }

  void default_construct_buffer(T * first, T * last) noexcept(
      std::is_nothrow_destructible<T>::value) {
    auto cur = first;
    try {
      for (; cur < last; ++cur) {
        new (cur) T(); // Placement new
      }
    } catch (...) {
      // If construction failed, destroy already constructed values
      destroy_buffer(first, cur);
      throw;
    }
  }

  void move_construct_buffer(T * first, T * last, T * d_first) noexcept(
      std::is_nothrow_move_constructible<T>::value) {
    T * d_cur = d_first;

    try {
      for (; first < last; ++first, ++d_cur) {
        new (d_cur) T(std::move(*first)); // Placement new
      }
    } catch (...) {
      destroy_buffer(d_first, d_cur);
      throw;
    }
  }

  void allocate() {
    assert(size_ >= 0);
    if (is_small())
      return;

    storage_.array = (T *)malloc(size_ * sizeof(T));
    if (!storage_.array) {
      throw std::bad_alloc();
    }
  }

  void deallocate() noexcept {
    if (!is_small()) {
      free(storage_.array);
    }
  }

public:
  using value_type = T;
  using iterator = value_type *;
  using const_iterator = const value_type *;
  using reference = value_type &;
  using const_reference = const value_type &;
  using size_type = ptrdiff_t;

  SmallDynamicArray(): size_(0) {}

  explicit SmallDynamicArray(size_t size): size_(size) {
    allocate();
    auto first = begin();
    try {
      default_construct_buffer(first, first + size_);
    } catch (...) {
      deallocate();
      throw;
    }
  }

  SmallDynamicArray(size_type size, const T & fill_value): size_(size) {
    allocate();
    try {
      std::uninitialized_fill_n(begin(), size_, fill_value);
    } catch (...) {
      deallocate();
      throw;
    }
  }

  SmallDynamicArray(const SmallDynamicArray & copy): size_(copy.size_) {
    allocate();
    try {
      std::uninitialized_copy_n(copy.begin(), size_, begin());
    } catch (...) {
      deallocate();
      throw;
    }
  }

  SmallDynamicArray(SmallDynamicArray && move) noexcept(
      std::is_nothrow_move_constructible<T>::value)
      : size_(move.size_) {
    if (!move.is_small()) {
      storage_.array = move.storage_.array;
      move.storage_.array = nullptr;
      move.storage_.size = 0;
      return;
    }

    auto m_first = move.begin();
    try {
      move_construct_buffer(m_first, m_first + size_, begin());
    } catch (...) {
      destroy_buffer(m_first, m_first + size_);
      move.size_ = 0;
      size_ = 0;
      throw;
    }

    destroy_buffer(&move.storage_.elements[0], &move.storage_.elements[size_]);
    move.size_ = 0;
  }

  ~SmallDynamicArray() { clear(); }

  SmallDynamicArray & operator=(const SmallDynamicArray & copy) {
    if (&copy == this)
      return *this;

    clear();

    size_ = copy.size;
    try {
      allocate();
    } catch (...) {
      size_ = 0;
      throw;
    }

    try {
      std::uninitialized_copy_n(copy.begin(), size_, begin());
    } catch (...) {
      deallocate();
      size_ = 0;
      throw;
    }
    return *this;
  }

  SmallDynamicArray & operator=(SmallDynamicArray && move) noexcept(
      std::is_nothrow_move_constructible<T>::value) {
    if (&move == this)
      return *this;

    if (!move.is_small()) {
      storage_.array = move.storage_.array;
      size_ = move.size_;
      move.storage_.array = nullptr;
      move.size_ = 0;
      return *this;
    }

    clear();

    size_ = move.size_;
    auto m_first = move.begin();
    try {
      move_construct_buffer(m_first, m_first + size_, begin());
    } catch (...) {
      destroy_buffer(m_first, m_first + size_);
      move.size_ = 0;
      size_ = 0;
      throw;
    }

    destroy_buffer(m_first, m_first + size_);
    move.size_ = 0;
    return *this;
  }

  void clear() noexcept(std::is_nothrow_destructible<T>::value) {
    auto first = begin();
    destroy_buffer(first, first + size_);
    deallocate();
    size_ = 0;
  }

  iterator begin() {
    return is_small() ? &storage_.elements[0] : storage_.array;
  }

  const_iterator cbegin() const {
    return is_small() ? &storage_.elements[0] : storage_.array;
  }

  iterator end() { return begin() + size_; }

  const_iterator cend() const { return cbegin() + size_; }

  const_iterator begin() const { return cbegin(); }

  const_iterator end() const { return cend(); }

  size_type size() const { return size_; }

  const_reference operator[](size_type idx) const {
    assert(0 <= idx && idx < size_);
    return begin()[idx];
  }

  reference operator[](size_type idx) {
    assert(0 <= idx && idx < size_);
    return begin()[idx];
  }
};

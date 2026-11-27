/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsUndoStack.h
 * @brief Records changes to variables so that they can be undone
 */
#ifndef UTIL_HIGHSUNDOSTACK_H_
#define UTIL_HIGHSUNDOSTACK_H_

#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

// Records the old values of variables before they are modified, so that
// the changes can be undone, most recent first.
//
// Usage:
//  - set(x, value) records the current value of x, then assigns value to x
//  - save(x) only records the current value of x, which the caller then
//    modifies
//  - m = mark() returns a handle to the current point; resetTo(m) undoes
//    only the changes recorded since then
//  - reset() undoes all changes; the destructor calls it
//
// Restrictions:
//  - the variables recorded must outlive the HighsUndoStack instance, since
//    references to them are stored
//  - each recorded value is copied, so avoid recording large objects
//  - the types recorded must be copy constructible and copy assignable
class HighsUndoStack {
 public:
  HighsUndoStack() = default;
  ~HighsUndoStack() { reset(); }
  HighsUndoStack(const HighsUndoStack&) = delete;
  HighsUndoStack& operator=(const HighsUndoStack&) = delete;

  template <typename T, typename U>
  void set(T& field, U&& value) {
    save(field);
    field = std::forward<U>(value);
  }

  template <typename T>
  void save(T& field) {
    const T old = field;
    undo_.push_back([&field, old]() { field = old; });
  }

  size_t mark() const { return undo_.size(); }

  void resetTo(const size_t m) {
    while (undo_.size() > m) {
      undo_.back()();
      undo_.pop_back();
    }
  }

  void reset() { resetTo(0); }

 private:
  std::vector<std::function<void()>> undo_;
};

#endif

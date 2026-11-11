/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HSet.cpp
 * @brief
 */
#include "util/HSet.h"

#include <cassert>

bool HSet::setup(const HighsInt size, const HighsInt max_entry,
                 const bool output_flag, FILE* log_file, const bool debug,
                 const bool allow_assert) {
  setup_ = false;
  if (size <= 0) return false;
  if (max_entry < min_entry) return false;
  max_entry_ = max_entry;
  debug_ = debug;
  allow_assert_ = allow_assert;
  output_flag_ = output_flag;
  log_file_ = log_file;
  entry_.resize(size);
  pointer_.assign(max_entry_ + 1, no_pointer);
  count_ = 0;
  setup_ = true;
  return true;
}

void HSet::clear() {
  if (!setup_) setup(1, 0);
  pointer_.assign(max_entry_ + 1, no_pointer);
  count_ = 0;
  if (debug_) debug();
}

bool HSet::add(const HighsInt entry) {
  if (entry < min_entry) return false;
  if (!setup_) setup(1, entry);
  if (entry > max_entry_) {
    // Entry exceeds what's allowable so far so can't be in the list
    pointer_.resize(entry + 1);
    for (HighsInt ix = max_entry_ + 1; ix < entry; ix++)
      pointer_[ix] = no_pointer;
    max_entry_ = entry;
  } else if (pointer_[entry] > no_pointer) {
    // Duplicate
    if (debug_) debug();
    return false;
  }
  // New entry
  HighsInt size = entry_.size();
  if (count_ == size) {
    size++;
    entry_.resize(size);
  }
  pointer_[entry] = count_;
  entry_[count_++] = entry;
  if (debug_) debug();
  return true;
}

bool HSet::remove(const HighsInt entry) {
  if (!setup_) {
    setup(1, 0);
    if (debug_) debug();
    return false;
  }
  if (entry < min_entry) return false;
  if (entry > max_entry_) return false;
  HighsInt pointer = pointer_[entry];
  if (pointer == no_pointer) return false;
  pointer_[entry] = no_pointer;
  if (pointer < count_ - 1) {
    HighsInt last_entry = entry_[count_ - 1];
    entry_[pointer] = last_entry;
    pointer_[last_entry] = pointer;
  }
  count_--;
  if (debug_) debug();
  return true;
}

bool HSet::in(const HighsInt entry) const {
  if (entry < min_entry) return false;
  if (entry > max_entry_) return false;
  return pointer_[entry] != no_pointer;
}

bool HSet::debug() const {
  if (!setup_) {
    if (output_flag_) fprintf(log_file_, "HSet: ERROR setup_ not called\n");
    if (allow_assert_) assert(setup_);
    return false;
  }
  bool max_entry_ok = max_entry_ >= min_entry;
  if (!max_entry_ok) {
    if (output_flag_) {
      fprintf(log_file_,
              "HSet: ERROR max_entry_ = %" HIGHSINT_FORMAT
              " < %" HIGHSINT_FORMAT "\n",
              max_entry_, min_entry);
      print();
    }
    if (allow_assert_) assert(max_entry_ok);
    return false;
  }
  HighsInt size = entry_.size();
  bool size_count_ok = size >= count_;
  if (!size_count_ok) {
    if (output_flag_) {
      fprintf(log_file_,
              "HSet: ERROR entry_.size() = %" HIGHSINT_FORMAT
              " is less than count_ = %" HIGHSINT_FORMAT "\n",
              size, count_);
      print();
    }
    if (allow_assert_) assert(size_count_ok);
    return false;
  }
  // Check pointer_ is consistent with count_ and entry_
  HighsInt count = 0;
  for (HighsInt ix = 0; ix <= max_entry_; ix++) {
    HighsInt pointer = pointer_[ix];
    if (pointer == no_pointer) continue;
    bool pointer_ok = pointer >= 0 && pointer < count_;
    if (!pointer_ok) {
      if (output_flag_) {
        fprintf(log_file_,
                "HSet: ERROR pointer_[%" HIGHSINT_FORMAT "] = %" HIGHSINT_FORMAT
                " is not in [0, %" HIGHSINT_FORMAT "]\n",
                ix, pointer, count_);
        print();
      }
      if (allow_assert_) assert(pointer_ok);
      return false;
    }
    count++;
    HighsInt entry = entry_[pointer];
    bool entry_ok = entry == ix;
    if (!entry_ok) {
      if (output_flag_) {
        fprintf(log_file_,
                "HSet: ERROR entry_[%" HIGHSINT_FORMAT "] is %" HIGHSINT_FORMAT
                ", not %" HIGHSINT_FORMAT "\n",
                pointer, entry, ix);
        print();
      }
      if (allow_assert_) assert(entry_ok);
      return false;
    }
  }
  bool count_ok = count == count_;
  if (!count_ok) {
    if (output_flag_) {
      fprintf(log_file_,
              "HSet: ERROR pointer_ has %" HIGHSINT_FORMAT
              " pointers, not %" HIGHSINT_FORMAT "\n",
              count, count_);
      print();
    }
    if (allow_assert_) assert(count_ok);
    return false;
  }
  return true;
}

void HSet::print() const {
  if (!setup_) return;
  if (log_file_ == NULL) return;
  HighsInt size = entry_.size();
  fprintf(log_file_, "\nSet(%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT "):\n",
          size, max_entry_);
  fprintf(log_file_, "Pointers: Pointers|");
  for (HighsInt ix = 0; ix <= max_entry_; ix++) {
    if (pointer_[ix] != no_pointer)
      fprintf(log_file_, " %4" HIGHSINT_FORMAT "", pointer_[ix]);
  }
  fprintf(log_file_, "\n");
  fprintf(log_file_, "          Entries |");
  for (HighsInt ix = 0; ix <= max_entry_; ix++) {
    if (pointer_[ix] != no_pointer)
      fprintf(log_file_, " %4" HIGHSINT_FORMAT "", ix);
  }
  fprintf(log_file_, "\n");
  fprintf(log_file_, "Entries:  Indices |");
  for (HighsInt ix = 0; ix < count_; ix++)
    fprintf(log_file_, " %4" HIGHSINT_FORMAT "", ix);
  fprintf(log_file_, "\n");
  fprintf(log_file_, "          Entries |");
  for (HighsInt ix = 0; ix < count_; ix++)
    fprintf(log_file_, " %4" HIGHSINT_FORMAT "", entry_[ix]);
  fprintf(log_file_, "\n");
}

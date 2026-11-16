/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HSet.h
 * @brief Set structure for HiGHS.
 */

// Maintains an unordered set of distinct non-negative integer entries,
// allowing entries to be removed from the set at cost O(1)
#ifndef UTIL_HSET_H_
#define UTIL_HSET_H_

#include <cstdio>
#include <vector>

#include "util/HighsType.h"

// #include <iostream>

using std::vector;

const HighsInt min_entry = 0;
const HighsInt no_pointer = min_entry - 1;
/**
 * @brief Class for the set structure for HiGHS
 */
class HSet {
 public:
  /**
   * @brief Initialise a set. Neither limit is binding, but more
   * efficient memory-wise if known in advance
   */
  bool setup(const HighsInt size,  //!< Dimension of the set to be initialised
             const HighsInt max_entry,  //!< Maximum entry to be in the set.
             const bool output_flag = false,  //!< Option for output
             FILE* log_stream = NULL,         //!< File stream for output
             const bool debug = false,        //!< Debug mode
             const bool allow_assert = true   //!< Allow asserts in debug
  );

  /**
   * @brief Clear the set
   */
  void clear();
  /**
   * @brief Add entry to the set
   */
  bool add(const HighsInt entry);
  /**
   * @brief Remove entry from the set
   */
  bool remove(const HighsInt entry);
  /**
   * @brief Returns whether entry is in the set
   */
  bool in(const HighsInt entry) const;
  /**
   * @brief Returns the number of entries in the set
   */
  const HighsInt& count() const { return count_; }
  /**
   * @brief Returns the entries in the set
   */
  const vector<HighsInt>& entry() const { return entry_; }
  /**
   * @brief Print out the set and pointer entries not set to no_pointer
   */
  void print() const;
  /**
   * @brief Debug the set
   */
  bool debug() const;

 private:
  HighsInt count_ = 0;      //!< Number of entries
  vector<HighsInt> entry_;  //!< Entries
  bool setup_ = false;
  bool debug_ = false;
  bool allow_assert_ = true;
  bool output_flag_ = false;
  FILE* log_file_;
  HighsInt max_entry_;        //!< Maximum entry to be in the set.
  vector<HighsInt> pointer_;  //!< Set of pointers into the set
};
#endif /* UTIL_HSET_H_ */

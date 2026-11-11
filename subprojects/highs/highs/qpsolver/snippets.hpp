/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_SNIPPETS_HPP__
#define __SRC_LIB_SNIPPETS_HPP__

#include <algorithm>
#include <vector>

template <typename T>
bool contains(const std::vector<T>& vec, const T& element) {
  return std::find(vec.begin(), vec.end(), element) != vec.end();
}

template <typename T>
bool remove(std::vector<T>& vec, const T& element) {
  auto rem = std::remove(vec.begin(), vec.end(), element);
  auto rem2 = vec.erase(rem, vec.end());
  return rem2 != vec.end();
}

template <typename T>
HighsInt indexof(const std::vector<T>& vec, const T& element) {
  auto it = std::find(vec.begin(), vec.end(), element);
  if (it != vec.end()) {
    return distance(vec.begin(), it);
  } else {
    return -1;
  }
}

#endif

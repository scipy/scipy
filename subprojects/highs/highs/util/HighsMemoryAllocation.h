/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsMemoryAllocation.h
 * @brief Utilities for memory allocation that return true if successful
 */

#ifndef UTIL_HIGHS_MEMORY_ALLOCATION_H_
#define UTIL_HIGHS_MEMORY_ALLOCATION_H_

#include <vector>

#include "util/HighsType.h"

template <typename T>
bool okResize(std::vector<T>& use_vector, HighsInt dimension, T value = T{}) {
  try {
    use_vector.resize(dimension, value);
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::okResize fails with %s\n", e.what());
    return false;
  }
  return true;
}

template <typename T>
bool okReserve(std::vector<T>& use_vector, HighsInt dimension) {
  try {
    use_vector.reserve(dimension);
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::okReserve fails with %s\n", e.what());
    return false;
  }
  return true;
}

template <typename T, typename T2>
bool okReserve(std::unordered_map<T, T2>& use_map, HighsInt dimension) {
  try {
    use_map.reserve(dimension);
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::okReserve fails with %s\n", e.what());
    return false;
  }
  return true;
}

template <typename T>
bool okAssign(std::vector<T>& use_vector, HighsInt dimension, T value = T{}) {
  try {
    use_vector.assign(dimension, value);
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::okAssign fails with %s\n", e.what());
    return false;
  }
  return true;
}

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsDebug.h"

#include <algorithm>  // For std::max
#include <cassert>    // For std::max

HighsStatus debugDebugToHighsStatus(const HighsDebugStatus debug_status) {
  switch (debug_status) {
    case HighsDebugStatus::NOT_CHECKED:
    case HighsDebugStatus::OK:
    case HighsDebugStatus::SMALL_ERROR:
      return HighsStatus::OK;
    case HighsDebugStatus::WARNING:
    case HighsDebugStatus::LARGE_ERROR:
      return HighsStatus::Warning;
    case HighsDebugStatus::ERROR:
    case HighsDebugStatus::EXCESSIVE_ERROR:
    case HighsDebugStatus::LOGICAL_ERROR:
      return HighsStatus::Error;
    default:
      return HighsStatus::OK;
  }
}

HighsDebugStatus debugWorseStatus(const HighsDebugStatus status0,
                                  const HighsDebugStatus status1) {
  return static_cast<HighsDebugStatus>(std::max((int)status0, (int)status1));
}

bool debugVectorRightSize(const std::vector<double> v, const int right_size) {
  const int v_size = v.size();
  const bool is_right_size = v_size == right_size;
  assert(is_right_size);
  return is_right_size;
}

bool debugVectorRightSize(const std::vector<int> v, const int right_size) {
  const int v_size = v.size();
  const bool is_right_size = v_size == right_size;
  assert(is_right_size);
  return is_right_size;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkInterface.cpp
 * @brief
 */

#include "lp_data/HighsLpUtils.h"
#include "simplex/HEkk.h"

void HEkk::appendColsToVectors(const HighsInt num_new_col,
                               const vector<double>& colCost,
                               const vector<double>& colLower,
                               const vector<double>& colUpper) {
  appendColsToLpVectors(lp_, num_new_col, colCost, colLower, colUpper);
}

void HEkk::appendRowsToVectors(const HighsInt num_new_row,
                               const vector<double>& rowLower,
                               const vector<double>& rowUpper) {
  appendRowsToLpVectors(lp_, num_new_row, rowLower, rowUpper);
}

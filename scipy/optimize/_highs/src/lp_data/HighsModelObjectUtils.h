/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsModelObjectUtil.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHSMODELOBJECTUTILS_H_
#define LP_DATA_HIGHSMODELOBJECTUTILS_H_

#include <cassert>
#include <iostream>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsModelObject.h"

void report_row_vec_sol(int nrow, vector<double>& XrowLower,
                        vector<double>& XrowUpper, vector<double>& XrowPrimal,
                        vector<double>& XrowDual, vector<int>& XrowStatus) {
  // Report the LP row data and solution passed to the method, where
  // XrowStatus is the SCIP-like basis status
  if (nrow <= 0) return;
  printf("Row    St      Primal       Lower       Upper        Dual\n");
  for (int row = 0; row < nrow; row++) {
    if (XrowStatus[row] == (int)HighsBasisStatus::BASIC)
      printf("%6d BC", row);
    else if (XrowStatus[row] == (int)HighsBasisStatus::ZERO)
      printf("%6d FR", row);
    else if (XrowStatus[row] == (int)HighsBasisStatus::LOWER) {
      if (XrowLower[row] == XrowUpper[row])
        printf("%6d FX", row);
      else
        printf("%6d LB", row);
    } else if (XrowStatus[row] == (int)HighsBasisStatus::UPPER)
      printf("%6d UB", row);
    else
      printf("%6d ??", row);
    printf(" %11g %11g %11g %11g\n", XrowPrimal[row], XrowLower[row],
           XrowUpper[row], XrowDual[row]);
  }
}

void report_row_matrix(int nrow, vector<int>& XARstart, vector<int>& XARindex,
                       vector<double>& XARvalue) {
  // Report the row-wise matrix passed to the method
  if (nrow <= 0) return;
  printf("Row    Index       Value\n");
  for (int row = 0; row < nrow; row++) {
    printf("%6d Start %8d\n", row, XARstart[row]);
    for (int el = XARstart[row]; el < XARstart[row + 1]; el++) {
      printf("      %6d %11g\n", XARindex[el], XARvalue[el]);
    }
  }
  printf("       Start %8d\n", XARstart[nrow]);
}

void report_col_vec_sol(int ncol, vector<double>& XcolCost,
                        vector<double>& XcolLower, vector<double>& XcolUpper,
                        vector<double>& XcolPrimal, vector<double>& XcolDual,
                        vector<int>& XcolStatus) {
  // Report the LP column data and solution passed to the method,
  // where XcolStatus is the SCIP-like basis status
  if (ncol <= 0) return;
  printf(
      "Col    St      Primal       Lower       Upper        Dual        "
      "Cost\n");
  for (int col = 0; col < ncol; col++) {
    if (XcolStatus[col] == (int)HighsBasisStatus::BASIC)
      printf("%6d BC", col);
    else if (XcolStatus[col] == (int)HighsBasisStatus::ZERO)
      printf("%6d FR", col);
    else if (XcolStatus[col] == (int)HighsBasisStatus::LOWER) {
      if (colLower[col] == XcolUpper[col])
        printf("%6d FX", col);
      else
        printf("%6d LB", col);
    } else if (XcolStatus[col] == (int)HighsBasisStatus::UPPER)
      printf("%6d UB", col);
    else
      printf("%6d ??", col);
    printf(" %11g %11g %11g %11g %11g\n", XcolPrimal[col], colLower[col],
           XcolUpper[col], XcolDual[col], XcolCost[col]);
  }
}

#endif

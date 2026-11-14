/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HFactorUtils.cpp
 * @brief Types of solution classes
 */
#include "util/HFactor.h"

void HFactor::invalidAMatrixAction() {
  this->a_matrix_valid = false;
  refactor_info_.clear();
}

void HFactor::reportLu(const HighsInt l_u_or_both, const bool full) const {
  if (l_u_or_both < kReportLuJustL || l_u_or_both > kReportLuBoth) return;
  if (l_u_or_both & 1) {
    printf("L");
    if (full) printf(" - full");
    printf(":\n");

    if (full) reportIntVector("l_pivot_lookup", l_pivot_lookup);
    if (full) reportIntVector("l_pivot_index", l_pivot_index);
    reportIntVector("l_start", l_start);
    reportIntVector("l_index", l_index);
    reportDoubleVector("l_value", l_value);
    if (full) {
      reportIntVector("lr_start", lr_start);
      reportIntVector("lr_index", lr_index);
      reportDoubleVector("lr_value", lr_value);
    }
  }
  if (l_u_or_both & 2) {
    printf("U");
    if (full) printf(" - full");
    printf(":\n");
    if (full) reportIntVector("u_pivot_lookup", u_pivot_lookup);
    reportIntVector("u_pivot_index", u_pivot_index);
    reportDoubleVector("u_pivot_value", u_pivot_value);
    reportIntVector("u_start", u_start);
    if (full) reportIntVector("u_last_p", u_last_p);
    reportIntVector("u_index", u_index);
    reportDoubleVector("u_value", u_value);
    if (full) {
      reportIntVector("ur_start", ur_start);
      reportIntVector("ur_lastp", ur_lastp);
      reportIntVector("ur_space", ur_space);
      for (size_t iRow = 0; iRow < ur_start.size(); iRow++) {
        const HighsInt start = ur_start[iRow];
        const HighsInt end = ur_lastp[iRow];
        if (start >= end) continue;
        printf("UR    Row %2d: ", (int)iRow);
        for (HighsInt iEl = start; iEl < end; iEl++)
          printf("%11d ", (int)ur_index[iEl]);
        printf("\n              ");
        for (HighsInt iEl = start; iEl < end; iEl++)
          printf("%11.4g ", ur_value[iEl]);
        printf("\n");
      }
      //      reportIntVector("ur_index", ur_index);
      //      reportDoubleVector("ur_value", ur_value);
    }
  }
  if (l_u_or_both == 3 && full) {
    reportDoubleVector("pf_pivot_value", pf_pivot_value);
    reportIntVector("pf_pivot_index", pf_pivot_index);
    reportIntVector("pf_start", pf_start);
    reportIntVector("pf_index", pf_index);
    reportDoubleVector("pf_value", pf_value);
  }
}

void HFactor::reportIntVector(const std::string& name,
                              const vector<HighsInt>& entry) const {
  const HighsInt num_en = entry.size();
  printf("%-12s: siz %4d; cap %4d: ", name.c_str(), (int)num_en,
         (int)entry.capacity());
  for (HighsInt iEn = 0; iEn < num_en; iEn++) {
    if (iEn > 0 && iEn % 10 == 0)
      printf("\n                                  ");
    printf("%11d ", (int)entry[iEn]);
  }
  printf("\n");
}
void HFactor::reportDoubleVector(const std::string& name,
                                 const vector<double>& entry) const {
  const HighsInt num_en = entry.size();
  printf("%-12s: siz %4d; cap %4d: ", name.c_str(), (int)num_en,
         (int)entry.capacity());
  for (HighsInt iEn = 0; iEn < num_en; iEn++) {
    if (iEn > 0 && iEn % 10 == 0)
      printf("\n                                  ");
    printf("%11.4g ", entry[iEn]);
  }
  printf("\n");
}

void HFactor::reportAsm() const {
  for (HighsInt count = 1; count <= num_row; count++) {
    if (col_link_first[count] < 0) continue;
    for (HighsInt j = col_link_first[count]; j != -1; j = col_link_next[j]) {
      double min_pivot = mc_min_pivot[j];
      HighsInt start = mc_start[j];
      HighsInt end = start + mc_count_a[j];
      printf("Col %4d: count = %2d; min_pivot = %10.4g; [%4d, %4d)\n", (int)j,
             (int)count, min_pivot, (int)start, (int)end);
      for (HighsInt k = start; k < end; k++) {
        HighsInt i = mc_index[k];
        double value = mc_value[k];
        //	if (abs_value < 1e-8) continue;
        HighsInt row_count = mr_count[i];
        double merit_local = 1.0 * (count - 1) * (row_count - 1);
        printf("   Row %4d; Count = %2d; Merit = %11.4g; Value = %11.4g: %s\n",
               (int)i, (int)row_count, merit_local, value,
               std::abs(value) >= min_pivot ? "OK" : "");
      }
    }
  }
}

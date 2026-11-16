/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HFactorRefactor.cpp
 * @brief Types of solution classes
 */
#include <cassert>
#include <cmath>
#include <iostream>

#include "util/HFactor.h"
#include "util/HVectorBase.h"

// std::max and std::min used in HFactor.h for local in-line
// functions, so HFactor.h has #include <algorithm>
using std::fabs;

void RefactorInfo::clear() {
  this->use = false;
  this->build_synthetic_tick = 0.0;
  this->pivot_var.clear();
  this->pivot_row.clear();
  this->pivot_type.clear();
}

HighsInt HFactor::rebuild(HighsTimerClock* factor_timer_clock_pointer) {
  const bool report_lu = false;
  // Check that the refactorization information should be used
  assert(refactor_info_.use);
  /**
   * 0. Clear L and U factor
   */
  luClear();

  nwork = 0;
  basis_matrix_num_el = 0;
  HighsInt stage = num_row;
  HighsInt rank_deficiency = 0;
  std::vector<HighsBool> has_pivot;
  has_pivot.assign(num_row, false);
  const bool report_unit = false;
  const bool report_singletons = false;
  const bool report_markowitz = false;
  const bool report_anything =
      report_unit || report_singletons || report_markowitz;
  if (report_anything) printf("\nRefactor\n");
  // Take build_synthetic_tick from the refactor info so that this
  // refactorization doesn't look unrealistically cheap.
  this->build_synthetic_tick = this->refactor_info_.build_synthetic_tick;
  // Check that the refactorization info has been set up
  assert((int)this->refactor_info_.pivot_row.size() >= num_row);
  assert((int)this->refactor_info_.pivot_var.size() >= num_row);
  assert((int)this->refactor_info_.pivot_type.size() >= num_row);
  for (HighsInt iK = 0; iK < num_row; iK++) {
    HighsInt iRow = this->refactor_info_.pivot_row[iK];
    HighsInt iVar = this->refactor_info_.pivot_var[iK];
    int8_t pivot_type = this->refactor_info_.pivot_type[iK];
    assert(!has_pivot[iRow]);

    if (pivot_type == kPivotLogical || pivot_type == kPivotUnit) {
      if (pivot_type == kPivotLogical) {
        //
        // 1.1 Logical column
        if (report_unit) printf("Stage %d: Logical\n", (int)iK);
        assert(iVar >= num_col);
        basis_matrix_num_el++;
      } else if (pivot_type == kPivotUnit) {
        //
        // 1.2 (Structural) unit column
        if (report_unit) printf("Stage %d: Unit\n", (int)iK);
        assert(iVar < num_col);
        HighsInt start = a_start[iVar];
        HighsInt count = a_start[iVar + 1] - start;
        assert(a_index[start] == iRow);
        assert(count == 1 && a_value[start] == 1);
        basis_matrix_num_el++;
      }
      // 1.3 Record unit column
      l_start.push_back(l_index.size());
      u_pivot_index.push_back(iRow);
      u_pivot_value.push_back(1);
      u_start.push_back(u_index.size());
    } else if (pivot_type == kPivotRowSingleton ||
               pivot_type == kPivotColSingleton) {
      //
      // Row or column singleton
      assert(iVar < num_col);
      const HighsInt start = a_start[iVar];
      const HighsInt end = a_start[iVar + 1];
      // Find where the pivot is
      HighsInt pivot_k = -1;
      for (HighsInt k = start; k < end; k++) {
        if (a_index[k] == iRow) {
          pivot_k = k;
          break;
        }
      }
      assert(pivot_k >= 0);
      // Check that the pivot isn't too small. Shouldn't happen since
      // this is refactorization
      double abs_pivot = std::fabs(a_value[pivot_k]);
      assert(abs_pivot >= pivot_tolerance);
      if (abs_pivot < pivot_tolerance) {
        rank_deficiency = nwork + 1;
        return rank_deficiency;
      }
      if (pivot_type == kPivotRowSingleton) {
        //
        // 2.2 Deal with row singleton
        const double pivot_multiplier = 1 / a_value[pivot_k];
        if (report_singletons)
          printf("Stage %d: Row singleton (%4d, %g)\n", (int)iK, (int)pivot_k,
                 pivot_multiplier);
        for (HighsInt section = 0; section < 2; section++) {
          HighsInt p0 = section == 0 ? start : pivot_k + 1;
          HighsInt p1 = section == 0 ? pivot_k : end;
          for (HighsInt k = p0; k < p1; k++) {
            HighsInt local_iRow = a_index[k];
            if (!has_pivot[local_iRow]) {
              if (report_singletons)
                printf("Row singleton: L En (%4d, %11.4g)\n", (int)local_iRow,
                       a_value[k] * pivot_multiplier);
              l_index.push_back(local_iRow);
              l_value.push_back(a_value[k] * pivot_multiplier);
            } else {
              if (report_singletons)
                printf("Row singleton: U En (%4d, %11.4g)\n", (int)local_iRow,
                       a_value[k]);
              u_index.push_back(local_iRow);
              u_value.push_back(a_value[k]);
            }
          }
        }
        l_start.push_back(l_index.size());
        if (report_singletons)
          printf("Row singleton: U Pv (%4d, %11.4g)\n", (int)iRow,
                 a_value[pivot_k]);
        u_pivot_index.push_back(iRow);
        u_pivot_value.push_back(a_value[pivot_k]);
        u_start.push_back(u_index.size());
      } else {
        //
        // 2.3 Deal with column singleton
        if (report_singletons) printf("Stage %d: Col singleton\n", (int)iK);
        for (HighsInt k = start; k < pivot_k; k++) {
          if (report_singletons)
            printf("Col singleton: U En (%4d, %11.4g)\n", (int)a_index[k],
                   a_value[k]);
          u_index.push_back(a_index[k]);
          u_value.push_back(a_value[k]);
        }
        for (HighsInt k = pivot_k + 1; k < end; k++) {
          if (report_singletons)
            printf("Col singleton: U En (%4d, %11.4g)\n", (int)a_index[k],
                   a_value[k]);
          u_index.push_back(a_index[k]);
          u_value.push_back(a_value[k]);
        }
        l_start.push_back(l_index.size());
        if (report_singletons)
          printf("Col singleton: U Pv (%4d, %11.4g)\n", (int)iRow,
                 a_value[pivot_k]);
        u_pivot_index.push_back(iRow);
        u_pivot_value.push_back(a_value[pivot_k]);
        u_start.push_back(u_index.size());
      }
    } else {
      assert(pivot_type == kPivotMarkowitz);
      stage = iK;
      break;
    }
    basic_index[iRow] = iVar;
    has_pivot[iRow] = true;
  }
  if (report_lu) {
    printf("\nAfter units and singletons\n");
    reportLu(kReportLuBoth, false);
  }
  if (stage < num_row) {
    // Handle the remaining Markowitz pivots
    //
    // First of all complete the L factor with identity columns so
    // that FtranL counts the RHS entries in rows that don't yet have
    // pivots by running to completion. In the hyper-sparse code,
    // these will HOPEFULLY be skipped
    //
    // There are already l_start entries for the first stage rows, but
    // l_pivot_index is not assigned, as u_pivot_index gets copied into it
    l_start.resize(num_row + 1);
    for (HighsInt iK = stage; iK < num_row; iK++) l_start[iK + 1] = l_start[iK];
    l_pivot_index.resize(num_row);
    for (HighsInt iK = 0; iK < num_row; iK++)
      l_pivot_index[iK] = this->refactor_info_.pivot_row[iK];
    // To do hyper-sparse FtranL operations, have to set up l_pivot_lookup.
    l_pivot_lookup.resize(num_row);
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow < stage) {
        if (l_pivot_lookup[l_pivot_index[iRow]] != iRow) {
          // printf("Strange: Thought that l_pivot_lookup[l_pivot_index[iRow]]
          // == iRow\n");
        }
      }
      l_pivot_lookup[l_pivot_index[iRow]] = iRow;
    }
    // Need to know whether to consider matrix entries for FtranL
    // operation. Initially these correspond to all the rows without
    // pivots
    std::vector<HighsBool> not_in_bump = has_pivot;
    // Monitor density of FtranL result to possibly switch from exploiting
    // hyper-sparsity
    double expected_density = 0.0;
    // Initialise a HVector in which the L and U entries of the
    // pivotal column will be formed
    HVector column;
    column.setup(num_row);
    for (HighsInt iK = stage; iK < num_row; iK++) {
      HighsInt iRow = this->refactor_info_.pivot_row[iK];
      HighsInt iVar = this->refactor_info_.pivot_var[iK];
      int8_t pivot_type = this->refactor_info_.pivot_type[iK];
      assert(!has_pivot[iRow]);
      assert(pivot_type == kPivotMarkowitz);
      // Set up the column for the FtranL. It contains the matrix
      // entries in rows without pivots, and the remaining entries
      // start forming the U column
      column.clear();
      HighsInt start = a_start[iVar];
      HighsInt end = a_start[iVar + 1];
      for (HighsInt iEl = start; iEl < end; iEl++) {
        HighsInt local_iRow = a_index[iEl];
        if (not_in_bump[local_iRow]) {
          u_index.push_back(local_iRow);
          u_value.push_back(a_value[iEl]);
        } else {
          column.index[column.count++] = local_iRow;
          column.array[local_iRow] = a_value[iEl];
        }
      }
      // Perform FtranL, but don't time it!
      ftranL(column, expected_density);
      // Update the running average density
      double local_density = (1.0 * column.count) / num_row;
      expected_density = kRunningAverageMultiplier * local_density +
                         (1 - kRunningAverageMultiplier) * expected_density;
      // Strip out small values
      column.tight();
      // Now form the column of L
      //
      // Find the pivot
      HighsInt pivot_k = -1;
      start = 0;
      end = column.count;
      for (HighsInt k = start; k < end; k++) {
        if (column.index[k] == iRow) {
          pivot_k = k;
          break;
        }
      }
      assert(pivot_k >= 0);
      // Check that the pivot isn't too small. Shouldn't happen since
      // this is refactorization
      double abs_pivot = std::fabs(column.array[iRow]);
      assert(abs_pivot >= pivot_tolerance);
      if (abs_pivot < pivot_tolerance) {
        rank_deficiency = num_row - iK;
        return rank_deficiency;
      }
      const double pivot_multiplier = 1 / column.array[iRow];
      for (HighsInt section = 0; section < 2; section++) {
        HighsInt p0 = section == 0 ? start : pivot_k + 1;
        HighsInt p1 = section == 0 ? pivot_k : end;
        for (HighsInt k = p0; k < p1; k++) {
          HighsInt local_iRow = column.index[k];
          if (!has_pivot[local_iRow]) {
            l_index.push_back(local_iRow);
            l_value.push_back(column.array[local_iRow] * pivot_multiplier);
          } else {
            u_index.push_back(local_iRow);
            u_value.push_back(column.array[local_iRow]);
          }
        }
      }
      l_start[iK + 1] = l_index.size();
      u_pivot_index.push_back(iRow);
      u_pivot_value.push_back(column.array[iRow]);
      u_start.push_back(u_index.size());
      basic_index[iRow] = iVar;
      has_pivot[iRow] = true;
      if (report_lu) {
        printf("\nAfter Markowitz %d\n", (int)(iK - stage));
        reportLu(kReportLuBoth, false);
      }
    }
  }
  if (report_lu) {
    printf("\nRefactored INVERT\n");
    reportLu(kReportLuBoth, false);
  }
  buildFinish();
  return 0;
}

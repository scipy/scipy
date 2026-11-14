/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HFactorExtend.cpp
 * @brief Types of solution classes
 */
#include <cassert>

#include "util/HFactor.h"
#include "util/HVectorBase.h"

using std::fabs;

void HFactor::addCols(const HighsInt num_new_col) {
  invalidAMatrixAction();
  num_col += num_new_col;
}

void HFactor::deleteNonbasicCols(const HighsInt num_deleted_col) {
  invalidAMatrixAction();
  num_col -= num_deleted_col;
}

void HFactor::addRows(const HighsSparseMatrix* ar_matrix) {
  invalidAMatrixAction();
  assert(kExtendInvertWhenAddingRows);
  HighsInt num_new_row = ar_matrix->num_row_;
  HighsInt new_num_row = num_row + num_new_row;
  printf("Adding %" HIGHSINT_FORMAT
         " new rows to HFactor instance: increasing dimension from "
         "%" HIGHSINT_FORMAT " to %" HIGHSINT_FORMAT " \n",
         num_new_row, num_row, new_num_row);

  // Need to know where (if) a column is basic
  vector<HighsInt> in_basis;
  in_basis.assign(num_col, -1);
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = basic_index[iRow];
    if (iVar < num_col) in_basis[iVar] = iRow;
  }
  for (HighsInt iRow = num_row; iRow < new_num_row; iRow++) {
    HighsInt iVar = basic_index[iRow];
    assert(iVar >= num_col);
  }
  //  reportLu(kReportLuBoth, true);
  //  reportLu(kReportLuJustL);

  // Create a row-wise sparse matrix containing the new rows of the L
  // matrix - so that a column-wise version can be created (after
  // inserting the rows into the LR matrix) allowing the cew column
  // entries to be inserted efficiently into the L matrix
  HighsSparseMatrix new_lr_rows;
  new_lr_rows.format_ = MatrixFormat::kRowwise;
  new_lr_rows.num_col_ = num_row;
  double expected_density = 0.0;
  HVector rhs;
  rhs.setup(num_row);
  this->lr_start.reserve(new_num_row + 1);
  for (HighsInt inewRow = 0; inewRow < num_new_row; inewRow++) {
    //    printf("\nFor new row %" HIGHSINT_FORMAT "\n", inewRow);
    // Prepare RHS for system U^T.v = r
    rhs.clear();
    rhs.packFlag = true;
    for (HighsInt iEl = ar_matrix->start_[inewRow];
         iEl < ar_matrix->start_[inewRow + 1]; iEl++) {
      HighsInt iCol = ar_matrix->index_[iEl];
      HighsInt basis_index = in_basis[iCol];
      if (basis_index >= 0) {
        rhs.array[basis_index] = ar_matrix->value_[iEl];
        rhs.index[rhs.count++] = basis_index;
      }
    }
    // Solve U^T.v = r
    btranU(rhs, expected_density);
    double local_density = (1.0 * rhs.count) / num_row;
    expected_density = kRunningAverageMultiplier * local_density +
                       (1 - kRunningAverageMultiplier) * expected_density;
    //    printf("New row btranU density: local = %11.4g; expected =  %11.4g\n",
    //    local_density, expected_density);
    rhs.tight();
    //
    // Append v to the matrix containing the new rows of L
    HighsInt rhs_num_nz = rhs.count;
    //    printf("New row has entries:\n");
    for (HighsInt iX = 0; iX < rhs_num_nz; iX++) {
      HighsInt iCol = rhs.index[iX];
      new_lr_rows.index_.push_back(iCol);
      new_lr_rows.value_.push_back(rhs.array[iCol]);
      //      printf("   %4d, %11.4g\n", (int)iCol, rhs.array[iCol]);
    }
    new_lr_rows.start_.push_back(new_lr_rows.index_.size());
    new_lr_rows.num_row_++;

    //
    // Append v to the L matrix
    for (HighsInt iX = 0; iX < rhs_num_nz; iX++) {
      HighsInt iCol = rhs.index[iX];
      lr_index.push_back(iCol);
      lr_value.push_back(rhs.array[iCol]);
    }
    lr_start.push_back(lr_index.size());
    //    reportLu(kReportLuJustL);
  }
  // Now create a column-wise copy of the new rows
  HighsSparseMatrix new_lr_cols = new_lr_rows.getColwise();
  assert(new_lr_cols.isColwise());
  // new_lr_cols.ensureColwise(); is now redundant

  //
  // Insert the column-wise copy into the L matrix
  //
  // Add pivot indices for the new columns
  this->l_pivot_index.resize(new_num_row);
  for (HighsInt iCol = num_row; iCol < new_num_row; iCol++)
    l_pivot_index[iCol] = iCol;
  //
  // Add starts for the identity columns
  HighsInt l_matrix_new_num_nz = lr_index.size();
  assert(l_matrix_new_num_nz ==
         (HighsInt)(l_index.size() + new_lr_cols.index_.size()));
  l_start.resize(new_num_row + 1);
  HighsInt to_el = l_matrix_new_num_nz;
  for (HighsInt iCol = num_row + 1; iCol < new_num_row + 1; iCol++)
    l_start[iCol] = l_matrix_new_num_nz;
  //
  // Insert the new entries, remembering to offset the index values by
  // num_row, since new_lr_cols only has the new rows
  l_index.resize(l_matrix_new_num_nz);
  l_value.resize(l_matrix_new_num_nz);
  for (HighsInt iCol = num_row - 1; iCol >= 0; iCol--) {
    const HighsInt from_el = l_start[iCol + 1];
    l_start[iCol + 1] = to_el;
    for (HighsInt iEl = new_lr_cols.start_[iCol + 1] - 1;
         iEl >= new_lr_cols.start_[iCol]; iEl--) {
      to_el--;
      l_index[to_el] = num_row + new_lr_cols.index_[iEl];
      l_value[to_el] = new_lr_cols.value_[iEl];
    }
    for (HighsInt iEl = from_el - 1; iEl >= l_start[iCol]; iEl--) {
      to_el--;
      l_index[to_el] = l_index[iEl];
      l_value[to_el] = l_value[iEl];
    }
  }
  assert(to_el == 0);
  this->l_pivot_lookup.resize(new_num_row);
  for (HighsInt iRow = num_row; iRow < new_num_row; iRow++)
    l_pivot_lookup[l_pivot_index[iRow]] = iRow;
  // Now update the U matrix with identity rows and columns
  // Allocate space for U factor
  //  reportLu(kReportLuJustL);
  //
  // Now add pivots corresponding to identity columns in U. The use of
  // most arrays in U is self-evident, except for U(R)lastp, which is
  // (one more than) the end of the index/value array that's applied
  // in the *tranU. Here it needs to be equal to the start since there
  // are no non-pivotal entries
  //

  HighsInt u_countX = u_index.size();
  HighsInt u_pivot_lookup_offset = u_pivot_index.size() - num_row;
  for (HighsInt iRow = num_row; iRow < new_num_row; iRow++) {
    u_pivot_lookup.push_back(u_pivot_lookup_offset + iRow);
    u_pivot_index.push_back(iRow);
    u_pivot_value.push_back(1);
    u_start.push_back(u_countX);
    u_last_p.push_back(u_countX);
  }

  // Now, to extend UR
  //
  // UR space
  //
  // Borrowing names from buildFinish()
  HighsInt ur_stuff_size = update_method == kUpdateMethodFt ? 5 : 0;
  HighsInt ur_size = ur_index.size();
  HighsInt ur_count_size = ur_size + ur_stuff_size * num_new_row;
  ur_index.resize(ur_count_size);
  ur_value.resize(ur_count_size);

  // Need to refer to just the new UR vectors
  HighsInt ur_cur_num_vec = ur_start.size();
  HighsInt ur_new_num_vec = ur_cur_num_vec + num_new_row;
  printf("\nUpdating UR vectors %d - %d\n", (int)ur_cur_num_vec,
         (int)ur_new_num_vec - 1);
  // UR pointer
  //
  // Allow space to the start of new rows, including the start for the
  // fictitious ur_new_num_vec'th row
  ur_start.resize(ur_new_num_vec + 1);
  for (HighsInt iRow = ur_cur_num_vec + 1; iRow < ur_new_num_vec + 1; iRow++) {
    ur_start[iRow] = ur_size;
  }
  // NB ur_temp plays the role of ur_lastp when it could be used as
  // temporary storage in buildFinish()
  vector<HighsInt> ur_temp;
  ur_temp.assign(ur_new_num_vec, 0);
  //
  // ur_space has its new entries assigned (to ur_stuff_size) as in
  // buildFinish()
  ur_space.resize(ur_new_num_vec);
  for (HighsInt iRow = ur_cur_num_vec; iRow < ur_new_num_vec; iRow++)
    ur_space[iRow] = ur_stuff_size;
  // Compute ur_temp exactly as in buildFinish()
  for (HighsInt k = 0; k < u_countX; k++) ur_temp[u_pivot_lookup[u_index[k]]]++;
  HighsInt iStart = ur_size;
  ur_start[ur_cur_num_vec] = iStart;
  for (HighsInt iRow = ur_cur_num_vec + 1; iRow < ur_new_num_vec + 1; iRow++) {
    HighsInt gap = ur_temp[iRow - 1] + ur_stuff_size;
    ur_start[iRow] = iStart + gap;
    iStart += gap;
    printf("ur_start[%d] = %d; gap = %d; iStart = %d\n", (int)iRow,
           (int)ur_start[iRow], (int)gap, (int)iStart);
  }
  printf("ur_count_size = %d; iStart%d\n", (int)ur_count_size, (int)iStart);
  // Lose the start for the fictitious ur_new_num_vec'th row
  ur_start.resize(ur_new_num_vec);
  // Resize ur_lastp and initialise its new entries to be the ur_start
  // values since the rows are empty
  ur_lastp.resize(ur_new_num_vec);
  for (HighsInt iRow = ur_cur_num_vec; iRow < ur_new_num_vec; iRow++)
    ur_lastp[iRow] = ur_start[iRow];
  //
  // Increase the number of rows in HFactor
  num_row += num_new_row;
  //  reportLu(kReportLuBoth, true);
}

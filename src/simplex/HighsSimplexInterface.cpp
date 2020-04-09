/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexInterface.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HighsSimplexInterface.h"

#include "HConfig.h"
#include "io/HMPSIO.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "simplex/HSimplex.h"
#include "util/HighsSort.h"
#include "util/HighsUtils.h"

HighsStatus HighsSimplexInterface::addCols(
    int XnumNewCol, const double* XcolCost, const double* XcolLower,
    const double* XcolUpper, int XnumNewNZ, const int* XAstart,
    const int* XAindex, const double* XAvalue) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
#ifdef HiGHSDEV
  printf("Called addCols(XnumNewCol=%d, XnumNewNZ = %d)\n", XnumNewCol,
         XnumNewNZ);
#endif
  HighsOptions& options = highs_model_object.options_;
  if (XnumNewCol < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewCol == 0) return HighsStatus::OK;
  if (XnumNewCol > 0)
    if (isColDataNull(options, XcolCost, XcolLower, XcolUpper))
      return HighsStatus::Error;
  if (XnumNewNZ > 0)
    if (isMatrixDataNull(options, XAstart, XAindex, XAvalue))
      return HighsStatus::Error;

  HighsLp& lp = highs_model_object.lp_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsScale& scale = highs_model_object.scale_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_basis = basis.valid_;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_basis = simplex_lp_status.has_basis;
  bool apply_row_scaling = scale.is_scaled_;

  // Check that if nonzeros are to be added then the model has a positive number
  // of rows
  if (lp.numRow_ <= 0 && XnumNewNZ > 0) return HighsStatus::Error;
  if (valid_simplex_lp && (simplex_lp.numRow_ <= 0 && XnumNewNZ > 0))
    return HighsStatus::Error;

  // Record the new number of columns
  int newNumCol = lp.numCol_ + XnumNewCol;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no basis, matrix or
  // scaling
  if (!valid_simplex_lp) {
    assert(!apply_row_scaling);
  }
#endif
  call_status = appendLpCols(options, lp, XnumNewCol, XcolCost, XcolLower,
                             XcolUpper, XnumNewNZ, XAstart, XAindex, XAvalue);
  return_status =
      interpretCallStatus(call_status, return_status, "appendLpCols");
  if (return_status == HighsStatus::Error) return return_status;

  if (valid_simplex_lp) {
    call_status =
        appendLpCols(options, simplex_lp, XnumNewCol, XcolCost, XcolLower,
                     XcolUpper, XnumNewNZ, XAstart, XAindex, XAvalue);
    return_status =
        interpretCallStatus(call_status, return_status, "appendLpCols");
    if (return_status == HighsStatus::Error) return return_status;
  }

  // Now consider scaling
  scale.col_.resize(newNumCol);
  for (int col = 0; col < XnumNewCol; col++)
    scale.col_[simplex_lp.numCol_ + col] = 1.0;

  if (apply_row_scaling) {
    // Determine scaling multipliers for this set of columns
    // Determine scale factors for this set of columns
    // Scale the simplex LP vectors for these columns
    // Scale the simplex LP matrix for these columns
  }

  // Update the basis correponding to new nonbasic columns
  if (valid_basis) append_nonbasic_cols_to_basis(lp, basis, XnumNewCol);
  if (valid_simplex_basis)
    append_nonbasic_cols_to_basis(simplex_lp, simplex_basis, XnumNewCol);

  // Deduce the consequences of adding new columns
  highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
  highs_model_object.unscaled_model_status_ =
      highs_model_object.scaled_model_status_;
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_COLS);

  // Increase the number of columns in the LPs
  lp.numCol_ += XnumNewCol;
  if (valid_simplex_lp) simplex_lp.numCol_ += XnumNewCol;

#ifdef HiGHSDEV
  if (valid_basis) {
    bool basis_ok = basisOk(options.logfile, lp, basis);
    if (!basis_ok) printf("HiGHS basis not OK in addCols\n");
    assert(basis_ok);
    report_basis(lp, basis);
  }
  if (valid_simplex_basis) {
    bool basis_ok = basisOk(options.logfile, simplex_lp, simplex_basis);
    if (!basis_ok) printf("Simplex basis not OK in addCols\n");
    assert(basis_ok);
    report_basis(simplex_lp, simplex_basis);
  }
#endif
  return return_status;
}

HighsStatus HighsSimplexInterface::deleteCols(int from_col, int to_col) {
  return deleteColsGeneral(true, from_col, to_col, false, 0, NULL, false, NULL);
}

HighsStatus HighsSimplexInterface::deleteCols(int num_set_entries,
                                              const int* col_set) {
  return deleteColsGeneral(false, 0, 0, true, num_set_entries, col_set, false,
                           NULL);
}

HighsStatus HighsSimplexInterface::deleteCols(int* col_mask) {
  return deleteColsGeneral(false, 0, 0, false, 0, NULL, true, col_mask);
}

HighsStatus HighsSimplexInterface::deleteColsGeneral(
    bool interval, int from_col, int to_col, bool set, int num_set_entries,
    const int* col_set, bool mask, int* col_mask) {
  HighsOptions& options = highs_model_object.options_;
  HighsLp& lp = highs_model_object.lp_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_simplex_lp = simplex_lp_status.valid;
  // Keep a copy of the original number of columns to check whether
  // any columns have been removed, and if there is mask to be updated
  int original_num_col = lp.numCol_;

  HighsStatus returnStatus;
  returnStatus = deleteLpCols(options, lp, interval, from_col, to_col, set,
                              num_set_entries, col_set, mask, col_mask);
  if (returnStatus != HighsStatus::OK) return returnStatus;
  assert(lp.numCol_ <= original_num_col);
  if (lp.numCol_ < original_num_col) {
    // Nontrivial deletion so reset the model_status and invalidate
    // the Highs basis
    highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
    highs_model_object.unscaled_model_status_ =
        highs_model_object.scaled_model_status_;
    basis.valid_ = false;
  }
  if (valid_simplex_lp) {
    HighsLp& simplex_lp = highs_model_object.simplex_lp_;
    //  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
    returnStatus = deleteLpCols(options, simplex_lp, interval, from_col, to_col,
                                set, num_set_entries, col_set, mask, col_mask);
    if (returnStatus != HighsStatus::OK) return returnStatus;
    //    HighsScale& scale = highs_model_object.scale_;
    //    for (int col = from_col; col < lp.numCol_ - numDeleteCol; col++)
    //    scale.col_[col] = scale.col_[col + numDeleteCol];
    assert(simplex_lp.numCol_ <= original_num_col);
    if (simplex_lp.numCol_ < original_num_col) {
      // Nontrivial deletion so invalidate all data relating to the simplex
      // basis
      invalidateSimplexLpBasis(simplex_lp_status);
    }
  }
  if (mask) {
    int new_col = 0;
    for (int col = 0; col < original_num_col; col++) {
      if (!col_mask[col]) {
        col_mask[col] = new_col;
        new_col++;
      } else {
        col_mask[col] = -1;
      }
    }
    assert(new_col == lp.numCol_);
  }
  return HighsStatus::OK;
}

// Get a single coefficient from the matrix
HighsStatus HighsSimplexInterface::getCoefficient(const int Xrow,
                                                  const int Xcol,
                                                  double& value) {
#ifdef HiGHSDEV
  printf("Called getCoeff(Xrow=%d, Xcol=%d)\n", Xrow, Xcol);
#endif
  HighsLp& lp = highs_model_object.lp_;
  if (Xrow < 0 || Xrow > lp.numRow_) return HighsStatus::Error;
  if (Xcol < 0 || Xcol > lp.numCol_) return HighsStatus::Error;
  value = 0;
  for (int el = lp.Astart_[Xcol]; el < lp.Astart_[Xcol + 1]; el++) {
    if (lp.Aindex_[el] == Xrow) {
      value = lp.Avalue_[el];
      break;
    }
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::addRows(int XnumNewRow,
                                           const double* XrowLower,
                                           const double* XrowUpper,
                                           int XnumNewNZ, const int* XARstart,
                                           const int* XARindex,
                                           const double* XARvalue) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
#ifdef HiGHSDEV
  printf("Called addRows(XnumNewRow=%d, XnumNewNZ = %d)\n", XnumNewRow,
         XnumNewNZ);
#endif
  HighsOptions& options = highs_model_object.options_;
  if (XnumNewRow < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewRow == 0) return HighsStatus::OK;
  if (XnumNewRow > 0)
    if (isRowDataNull(options, XrowLower, XrowUpper)) return HighsStatus::Error;
  if (XnumNewNZ > 0)
    if (isMatrixDataNull(options, XARstart, XARindex, XARvalue))
      return HighsStatus::Error;

  HighsLp& lp = highs_model_object.lp_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsScale& scale = highs_model_object.scale_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_basis = basis.valid_;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_basis = simplex_lp_status.has_basis;
  bool apply_row_scaling = scale.is_scaled_;

  // Check that if nonzeros are to be added then the model has a positive number
  // of columns
  if (lp.numCol_ <= 0 && XnumNewNZ > 0) return HighsStatus::Error;
  if (valid_simplex_lp && (simplex_lp.numCol_ <= 0 && XnumNewNZ > 0))
    return HighsStatus::Error;

  // Record the new number of rows
  int newNumRow = lp.numRow_ + XnumNewRow;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no basis, matrix or
  // scaling
  if (!valid_simplex_lp) {
    assert(!apply_row_scaling);
  }
#endif
  // Assess the bounds and matrix indices, returning on error
  bool normalise = false;
  call_status = assessBounds(options, "Row", lp.numRow_, XnumNewRow, true, 0,
                             XnumNewRow - 1, false, 0, NULL, false, NULL,
                             (double*)XrowLower, (double*)XrowUpper,
                             options.infinite_bound, normalise);
  return_status =
      interpretCallStatus(call_status, return_status, "assessBounds");
  if (return_status == HighsStatus::Error) return return_status;

  if (XnumNewNZ) {
    call_status = assessMatrix(
        options, lp.numCol_, 0, XnumNewRow - 1, XnumNewRow, XnumNewNZ,
        (int*)XARstart, (int*)XARindex, (double*)XARvalue,
        options.small_matrix_value, options.large_matrix_value, normalise);
    return_status =
        interpretCallStatus(call_status, return_status, "assessMatrix");
    if (return_status == HighsStatus::Error) return return_status;
  }

  // Append the columns to the LP vectors and matrix
  appendRowsToLpVectors(lp, XnumNewRow, XrowLower, XrowUpper);

  // Normalise the LP row bounds
  normalise = true;
  call_status =
      assessBounds(options, "Row", lp.numRow_, newNumRow, true, 0,
                   newNumRow - 1, false, 0, NULL, false, NULL, &lp.rowLower_[0],
                   &lp.rowUpper_[0], options.infinite_bound, normalise);
  return_status =
      interpretCallStatus(call_status, return_status, "assessBounds");
  if (return_status == HighsStatus::Error) return return_status;

  int lc_XnumNewNZ = XnumNewNZ;
  int* lc_XARstart = (int*)malloc(sizeof(int) * XnumNewRow);
  int* lc_XARindex = (int*)malloc(sizeof(int) * XnumNewNZ);
  double* lc_XARvalue = (double*)malloc(sizeof(double) * XnumNewNZ);
  if (XnumNewNZ) {
    // Copy the new row-wise matrix into a local copy that can be normalised
    std::memcpy(lc_XARstart, XARstart, sizeof(int) * XnumNewRow);
    std::memcpy(lc_XARindex, XARindex, sizeof(int) * XnumNewNZ);
    std::memcpy(lc_XARvalue, XARvalue, sizeof(double) * XnumNewNZ);
    // Normalise the new matrix columns
    normalise = true;
    call_status = assessMatrix(
        options, lp.numCol_, 0, XnumNewRow - 1, XnumNewRow, lc_XnumNewNZ,
        lc_XARstart, lc_XARindex, lc_XARvalue, options.small_matrix_value,
        options.large_matrix_value, normalise);
    if (lc_XnumNewNZ) {
      // Append rows to LP matrix
      appendRowsToLpMatrix(lp, XnumNewRow, lc_XnumNewNZ, lc_XARstart,
                           lc_XARindex, lc_XARvalue);
    }
  }

  if (valid_simplex_lp) {
    appendRowsToLpVectors(simplex_lp, XnumNewRow, XrowLower, XrowUpper);
    call_status = assessBounds(
        options, "Row", simplex_lp.numRow_, newNumRow, true, 0, newNumRow - 1,
        false, 0, NULL, false, NULL, &simplex_lp.rowLower_[0],
        &simplex_lp.rowUpper_[0], options.infinite_bound, normalise);
    return_status =
        interpretCallStatus(call_status, return_status, "assessBounds");
    if (return_status == HighsStatus::Error) return return_status;
  }
  if (lc_XnumNewNZ) {
    appendRowsToLpMatrix(simplex_lp, XnumNewRow, lc_XnumNewNZ, lc_XARstart,
                         lc_XARindex, lc_XARvalue);
  }

  // Now consider scaling
  scale.row_.resize(newNumRow);
  for (int row = 0; row < XnumNewRow; row++) scale.row_[lp.numRow_ + row] = 1.0;

  if (apply_row_scaling) {
    // Determine scaling multipliers for this set of rows
    // Determine scale factors for this set of rows
    // Scale the simplex LP vectors for these rows
    // Scale the simplex LP matrix for these rows
  }

  // Update the basis correponding to new basic rows
  if (valid_basis) append_basic_rows_to_basis(lp, basis, XnumNewRow);
  if (valid_simplex_basis)
    append_basic_rows_to_basis(simplex_lp, simplex_basis, XnumNewRow);

  // Deduce the consequences of adding new rows
  highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
  highs_model_object.unscaled_model_status_ =
      highs_model_object.scaled_model_status_;
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_ROWS);

  // Increase the number of rows in the LPs
  lp.numRow_ += XnumNewRow;
  if (valid_simplex_lp) simplex_lp.numRow_ += XnumNewRow;

#ifdef HiGHSDEV
  if (valid_basis) {
    bool basis_ok = basisOk(options.logfile, lp, basis);
    if (!basis_ok) printf("HiGHS basis not OK in addRows\n");
    assert(basis_ok);
    report_basis(lp, basis);
  }
  if (valid_simplex_basis) {
    bool basis_ok = basisOk(options.logfile, simplex_lp, simplex_basis);
    if (!basis_ok) printf("Simplex basis not OK in addRows\n");
    assert(basis_ok);
    report_basis(simplex_lp, simplex_basis);
  }
#endif
  free(lc_XARstart);
  free(lc_XARindex);
  free(lc_XARvalue);
  return return_status;
}

HighsStatus HighsSimplexInterface::deleteRows(int from_row, int to_row) {
  return deleteRowsGeneral(true, from_row, to_row, false, 0, NULL, false, NULL);
}

HighsStatus HighsSimplexInterface::deleteRows(int num_set_entries,
                                              const int* row_set) {
  return deleteRowsGeneral(false, 0, 0, true, num_set_entries, row_set, false,
                           NULL);
}

HighsStatus HighsSimplexInterface::deleteRows(int* row_mask) {
  return deleteRowsGeneral(false, 0, 0, false, 0, NULL, true, row_mask);
}

HighsStatus HighsSimplexInterface::deleteRowsGeneral(
    bool interval, int from_row, int to_row, bool set, int num_set_entries,
    const int* row_set, bool mask, int* row_mask) {
#ifdef HiGHSDEV
  printf("Called model.util_deleteRows(from_row=%d, to_row=%d)\n", from_row,
         to_row);
#endif
  HighsOptions& options = highs_model_object.options_;
  HighsLp& lp = highs_model_object.lp_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_simplex_lp = simplex_lp_status.valid;
  // Keep a copy of the original number of rows to check whether
  // any rows have been removed, and if there is mask to be updated
  int original_num_row = lp.numRow_;

  HighsStatus returnStatus;
  returnStatus = deleteLpRows(options, lp, interval, from_row, to_row, set,
                              num_set_entries, row_set, mask, row_mask);
  if (returnStatus != HighsStatus::OK) return returnStatus;
  assert(lp.numRow_ <= original_num_row);
  if (lp.numRow_ < original_num_row) {
    // Nontrivial deletion so reset the model_status and invalidate
    // the Highs basis
    highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
    highs_model_object.unscaled_model_status_ =
        highs_model_object.scaled_model_status_;
    basis.valid_ = false;
  }
  if (valid_simplex_lp) {
    HighsLp& simplex_lp = highs_model_object.simplex_lp_;
    //    SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
    returnStatus = deleteLpRows(options, simplex_lp, interval, from_row, to_row,
                                set, num_set_entries, row_set, mask, row_mask);
    if (returnStatus != HighsStatus::OK) return returnStatus;
    //    HighsScale& scale = highs_model_object.scale_;
    //    for (int row = from_row; row < lp.numRow_ - numDeleteRow; row++)
    //    scale.row_[row] = scale.row_[row + numDeleteRow];
    assert(simplex_lp.numRow_ <= original_num_row);
    if (simplex_lp.numRow_ < original_num_row) {
      // Nontrivial deletion so invalidate all data relating to the simplex
      // basis
      invalidateSimplexLpBasis(simplex_lp_status);
    }
  }
  if (mask) {
    int new_row = 0;
    for (int row = 0; row < original_num_row; row++) {
      if (!row_mask[row]) {
        row_mask[row] = new_row;
        new_row++;
      } else {
        row_mask[row] = -1;
      }
    }
    assert(new_row == lp.numRow_);
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::getCols(const int from_col, const int to_col,
                                           int& num_col, double* col_cost,
                                           double* col_lower, double* col_upper,
                                           int& num_nz, int* col_matrix_start,
                                           int* col_matrix_index,
                                           double* col_matrix_value) {
  return getColsGeneral(true, from_col, to_col, false, 0, NULL, false, NULL,
                        num_col, col_cost, col_lower, col_upper, num_nz,
                        col_matrix_start, col_matrix_index, col_matrix_value);
}

HighsStatus HighsSimplexInterface::getCols(
    const int num_set_entries, const int* col_set, int& num_col,
    double* col_cost, double* col_lower, double* col_upper, int& num_nz,
    int* col_matrix_start, int* col_matrix_index, double* col_matrix_value) {
  return getColsGeneral(false, 0, 0, true, num_set_entries, col_set, false,
                        NULL, num_col, col_cost, col_lower, col_upper, num_nz,
                        col_matrix_start, col_matrix_index, col_matrix_value);
}

HighsStatus HighsSimplexInterface::getCols(const int* col_mask, int& num_col,
                                           double* col_cost, double* col_lower,
                                           double* col_upper, int& num_nz,
                                           int* col_matrix_start,
                                           int* col_matrix_index,
                                           double* col_matrix_value) {
  return getColsGeneral(false, 0, 0, false, 0, NULL, true, col_mask, num_col,
                        col_cost, col_lower, col_upper, num_nz,
                        col_matrix_start, col_matrix_index, col_matrix_value);
}

HighsStatus HighsSimplexInterface::getColsGeneral(
    const bool interval, const int from_col, const int to_col, const bool set,
    const int num_set_entries, const int* col_set, const bool mask,
    const int* col_mask, int& num_col, double* col_cost, double* col_lower,
    double* col_upper, int& num_nz, int* col_matrix_start,
    int* col_matrix_index, double* col_matrix_value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  int from_k;
  int to_k;
  HighsLp& lp = highs_model_object.lp_;
  HighsOptions& options = highs_model_object.options_;
  call_status = assessIntervalSetMask(options, lp.numCol_, interval, from_col,
                                      to_col, set, num_set_entries, col_set,
                                      mask, col_mask, from_k, to_k);
  return_status =
      interpretCallStatus(call_status, return_status, "assessIntervalSetMask");
  if (return_status == HighsStatus::Error) return return_status;
  if (from_k < 0 || to_k > lp.numCol_) {
    call_status = HighsStatus::Error;
    return_status =
        interpretCallStatus(call_status, return_status, "getColsGeneral");
    return return_status;
  }
  num_col = 0;
  num_nz = 0;
  if (from_k > to_k) {
    call_status = HighsStatus::Error;
    return_status =
        interpretCallStatus(call_status, return_status, "getColsGeneral");
    return return_status;
  }
  int out_from_col;
  int out_to_col;
  int in_from_col;
  int in_to_col = -1;
  int current_set_entry = 0;
  int col_dim = lp.numCol_;
  for (int k = from_k; k <= to_k; k++) {
    updateOutInIx(col_dim, interval, from_col, to_col, set, num_set_entries,
                  col_set, mask, col_mask, out_from_col, out_to_col,
                  in_from_col, in_to_col, current_set_entry);
    assert(out_to_col < col_dim);
    assert(in_to_col < col_dim);
    for (int col = out_from_col; col <= out_to_col; col++) {
      if (col_cost != NULL) col_cost[num_col] = lp.colCost_[col];
      if (col_lower != NULL) col_lower[num_col] = lp.colLower_[col];
      if (col_upper != NULL) col_upper[num_col] = lp.colUpper_[col];
      if (col_matrix_start != NULL)
        col_matrix_start[num_col] =
            num_nz + lp.Astart_[col] - lp.Astart_[out_from_col];
      num_col++;
    }
    if (col_matrix_index != NULL || col_matrix_value != NULL) {
      for (int el = lp.Astart_[out_from_col]; el < lp.Astart_[out_to_col + 1];
           el++) {
        if (col_matrix_index != NULL) col_matrix_index[num_nz] = lp.Aindex_[el];
        if (col_matrix_value != NULL) col_matrix_value[num_nz] = lp.Avalue_[el];
        num_nz++;
      }
    }
    if (out_to_col == col_dim - 1 || in_to_col == col_dim - 1) break;
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::getRows(const int from_row, const int to_row,
                                           int& num_row, double* row_lower,
                                           double* row_upper, int& num_nz,
                                           int* row_matrix_start,
                                           int* row_matrix_index,
                                           double* row_matrix_value) {
  return getRowsGeneral(true, from_row, to_row, false, 0, NULL, false, NULL,
                        num_row, row_lower, row_upper, num_nz, row_matrix_start,
                        row_matrix_index, row_matrix_value);
}

HighsStatus HighsSimplexInterface::getRows(const int num_set_entries,
                                           const int* row_set, int& num_row,
                                           double* row_lower, double* row_upper,
                                           int& num_nz, int* row_matrix_start,
                                           int* row_matrix_index,
                                           double* row_matrix_value) {
  return getRowsGeneral(false, 0, 0, true, num_set_entries, row_set, false,
                        NULL, num_row, row_lower, row_upper, num_nz,
                        row_matrix_start, row_matrix_index, row_matrix_value);
}

HighsStatus HighsSimplexInterface::getRows(const int* row_mask, int& num_row,
                                           double* row_lower, double* row_upper,
                                           int& num_nz, int* row_matrix_start,
                                           int* row_matrix_index,
                                           double* row_matrix_value) {
  return getRowsGeneral(false, 0, 0, false, 0, NULL, true, row_mask, num_row,
                        row_lower, row_upper, num_nz, row_matrix_start,
                        row_matrix_index, row_matrix_value);
}

HighsStatus HighsSimplexInterface::getRowsGeneral(
    const bool interval, const int from_row, const int to_row, const bool set,
    const int num_set_entries, const int* row_set, const bool mask,
    const int* row_mask, int& num_row, double* row_lower, double* row_upper,
    int& num_nz, int* row_matrix_start, int* row_matrix_index,
    double* row_matrix_value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  int from_k;
  int to_k;
  HighsLp& lp = highs_model_object.lp_;
  HighsOptions& options = highs_model_object.options_;
  call_status = assessIntervalSetMask(options, lp.numRow_, interval, from_row,
                                      to_row, set, num_set_entries, row_set,
                                      mask, row_mask, from_k, to_k);
  return_status =
      interpretCallStatus(call_status, return_status, "assessIntervalSetMask");
  if (return_status == HighsStatus::Error) return return_status;
  if (from_k < 0 || to_k > lp.numRow_) {
    call_status = HighsStatus::Error;
    return_status =
        interpretCallStatus(call_status, return_status, "getColsGeneral");
    return return_status;
  }
  num_row = 0;
  num_nz = 0;
  if (from_k > to_k) {
    call_status = HighsStatus::Error;
    return_status =
        interpretCallStatus(call_status, return_status, "getColsGeneral");
    return return_status;
  }
  // "Out" means not in the set to be extrated
  // "In" means in the set to be extrated
  int out_from_row;
  int out_to_row;
  int in_from_row;
  int in_to_row = -1;
  int current_set_entry = 0;
  int row_dim = lp.numRow_;
  // Set up a row mask so that entries to be got from the column-wise
  // matrix can be identified and have their correct row index.
  int* new_index = (int*)malloc(sizeof(int) * lp.numRow_);

  if (!mask) {
    out_to_row = -1;
    current_set_entry = 0;
    for (int k = from_k; k <= to_k; k++) {
      updateOutInIx(row_dim, interval, from_row, to_row, set, num_set_entries,
                    row_set, mask, row_mask, in_from_row, in_to_row,
                    out_from_row, out_to_row, current_set_entry);
      if (k == from_k) {
        // Account for any initial rows not being extracted
        for (int row = 0; row < in_from_row; row++) {
          new_index[row] = -1;
        }
      }
      for (int row = in_from_row; row <= in_to_row; row++) {
        new_index[row] = num_row;
        num_row++;
      }
      for (int row = out_from_row; row <= out_to_row; row++) {
        new_index[row] = -1;
      }
      if (out_to_row >= row_dim - 1) break;
    }
  } else {
    for (int row = 0; row < lp.numRow_; row++) {
      if (row_mask[row]) {
        new_index[row] = num_row;
        num_row++;
      } else {
        new_index[row] = -1;
      }
    }
  }

  // Bail out if no rows are to be extracted
  if (num_row == 0) {
    free(new_index);
    return HighsStatus::OK;
  }

  // Allocate an array of lengths for the row-wise matrix to be extracted
  int* row_matrix_length = (int*)malloc(sizeof(int) * num_row);

  for (int row = 0; row < lp.numRow_; row++) {
    int new_row = new_index[row];
    if (new_row >= 0) {
      assert(new_row < num_row);
      if (row_lower != NULL) row_lower[new_row] = lp.rowLower_[row];
      if (row_upper != NULL) row_upper[new_row] = lp.rowUpper_[row];
      row_matrix_length[new_row] = 0;
    }
  }
  // Identify the lengths of the rows in the row-wise matrix to be extracted
  for (int col = 0; col < lp.numCol_; col++) {
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      int new_row = new_index[row];
      if (new_row >= 0) row_matrix_length[new_row]++;
    }
  }

  if (row_matrix_start == NULL) {
    // If the matrix start vector is null then don't get values of
    // indices, otherwise both are meaningless
    if (row_matrix_index != NULL || row_matrix_value != NULL) {
      HighsLogMessage(highs_model_object.options_.logfile,
                      HighsMessageType::ERROR,
                      "Cannot supply meaningful row matrix indices/values with "
                      "null starts");
      free(new_index);
      free(row_matrix_length);
      return HighsStatus::Error;
    }
  } else {
    row_matrix_start[0] = 0;
    for (int row = 0; row < num_row - 1; row++) {
      row_matrix_start[row + 1] =
          row_matrix_start[row] + row_matrix_length[row];
    }

    // Fill the row-wise matrix with indices and values
    for (int col = 0; col < lp.numCol_; col++) {
      for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
        int row = lp.Aindex_[el];
        int new_row = new_index[row];
        if (new_row >= 0) {
          int row_el = row_matrix_start[new_row];
          if (row_matrix_index != NULL) row_matrix_index[row_el] = col;
          if (row_matrix_value != NULL)
            row_matrix_value[row_el] = lp.Avalue_[el];
          row_matrix_start[new_row]++;
        }
      }
    }
    // Restore the starts of the row-wise matrix and count the number of
    // nonzeros in it
    num_nz = 0;
    row_matrix_start[0] = 0;
    for (int row = 0; row < num_row - 1; row++) {
      row_matrix_start[row + 1] =
          row_matrix_start[row] + row_matrix_length[row];
      num_nz += row_matrix_length[row];
    }
    num_nz += row_matrix_length[num_row - 1];
  }
  free(new_index);
  free(row_matrix_length);
  return HighsStatus::OK;
}

// Change a single coefficient in the matrix
HighsStatus HighsSimplexInterface::changeCoefficient(const int Xrow,
                                                     const int Xcol,
                                                     const double XnewValue) {
#ifdef HiGHSDEV
  printf("Called changeCoeff(Xrow=%d, Xcol=%d, XnewValue=%g)\n", Xrow, Xcol,
         XnewValue);
#endif
  HighsLp& lp = highs_model_object.lp_;
  if (Xrow < 0 || Xrow > lp.numRow_) return HighsStatus::Error;
  if (Xcol < 0 || Xcol > lp.numCol_) return HighsStatus::Error;
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n",
  //  Xrow, Xcol, XnewValue);
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  bool valid_simplex_lp = simplex_lp_status.valid;
#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no matrix or scaling
  if (!valid_simplex_lp) {
    assert(!simplex_lp_status.has_matrix_col_wise);
    //    assert(!apply_row_scaling);
  }
#endif
  changeLpMatrixCoefficient(lp, Xrow, Xcol, XnewValue);
  if (valid_simplex_lp) {
    HighsLp& simplex_lp = highs_model_object.simplex_lp_;
    HighsScale& scale = highs_model_object.scale_;
    double scaledXnewValue = XnewValue * scale.row_[Xrow] * scale.col_[Xcol];
    changeLpMatrixCoefficient(simplex_lp, Xrow, Xcol, scaledXnewValue);
  }
  // simplex_lp.reportLp();
  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if it's a new row
  highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
  highs_model_object.unscaled_model_status_ =
      highs_model_object.scaled_model_status_;
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_ROWS);
  //  simplex_lp.reportLp();
  return HighsStatus::OK;
}

void HighsSimplexInterface::shiftObjectiveValue(const double Xshift) {
  printf(
      "Where is shiftObjectiveValue required - so I can interpret what's "
      "required\n");
  // Update the LP objective value with the shift
  highs_model_object.simplex_info_.dual_objective_value += Xshift;
  // Update the LP offset with the shift
  highs_model_object.lp_.offset_ += Xshift;
  if (highs_model_object.simplex_lp_status_.valid) {
    // Update the simplex LP offset with the shift
    highs_model_object.simplex_lp_.offset_ += Xshift;
  }
}

HighsStatus HighsSimplexInterface::changeObjectiveSense(const ObjSense Xsense) {
  HighsLp& lp = highs_model_object.lp_;
  if ((Xsense == ObjSense::MINIMIZE) != (lp.sense_ == ObjSense::MINIMIZE)) {
    // Flip the LP objective sense
    lp.sense_ = Xsense;
  }
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.valid) {
    HighsLp& simplex_lp = highs_model_object.simplex_lp_;
    if ((Xsense == ObjSense::MINIMIZE) !=
        (simplex_lp.sense_ == ObjSense::MINIMIZE)) {
      // Flip the objective sense
      simplex_lp.sense_ = Xsense;
      highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
      highs_model_object.unscaled_model_status_ =
          highs_model_object.scaled_model_status_;
    }
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::changeCosts(int from_col, int to_col,
                                               const double* usr_col_cost) {
  return changeCostsGeneral(true, from_col, to_col, false, 0, NULL, false, NULL,
                            usr_col_cost);
}

HighsStatus HighsSimplexInterface::changeCosts(int num_set_entries,
                                               const int* col_set,
                                               const double* usr_col_cost) {
  return changeCostsGeneral(false, 0, 0, true, num_set_entries, col_set, false,
                            NULL, usr_col_cost);
}

HighsStatus HighsSimplexInterface::changeCosts(const int* col_mask,
                                               const double* usr_col_cost) {
  return changeCostsGeneral(false, 0, 0, false, 0, NULL, true, col_mask,
                            usr_col_cost);
}

HighsStatus HighsSimplexInterface::changeCostsGeneral(
    bool interval, int from_col, int to_col, bool set, int num_set_entries,
    const int* col_set, bool mask, const int* col_mask,
    const double* usr_col_cost) {
  bool null_data = false;
  if (usr_col_cost == NULL) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::ERROR,
                    "User-supplied column costs are NULL");
    null_data = true;
  }
  if (null_data) return HighsStatus::Error;
  int* use_set;
  double* use_cost;
  if (set) {
    // Changing the costs for a set of columns, so ensure that the
    // set and data are in ascending order
    use_set = (int*)malloc(sizeof(int) * num_set_entries);
    use_cost = (double*)malloc(sizeof(double) * num_set_entries);
    sortSetData(num_set_entries, col_set, usr_col_cost, NULL, NULL, use_set,
                use_cost, NULL, NULL);
  } else {
    use_set = (int*)col_set;
    use_cost = (double*)usr_col_cost;
  }
  HighsStatus call_status = changeLpCosts(
      highs_model_object.options_, highs_model_object.lp_, interval, from_col,
      to_col, set, num_set_entries, use_set, mask, col_mask, use_cost,
      highs_model_object.options_.infinite_cost);
  if (call_status == HighsStatus::Error) return HighsStatus::Error;
  // Deduce the consequences of new costs
  highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
  highs_model_object.unscaled_model_status_ =
      highs_model_object.scaled_model_status_;
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                        LpAction::NEW_COSTS);
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::changeColBounds(
    int from_col, int to_col, const double* usr_col_lower,
    const double* usr_col_upper) {
  return changeColBoundsGeneral(true, from_col, to_col, false, 0, NULL, false,
                                NULL, usr_col_lower, usr_col_upper);
}

HighsStatus HighsSimplexInterface::changeColBounds(
    int num_set_entries, const int* col_set, const double* usr_col_lower,
    const double* usr_col_upper) {
  return changeColBoundsGeneral(false, 0, 0, true, num_set_entries, col_set,
                                false, NULL, usr_col_lower, usr_col_upper);
}

HighsStatus HighsSimplexInterface::changeColBounds(
    const int* col_mask, const double* usr_col_lower,
    const double* usr_col_upper) {
  return changeColBoundsGeneral(false, 0, 0, false, 0, NULL, true, col_mask,
                                usr_col_lower, usr_col_upper);
}

HighsStatus HighsSimplexInterface::changeColBoundsGeneral(
    bool interval, int from_col, int to_col, bool set, int num_set_entries,
    const int* col_set, bool mask, const int* col_mask,
    const double* usr_col_lower, const double* usr_col_upper) {
  bool null_data = false;
  if (usr_col_lower == NULL) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::ERROR,
                    "User-supplied column lower bounds are NULL");
    null_data = true;
  }
  if (usr_col_upper == NULL) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::ERROR,
                    "User-supplied column upper bounds are NULL");
    null_data = true;
  }
  if (null_data) return HighsStatus::Error;
  int* use_set;
  double* use_lower;
  double* use_upper;
  if (set) {
    // Changing the bounds for a set of columns, so ensure that the set
    // and data are in ascending order
    use_set = (int*)malloc(sizeof(int) * num_set_entries);
    use_lower = (double*)malloc(sizeof(double) * num_set_entries);
    use_upper = (double*)malloc(sizeof(double) * num_set_entries);
    sortSetData(num_set_entries, col_set, usr_col_lower, usr_col_upper, NULL,
                use_set, use_lower, use_upper, NULL);
  } else {
    use_set = (int*)col_set;
    use_lower = (double*)usr_col_lower;
    use_upper = (double*)usr_col_upper;
  }
  HighsStatus call_status = changeLpColBounds(
      highs_model_object.options_, highs_model_object.lp_, interval, from_col,
      to_col, set, num_set_entries, use_set, mask, col_mask, use_lower,
      use_upper, highs_model_object.options_.infinite_bound);
  if (call_status == HighsStatus::Error) return HighsStatus::Error;

  if (highs_model_object.simplex_lp_status_.valid) {
    // Also change the simplex LP's column bounds
    assert(highs_model_object.lp_.numCol_ ==
           highs_model_object.simplex_lp_.numCol_);
    assert(highs_model_object.lp_.numRow_ ==
           highs_model_object.simplex_lp_.numRow_);

    call_status = changeLpColBounds(
        highs_model_object.options_, highs_model_object.simplex_lp_, interval,
        from_col, to_col, set, num_set_entries, use_set, mask, col_mask,
        use_lower, use_upper, highs_model_object.options_.infinite_bound);
    if (call_status == HighsStatus::Error) return HighsStatus::Error;
    if (highs_model_object.scale_.is_scaled_) {
      scaleLpColBounds(highs_model_object.options_,
                       highs_model_object.simplex_lp_,
                       highs_model_object.scale_.col_, interval, from_col,
                       to_col, set, num_set_entries, use_set, mask, col_mask);
    }
    // Deduce the consequences of new col bounds
    highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
    highs_model_object.unscaled_model_status_ =
        highs_model_object.scaled_model_status_;
    updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                          LpAction::NEW_BOUNDS);
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::changeRowBounds(
    int from_row, int to_row, const double* usr_row_lower,
    const double* usr_row_upper) {
  return changeRowBoundsGeneral(true, from_row, to_row, false, 0, NULL, false,
                                NULL, usr_row_lower, usr_row_upper);
}

HighsStatus HighsSimplexInterface::changeRowBounds(
    int num_set_entries, const int* row_set, const double* usr_row_lower,
    const double* usr_row_upper) {
  return changeRowBoundsGeneral(false, 0, 0, true, num_set_entries, row_set,
                                false, NULL, usr_row_lower, usr_row_upper);
}

HighsStatus HighsSimplexInterface::changeRowBounds(
    const int* row_mask, const double* usr_row_lower,
    const double* usr_row_upper) {
  return changeRowBoundsGeneral(false, 0, 0, false, 0, NULL, true, row_mask,
                                usr_row_lower, usr_row_upper);
}

HighsStatus HighsSimplexInterface::changeRowBoundsGeneral(
    bool interval, int from_row, int to_row, bool set, int num_set_entries,
    const int* row_set, bool mask, const int* row_mask,
    const double* usr_row_lower, const double* usr_row_upper) {
  bool null_data = false;
  if (usr_row_lower == NULL) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::ERROR,
                    "User-supplied row lower bounds are NULL");
    null_data = true;
  }
  if (usr_row_upper == NULL) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::ERROR,
                    "User-supplied row upper bounds are NULL");
    null_data = true;
  }
  if (null_data) return HighsStatus::Error;
  int* use_set;
  double* use_lower;
  double* use_upper;
  if (set) {
    // Changing the bounds for a set of rows, so ensure that the set
    // and data are in ascending order
    use_set = (int*)malloc(sizeof(int) * num_set_entries);
    use_lower = (double*)malloc(sizeof(double) * num_set_entries);
    use_upper = (double*)malloc(sizeof(double) * num_set_entries);
    sortSetData(num_set_entries, row_set, usr_row_lower, usr_row_upper, NULL,
                use_set, use_lower, use_upper, NULL);
  } else {
    use_set = (int*)row_set;
    use_lower = (double*)usr_row_lower;
    use_upper = (double*)usr_row_upper;
  }
  HighsStatus call_status = changeLpRowBounds(
      highs_model_object.options_, highs_model_object.lp_, interval, from_row,
      to_row, set, num_set_entries, use_set, mask, row_mask, use_lower,
      use_upper, highs_model_object.options_.infinite_bound);
  if (call_status == HighsStatus::Error) return HighsStatus::Error;
  if (highs_model_object.simplex_lp_status_.valid) {
    // Also change the simplex LP's column bounds
    assert(highs_model_object.lp_.numCol_ ==
           highs_model_object.simplex_lp_.numCol_);
    assert(highs_model_object.lp_.numRow_ ==
           highs_model_object.simplex_lp_.numRow_);
    call_status = changeLpRowBounds(
        highs_model_object.options_, highs_model_object.simplex_lp_, interval,
        from_row, to_row, set, num_set_entries, use_set, mask, row_mask,
        use_lower, use_upper, highs_model_object.options_.infinite_bound);
    if (call_status == HighsStatus::Error) return HighsStatus::Error;
    if (highs_model_object.scale_.is_scaled_) {
      scaleLpRowBounds(highs_model_object.options_,
                       highs_model_object.simplex_lp_,
                       highs_model_object.scale_.row_, interval, from_row,
                       to_row, set, num_set_entries, row_set, mask, row_mask);
    }
    // Deduce the consequences of new row bounds
    highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
    highs_model_object.unscaled_model_status_ =
        highs_model_object.scaled_model_status_;
    updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                          LpAction::NEW_BOUNDS);
  }
  return HighsStatus::OK;
}

// Solve (transposed) system involving the basis matrix

HighsStatus HighsSimplexInterface::basisSolve(const vector<double>& rhs,
                                              double* solution_vector,
                                              int* solution_num_nz,
                                              int* solution_indices,
                                              bool transpose) {
  HVector solve_vector;
  int numRow = highs_model_object.simplex_lp_.numRow_;
  int numCol = highs_model_object.simplex_lp_.numCol_;
  HighsScale& scale = highs_model_object.scale_;
  // Set up solve vector with suitably scaled RHS
  solve_vector.setup(numRow);
  solve_vector.clear();
  int rhs_num_nz = 0;
  if (transpose) {
    for (int row = 0; row < numRow; row++) {
      if (rhs[row]) {
        solve_vector.index[rhs_num_nz++] = row;
        double rhs_value = rhs[row];
        int col = highs_model_object.simplex_basis_.basicIndex_[row];
        if (col < numCol) {
          //	  printf("RHS row %2d: col %2d: scale rhs[row] = %11.4g by
          //%11.4g\n", row, col, rhs_value, scale.col_[col]);
          rhs_value *= scale.col_[col];
        } else {
          double scale_value = scale.row_[col - numCol];
          //	  printf("RHS row %2d: row %2d: scale rhs[row] = %11.4g by
          //%11.4g\n", row, col - numCol, rhs_value, scale_value);
          rhs_value /= scale_value;
        }
        solve_vector.array[row] = rhs_value;
      }
    }
  } else {
    for (int row = 0; row < numRow; row++) {
      if (rhs[row]) {
        solve_vector.index[rhs_num_nz++] = row;
        //	printf("RHS row %2d: scale rhs[row] = %11.4g by scale.row_[row]
        //= %11.4g\n", row, rhs[row], scale.row_[row]);
        solve_vector.array[row] = rhs[row] * scale.row_[row];
      }
    }
  }
  solve_vector.count = rhs_num_nz;
  //  printf("RHS has %d nonzeros\n", rhs_num_nz);
  //
  // Note that solve_vector.count is just used to determine whether
  // hyper-sparse solves should be used. The indices of the nonzeros
  // in the solution are always accumulated. There's no switch (such
  // as setting solve_vector.count = numRow+1) to not do this.
  //
  // Get hist_dsty from analysis during simplex solve.
  double hist_dsty = 1;
  if (transpose) {
    highs_model_object.factor_.btran(solve_vector, hist_dsty);
  } else {
    highs_model_object.factor_.ftran(solve_vector, hist_dsty);
  }
  //  printf("After solve: solve_vector.count = %d\n", solve_vector.count);
  // Extract the solution
  if (solution_indices == NULL) {
    // Nonzeros in the solution not required
    if (solve_vector.count > numRow) {
      // Solution nonzeros not known
      for (int row = 0; row < numRow; row++) {
        solution_vector[row] = solve_vector.array[row];
        //	printf("Solution vector[%2d] = solve_vector.array[row] =
        //%11.4g\n", row, solution_vector[row]);
      }
    } else {
      // Solution nonzeros are known
      for (int row = 0; row < numRow; row++) solution_vector[row] = 0;
      for (int ix = 0; ix < solve_vector.count; ix++) {
        int row = solve_vector.index[ix];
        solution_vector[row] = solve_vector.array[row];
        //	printf("Solution vector[%2d] = solve_vector.array[row] = %11.4g
        // from index %2d\n", row, solution_vector[row], ix);
      }
    }
  } else {
    // Nonzeros in the solution are required
    if (solve_vector.count > numRow) {
      // Solution nonzeros not known
      solution_num_nz = 0;
      for (int row = 0; row < numRow; row++) {
        solution_vector[row] = 0;
        if (solve_vector.array[row]) {
          solution_vector[row] = solve_vector.array[row];
          solution_indices[*solution_num_nz++] = row;
          //	  printf("Solution vector[%2d] = solve_vector.array[row] =
          //%11.4g from index %2d\n", row, solution_vector[row],
          // solution_num_nz-1);
        }
      }
    } else {
      // Solution nonzeros are known
      for (int row = 0; row < numRow; row++) solution_vector[row] = 0;
      for (int ix = 0; ix < solve_vector.count; ix++) {
        int row = solve_vector.index[ix];
        solution_vector[row] = solve_vector.array[row];
        solution_indices[ix] = row;
        //	printf("Solution vector[%2d] = solve_vector.array[row] = %11.4g
        // from index %2d\n", row, solution_vector[row], ix);
      }
      *solution_num_nz = solve_vector.count;
    }
  }
  // Scale the solution
  if (transpose) {
    if (solve_vector.count > numRow) {
      // Solution nonzeros not known
      for (int row = 0; row < numRow; row++) {
        double scale_value = scale.row_[row];
        solution_vector[row] *= scale_value;
        //	printf("Row %2d so scale by %11.4g to give %11.4g\n", row,
        // scale_value, solution_vector[row]);
      }
    } else {
      for (int ix = 0; ix < solve_vector.count; ix++) {
        int row = solve_vector.index[ix];
        double scale_value = scale.row_[row];
        solution_vector[row] *= scale_value;
        //	printf("Row %2d so scale by %11.4g to give %11.4g\n", row,
        // scale_value, solution_vector[row]);
      }
    }
  } else {
    if (solve_vector.count > numRow) {
      // Solution nonzeros not known
      for (int row = 0; row < numRow; row++) {
        int col = highs_model_object.simplex_basis_.basicIndex_[row];
        if (col < numCol) {
          solution_vector[row] *= scale.col_[col];
          //	  printf("Col %2d so scale by %11.4g to give %11.4g\n", col,
          // scale.col_[col], solution_vector[row]);
        } else {
          double scale_value = scale.row_[col - numCol];
          solution_vector[row] /= scale_value;
          //	  printf("Row %2d so scale by %11.4g to give %11.4g\n", col -
          // numCol, scale_value, solution_vector[row]);
        }
      }
    } else {
      for (int ix = 0; ix < solve_vector.count; ix++) {
        int row = solve_vector.index[ix];
        int col = highs_model_object.simplex_basis_.basicIndex_[row];
        if (col < numCol) {
          solution_vector[row] *= scale.col_[col];
          //	  printf("Col %2d so scale by %11.4g to give %11.4g\n", col,
          // scale.col_[col], solution_vector[row]);
        } else {
          double scale_value = scale.row_[col - numCol];
          solution_vector[row] /= scale_value;
          //	  printf("Row %2d so scale by %11.4g to give %11.4g\n", col -
          // numCol, scale_value, solution_vector[row]);
        }
      }
    }
  }
  //  for (int row = 0; row < numRow; row++) printf("Solution vector[%2d] =
  //  %11.4g\n", row, solution_vector[row]);
  return HighsStatus::OK;
}

#ifdef HiGHSDEV
void HighsSimplexInterface::change_update_method(int updateMethod) {
  highs_model_object.factor_.change(updateMethod);
}
#endif

// Utilities to convert model basic/nonbasic status to/from SCIP-like status
int HighsSimplexInterface::convertBaseStatToHighsBasis(const int* cstat,
                                                       const int* rstat) {
  HighsBasis& basis = highs_model_object.basis_;
  HighsLp& lp = highs_model_object.lp_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  //  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;

  int numBasic = 0;
  bool error_found = false;
  basis.valid_ = false;
  for (int col = 0; col < lp.numCol_; col++) {
    if (cstat[col] == (int)HighsBasisStatus::BASIC) {
      basis.col_status[col] = HighsBasisStatus::BASIC;
      numBasic++;
      continue;
    }
    if (cstat[col] == (int)HighsBasisStatus::LOWER) {
      // Supplied basis has this column nonbasic at its lower bound: check that
      // the lower bound is finite
      error_found = highs_isInfinity(-lp.colLower_[col]);
      basis.col_status[col] = HighsBasisStatus::LOWER;
    } else if (cstat[col] == (int)HighsBasisStatus::UPPER) {
      // Supplied basis has this column nonbasic at its upper bound: check that
      // the upper bound is finite
      error_found = highs_isInfinity(lp.colUpper_[col]);
      basis.col_status[col] = HighsBasisStatus::UPPER;
    } else if (cstat[col] == (int)HighsBasisStatus::ZERO) {
      // Supplied basis has this column nonbasic at zero so free: check that
      // neither bound is finite
      error_found = !highs_isInfinity(-lp.colLower_[col]) ||
                    !highs_isInfinity(lp.colUpper_[col]);
      basis.col_status[col] = HighsBasisStatus::UPPER;
    } else {
      error_found = true;
    }
    if (error_found) {
#ifdef HiGHSDEV
      printf("Invalid basis status: col=%d, cstat=%d, lower=%g, upper=%g\n",
             col, cstat[col], lp.colLower_[col], lp.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  for (int row = 0; row < lp.numRow_; row++) {
    if (rstat[row] == (int)HighsBasisStatus::BASIC) {
      basis.row_status[row] = HighsBasisStatus::BASIC;
      numBasic++;
      continue;
    }
    if (rstat[row] == (int)HighsBasisStatus::LOWER) {
      // Supplied basis has this row nonbasic at its lower bound: check that the
      // lower bound is finite
      error_found = highs_isInfinity(-lp.rowLower_[row]);
      basis.row_status[row] = HighsBasisStatus::LOWER;
    } else if (rstat[row] == (int)HighsBasisStatus::UPPER) {
      // Supplied basis has this row nonbasic at its upper bound: check that the
      // upper bound is finite
      error_found = highs_isInfinity(lp.rowUpper_[row]);
      basis.row_status[row] = HighsBasisStatus::UPPER;
    } else if (rstat[row] == (int)HighsBasisStatus::ZERO) {
      // Supplied basis has this row nonbasic at zero so free: check that
      // neither bound is finite
      error_found = !highs_isInfinity(-lp.rowLower_[row]) ||
                    !highs_isInfinity(lp.rowUpper_[row]);
      basis.row_status[row] = HighsBasisStatus::UPPER;
    } else {
      error_found = true;
    }
    if (error_found) {
#ifdef HiGHSDEV
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], lp.rowLower_[row], lp.rowUpper_[row]);
#endif
      return -(row + 1);
    }
  }
  assert(numBasic = lp.numRow_);
  basis.valid_ = true;
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_BASIS);
  return 0;
}

int HighsSimplexInterface::convertHighsBasisToBaseStat(int* cstat, int* rstat) {
  HighsBasis& basis = highs_model_object.basis_;
  HighsLp& lp = highs_model_object.lp_;
  if (cstat != NULL) {
    for (int col = 0; col < lp.numCol_; col++)
      cstat[col] = (int)basis.col_status[col];
  }
  printf("NB SCIP has row bounds [-u, -l]\n");
  if (rstat != NULL) {
    for (int row = 0; row < lp.numRow_; row++)
      rstat[row] = (int)basis.row_status[row];
  }
  return 0;
}

void HighsSimplexInterface::convertSimplexToHighsBasis() {
  HighsBasis& basis = highs_model_object.basis_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsLp& lp = highs_model_object.lp_;
  //  HighsSimplexLpStatus &simplex_lp_status =
  //  highs_model_object.simplex_lp_status_; HighsSimplexInfo &simplex_info =
  //  highs_model_object.simplex_info_;
  basis.col_status.resize(lp.numCol_);
  basis.row_status.resize(lp.numRow_);

  assert(highs_model_object.simplex_lp_status_.has_basis);
  bool permuted = highs_model_object.simplex_lp_status_.is_permuted;
  int* numColPermutation =
      &highs_model_object.simplex_info_.numColPermutation_[0];
  // numColPermutation[iCol] is the true column in column iCol
  const bool optimal_basis =
      highs_model_object.scaled_model_status_ == HighsModelStatus::OPTIMAL;
  bool error_found = false;
  basis.valid_ = false;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    int simplex_var = iCol;
    int lp_col = iCol;
    const double lower = lp.colLower_[lp_col];
    const double upper = lp.colUpper_[lp_col];
    HighsBasisStatus basis_status;
    if (permuted) lp_col = numColPermutation[iCol];
    if (!simplex_basis.nonbasicFlag_[simplex_var]) {
      basis_status = HighsBasisStatus::BASIC;
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_UP) {
      // Nonbasic and free to move up so should be OK to give status LOWER
      if (optimal_basis) {
#ifdef HiGHSDEV
        // Check that the lower bound isn't infinite
        error_found = highs_isInfinity(-lower);
#endif
        basis_status = HighsBasisStatus::LOWER;
      } else {
        basis_status = checkedVarHighsNonbasicStatus(HighsBasisStatus::LOWER,
                                                     lower, upper);
      }
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_DN) {
      // Nonbasic and free to move down so should be OK to give status UPPER
      if (optimal_basis) {
#ifdef HiGHSDEV
        // Check that the upper bound isn't infinite
        error_found = highs_isInfinity(upper);
#endif
        basis_status = HighsBasisStatus::UPPER;
      } else {
        basis_status = checkedVarHighsNonbasicStatus(HighsBasisStatus::UPPER,
                                                     lower, upper);
      }
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_ZE) {
      // Column is either fixed or free, depending on the bounds
      if (lower == upper) {
        // Equal bounds so should be OK to give status LOWER
        if (optimal_basis) {
#ifdef HiGHSDEV
          // Check that the lower bound isn't infinite
          error_found = highs_isInfinity(-lower);
#endif
          basis_status = HighsBasisStatus::LOWER;
        } else {
          basis_status = checkedVarHighsNonbasicStatus(HighsBasisStatus::LOWER,
                                                       lower, upper);
        }
      } else {
        // Unequal bounds so should be OK to give status ZERO
        if (optimal_basis) {
#ifdef HiGHSDEV
          // Check that neither of the bounds is finite
          error_found = !highs_isInfinity(-lower) || !highs_isInfinity(upper);
#endif
          basis_status = HighsBasisStatus::ZERO;
        } else {
          basis_status = checkedVarHighsNonbasicStatus(HighsBasisStatus::ZERO,
                                                       lower, upper);
        }
      }
    } else {
      error_found = true;
    }
    if (error_found) {
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: col=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          lp_col, simplex_basis.nonbasicFlag_[simplex_var],
          simplex_basis.nonbasicMove_[simplex_var], lower, upper);
#endif
      assert(!error_found);
      return;
    } else {
      basis.col_status[lp_col] = basis_status;
    }
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    int simplex_var = lp.numCol_ + iRow;
    int lp_row = iRow;
    const double lower = lp.rowLower_[lp_row];
    const double upper = lp.rowUpper_[lp_row];
    HighsBasisStatus basis_status;
    if (!simplex_basis.nonbasicFlag_[simplex_var]) {
      basis_status = HighsBasisStatus::BASIC;
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_UP) {
      // Nonbasic and free to move up so should be OK to give status UPPER -
      // since simplex row bounds are flipped and negated
      if (optimal_basis) {
#ifdef HiGHSDEV
        // Check that the upper bound isn't infinite
        error_found = highs_isInfinity(upper);
#endif
        basis_status = HighsBasisStatus::UPPER;
      } else {
        basis_status = checkedVarHighsNonbasicStatus(HighsBasisStatus::UPPER,
                                                     lower, upper);
      }
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_DN) {
      // Nonbasic and free to move down so should be OK to give status
      // LOWER - since simplex row bounds are flipped and negated
      if (optimal_basis) {
#ifdef HiGHSDEV
        // Check that the lower bound isn't infinite
        error_found = highs_isInfinity(-lower);
#endif
        basis_status = HighsBasisStatus::LOWER;
      } else {
        basis_status = checkedVarHighsNonbasicStatus(HighsBasisStatus::LOWER,
                                                     lower, upper);
      }
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_ZE) {
      // Row is either fixed or free, depending on the bounds
      if (lower == upper) {
        // Equal bounds so should be OK to give status LOWER
        if (optimal_basis) {
#ifdef HiGHSDEV
          // Check that the lower bound isn't infinite
          error_found = highs_isInfinity(-lower);
#endif
          basis_status = HighsBasisStatus::LOWER;
        } else {
          basis_status = checkedVarHighsNonbasicStatus(HighsBasisStatus::LOWER,
                                                       lower, upper);
        }
      } else {
        // Unequal bounds so should be OK to give status ZERO
        if (optimal_basis) {
#ifdef HiGHSDEV
          // Check that neither of the bounds is finite
          error_found = !highs_isInfinity(-lower) || !highs_isInfinity(upper);
#endif
          basis_status = HighsBasisStatus::ZERO;
        } else {
          basis_status = checkedVarHighsNonbasicStatus(HighsBasisStatus::ZERO,
                                                       lower, upper);
        }
      }
    } else {
      error_found = true;
    }
    if (error_found) {
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: row=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          lp_row, simplex_basis.nonbasicFlag_[simplex_var],
          simplex_basis.nonbasicMove_[simplex_var], lower, upper);
#endif
      assert(!error_found);
      return;
    } else {
      basis.row_status[lp_row] = basis_status;
    }
  }
  basis.valid_ = true;
}

void HighsSimplexInterface::convertHighsToSimplexBasis() {
  HighsBasis& basis = highs_model_object.basis_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsLp& lp = highs_model_object.lp_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  bool error_found = false;
  bool permuted = highs_model_object.simplex_lp_status_.is_permuted;
  int* numColPermutation =
      &highs_model_object.simplex_info_.numColPermutation_[0];
  // numColPermutation[iCol] is the true column in column iCol
  int num_basic = 0;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    int simplex_var = iCol;
    int lp_col = iCol;
    if (permuted) lp_col = numColPermutation[iCol];
    if (basis.col_status[lp_col] == HighsBasisStatus::BASIC) {
      // Basic
      simplex_basis.nonbasicFlag_[simplex_var] = NONBASIC_FLAG_FALSE;
      simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
      simplex_basis.basicIndex_[num_basic] = simplex_var;
      num_basic++;
    } else {
      // Nonbasic
      simplex_basis.nonbasicFlag_[simplex_var] = NONBASIC_FLAG_TRUE;
      if (basis.col_status[lp_col] == HighsBasisStatus::LOWER) {
        // HighsBasisStatus::LOWER includes fixed variables
#ifdef HiGHSDEV
        // Check that the lower bound isn't infinite
        error_found = highs_isInfinity(-lp.colLower_[lp_col]);
#endif
        if (lp.colLower_[lp_col] == lp.colUpper_[lp_col]) {
          // Equal bounds so indicate that the column can't move
#ifdef HiGHSDEV
          // Check that the upper bound isn't infinite
          error_found = highs_isInfinity(lp.colUpper_[lp_col]);
#endif
          simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
        } else {
          // unequal bounds so indicate that the column can only move up
          simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_UP;
        }
      } else if (basis.col_status[lp_col] == HighsBasisStatus::UPPER) {
        // HighsBasisStatus::UPPER includes only variables at their upper bound
#ifdef HiGHSDEV
        // Check that the upper bound isn't infinite
        error_found = highs_isInfinity(lp.colUpper_[lp_col]);
#endif
        simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_DN;
      } else if (basis.col_status[lp_col] == HighsBasisStatus::ZERO) {
        // HighsBasisStatus::ZERO implies a free variable
#ifdef HiGHSDEV
        // Check that neither bound is finite
        error_found = !highs_isInfinity(-lp.colLower_[lp_col]) ||
                      !highs_isInfinity(lp.colUpper_[lp_col]);
#endif
        simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
      } else {
        error_found = true;
      }
    }
#ifdef HiGHSDEV
    if (error_found)
      printf(
          "Invalid basis status: col=%d, basis.col_status=%d, lower=%g, "
          "upper=%g\n",
          lp_col, (int)basis.col_status[lp_col], lp.colLower_[lp_col],
          lp.colUpper_[lp_col]);
#endif
    assert(!error_found);
    if (error_found) return;
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    int simplex_var = lp.numCol_ + iRow;
    int lp_row = iRow;
    if (basis.row_status[lp_row] == HighsBasisStatus::BASIC) {
      // Basic
      simplex_basis.nonbasicFlag_[simplex_var] = NONBASIC_FLAG_FALSE;
      simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
      simplex_basis.basicIndex_[num_basic] = simplex_var;
      num_basic++;
    } else {
      // Nonbasic
      simplex_basis.nonbasicFlag_[simplex_var] = NONBASIC_FLAG_TRUE;
      if (basis.row_status[lp_row] == HighsBasisStatus::LOWER) {
        // HighsBasisStatus::LOWER includes fixed variables
#ifdef HiGHSDEV
        // Check that the lower bound isn't infinite
        error_found = highs_isInfinity(-lp.rowLower_[lp_row]);
#endif
        if (lp.rowLower_[lp_row] == lp.rowUpper_[lp_row]) {
          // Equal bounds so indicate that the row can't move
#ifdef HiGHSDEV
          // Check that the upper bound isn't infinite
          error_found = highs_isInfinity(lp.rowUpper_[lp_row]);
#endif
          simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
        } else {
          // Unequal bounds so indicate that the row can only move
          // down - since simplex row bounds are flipped and negated
          simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_DN;
        }
      } else if (basis.row_status[lp_row] == HighsBasisStatus::UPPER) {
        // HighsBasisStatus::UPPER includes only variables at their upper bound
#ifdef HiGHSDEV
        // Check that the upper bound isn't infinite
        error_found = highs_isInfinity(lp.rowUpper_[lp_row]);
#endif
        // Upper bounded so indicate that the row can only move
        // up - since simplex row bounds are flipped and negated
        simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_UP;
      } else if (basis.row_status[lp_row] == HighsBasisStatus::ZERO) {
        // HighsBasisStatus::ZERO implies a free variable
#ifdef HiGHSDEV
        // Check that neither bound is finite
        error_found = !highs_isInfinity(-lp.rowLower_[lp_row]) ||
                      !highs_isInfinity(lp.rowUpper_[lp_row]);
#endif
        simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
      } else {
        error_found = true;
      }
    }
#ifdef HiGHSDEV
    if (error_found)
      printf(
          "Invalid basis status: row=%d, basis.row_status=%d, lower=%g, "
          "upper=%g\n",
          lp_row, (int)basis.row_status[lp_row], lp.rowLower_[lp_row],
          lp.rowUpper_[lp_row]);
#endif
    assert(!error_found);
    if (error_found) return;
  }
  assert(num_basic = lp.numRow_);
  //  populate_work_arrays(highs_model_object); // Why might this have been done
  //  here?
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_BASIS);
  simplex_lp_status.has_basis = true;
}

void HighsSimplexInterface::convertSimplexToHighsSolution() {
  HighsSolution& solution = highs_model_object.solution_;
  HighsScale& scale = highs_model_object.scale_;
  SimplexBasis& basis = highs_model_object.simplex_basis_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;

  // Take primal solution
  vector<double> value = simplex_info.workValue_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
    value[basis.basicIndex_[iRow]] = simplex_info.baseValue_[iRow];
  // Take dual solution
  vector<double> dual = simplex_info.workDual_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
    dual[basis.basicIndex_[iRow]] = 0;
  // Scale back
  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    value[iCol] *= scale.col_[iCol];
    dual[iCol] /= (scale.col_[iCol] / scale.cost_);
  }
  for (int iRow = 0, iTot = simplex_lp.numCol_; iRow < simplex_lp.numRow_;
       iRow++, iTot++) {
    value[iTot] /= scale.row_[iRow];
    dual[iTot] *= (scale.row_[iRow] * scale.cost_);
  }

  // Now we can get the solution
  solution.col_value.resize(simplex_lp.numCol_);
  solution.col_dual.resize(simplex_lp.numCol_);
  solution.row_value.resize(simplex_lp.numRow_);
  solution.row_dual.resize(simplex_lp.numRow_);

  if (highs_model_object.simplex_lp_status_.is_permuted) {
    const int* numColPermutation =
        &highs_model_object.simplex_info_.numColPermutation_[0];
    // numColPermutation[i] is the true column in column i
    for (int i = 0; i < simplex_lp.numCol_; i++) {
      int iCol = numColPermutation[i];
      solution.col_value[iCol] = value[i];
      solution.col_dual[iCol] = (int)simplex_lp.sense_ * dual[i];
    }
  } else {
    for (int i = 0; i < simplex_lp.numCol_; i++) {
      int iCol = i;
      solution.col_value[iCol] = value[i];
      solution.col_dual[iCol] = (int)simplex_lp.sense_ * dual[i];
    }
  }
  int row_dual_sign = 1;  //-1;
  for (int i = 0; i < simplex_lp.numRow_; i++) {
    solution.row_value[i] = -value[i + simplex_lp.numCol_];
    solution.row_dual[i] =
        row_dual_sign * (int)simplex_lp.sense_ * dual[i + simplex_lp.numCol_];
  }
}

int HighsSimplexInterface::get_basic_indices(int* bind) {
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = simplex_basis.basicIndex_[row];
    if (var >= simplex_lp.numCol_)
      bind[row] = -(1 + var - simplex_lp.numCol_);
    else
      bind[row] = var;
  }
  return 0;
}

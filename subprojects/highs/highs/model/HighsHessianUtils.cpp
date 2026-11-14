/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsHessianUtils.cpp
 * @brief
 */
#include "model/HighsHessianUtils.h"

#include <algorithm>
#include <cmath>

#include "lp_data/HighsModelUtils.h"
#include "util/HighsMatrixUtils.h"
#include "util/HighsSort.h"

using std::fabs;

HighsStatus assessHessian(HighsHessian& hessian, const HighsOptions& options,
                          const std::string* col_names) {
  if (hessian.isOracle()) return HighsStatus::kOk;
  const HighsLogOptions& log_options = options.log_options;
  if (!hessian.legalDimensions(log_options)) return HighsStatus::kError;
  HighsStatus return_status = hessian.trimDimensions(log_options);

  // If the Hessian has no columns there is nothing left to test
  if (hessian.dim_ == 0) {
    hessian.clear();
    return return_status;
  }
  // Assess the Hessian matrix
  //
  // The start of column 0 must be zero.
  if (hessian.start_[0]) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Hessian has nonzero value (%" HIGHSINT_FORMAT
                 ") for the start of column 0\n",
                 hessian.start_[0]);
    return HighsStatus::kError;
  }
  // Assess Q, summing duplicates, but deferring the assessment of
  // values (other than those which are identically zero)
  bool sum_duplicates = true;
  HighsStatus call_status =
      assessMatrix(log_options, "HighsHessian", hessian.name_, hessian.dim_,
                   hessian.dim_, hessian.start_, hessian.index_, hessian.value_,
                   0, kHighsInf, sum_duplicates, col_names, col_names);
  return_status = interpretCallStatus(log_options, call_status, return_status,
                                      "assessMatrix");
  if (return_status == HighsStatus::kError) return return_status;
  // Transform the Hessian to pure column-wise lower triangle format
  call_status = normaliseHessian(options, hessian);
  return_status = interpretCallStatus(log_options, call_status, return_status,
                                      "normaliseHessian");
  if (return_status == HighsStatus::kError) return return_status;
  // Assess values in Q
  sum_duplicates = false;
  call_status =
      assessMatrix(log_options, "HighsHessian", hessian.name_, hessian.dim_,
                   hessian.dim_, hessian.start_, hessian.index_, hessian.value_,
                   options.small_matrix_value, options.large_matrix_value,
                   sum_duplicates, col_names, col_names);
  return_status = interpretCallStatus(log_options, call_status, return_status,
                                      "assessMatrix");
  if (return_status == HighsStatus::kError) return return_status;

  HighsInt hessian_num_nz = hessian.numNz();
  // If the Hessian has nonzeros, complete its diagonal with explicit
  // zeros if necessary
  if (hessian_num_nz) {
    completeHessianDiagonal(options, hessian);
    hessian_num_nz = hessian.numNz();
  }
  // If entries have been removed from the matrix, resize the index
  // and value vectors
  if ((HighsInt)hessian.index_.size() > hessian_num_nz)
    hessian.index_.resize(hessian_num_nz);
  if ((HighsInt)hessian.value_.size() > hessian_num_nz)
    hessian.value_.resize(hessian_num_nz);

  //  if (return_status != HighsStatus::kError) return_status =
  //  HighsStatus::kOk;
  if (return_status != HighsStatus::kOk)
    highsLogDev(log_options, HighsLogType::kInfo,
                "assessHessian returns HighsStatus = %s\n",
                highsStatusToString(return_status).c_str());
  return return_status;
}

void completeHessianDiagonal(const HighsOptions& options,
                             HighsHessian& hessian) {
  assert(!hessian.isOracle());
  // Count the number of missing diagonal entries
  HighsInt num_missing_diagonal_entries = 0;
  const HighsInt dim = hessian.dim_;
  const HighsInt num_nz = hessian.numNz();
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    HighsInt iEl = hessian.start_[iCol];
    if (iEl < num_nz) {
      if (hessian.index_[iEl] != iCol) num_missing_diagonal_entries++;
    } else {
      num_missing_diagonal_entries++;
    }
  }
  if (num_missing_diagonal_entries > 0)
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "Hessian has dimension %d and %d nonzeros: inserting %d zeros "
                "onto the diagonal\n",
                int(dim), int(num_nz), int(num_missing_diagonal_entries));
  assert(num_missing_diagonal_entries >= dim - num_nz);
  if (!num_missing_diagonal_entries) return;
  // There are missing diagonal entries to be inserted as explicit zeros
  const HighsInt new_num_nz = hessian.numNz() + num_missing_diagonal_entries;
  HighsInt to_iEl = new_num_nz;
  hessian.index_.resize(new_num_nz);
  hessian.value_.resize(new_num_nz);
  HighsInt next_start = hessian.numNz();
  hessian.start_[dim] = to_iEl;
  HighsInt num_missing_diagonal_entries_added = 0;
  for (HighsInt iCol = dim - 1; iCol >= 0; iCol--) {
    // Shift the entries that are sure to be off-diagonal
    for (HighsInt iEl = next_start - 1; iEl > hessian.start_[iCol]; iEl--) {
      assert(hessian.index_[iEl] != iCol);
      to_iEl--;
      hessian.index_[to_iEl] = hessian.index_[iEl];
      hessian.value_[to_iEl] = hessian.value_[iEl];
    }
    // Now consider any first entry. If there is none, or if it's not
    // the diagonal, then there is no diagonal entry for this column
    bool no_diagonal_entry;
    if (hessian.start_[iCol] < next_start) {
      const HighsInt iEl = hessian.start_[iCol];
      // Copy the first entry
      to_iEl--;
      hessian.index_[to_iEl] = hessian.index_[iEl];
      hessian.value_[to_iEl] = hessian.value_[iEl];
      // If the first entry isn't the diagonal, then there is no
      // diagonal entry for this column
      no_diagonal_entry = hessian.index_[iEl] != iCol;
    } else {
      no_diagonal_entry = true;
    }
    if (no_diagonal_entry) {
      // There is no diagonal entry, so have insert an explicit zero
      to_iEl--;
      hessian.index_[to_iEl] = iCol;
      hessian.value_[to_iEl] = 0;
      num_missing_diagonal_entries_added++;
      assert(num_missing_diagonal_entries_added <=
             num_missing_diagonal_entries);
    }
    next_start = hessian.start_[iCol];
    hessian.start_[iCol] = to_iEl;
  }
  assert(to_iEl == 0);
}

bool okHessianDiagonal(const HighsOptions& options, HighsHessian& hessian,
                       const ObjSense sense) {
  double min_diagonal_value = kHighsInf;
  double max_diagonal_value = -kHighsInf;
  const HighsInt dim = hessian.isOracle() ? hessian.oracle_.dim_ : hessian.dim_;
  const HighsInt sense_sign = (HighsInt)sense;
  HighsInt num_illegal_diagonal_value = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    double diagonal_value = sense_sign * hessian.diag(iCol);
    min_diagonal_value = std::min(diagonal_value, min_diagonal_value);
    max_diagonal_value = std::max(diagonal_value, max_diagonal_value);
    // Diagonal entries signed by sense must be non-negative
    if (diagonal_value < 0) num_illegal_diagonal_value++;
  }

  const bool certainly_not_positive_semidefinite =
      num_illegal_diagonal_value > 0;
  if (certainly_not_positive_semidefinite) {
    if (sense == ObjSense::kMinimize) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Hessian has %" HIGHSINT_FORMAT
                   " diagonal entries in [%g, 0) so is not positive "
                   "semidefinite for minimization\n",
                   num_illegal_diagonal_value, min_diagonal_value);
    } else {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Hessian has %" HIGHSINT_FORMAT
                   " diagonal entries in (0, %g] so is not negative "
                   "semidefinite for maximization\n",
                   num_illegal_diagonal_value, -min_diagonal_value);
    }
  }
  return !certainly_not_positive_semidefinite;
}

HighsStatus extractTriangularHessian(const HighsOptions& options,
                                     HighsHessian& hessian) {
  assert(!hessian.isOracle());
  if (hessian.isOracle()) return HighsStatus::kError;
  // Viewing the Hessian column-wise, remove any entries in the strict
  // upper triangle
  HighsStatus return_status = HighsStatus::kOk;
  const HighsInt dim = hessian.dim_;
  HighsInt nnz = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    const HighsInt nnz0 = nnz;
    for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol + 1];
         iEl++) {
      HighsInt iRow = hessian.index_[iEl];
      if (iRow < iCol) continue;
      hessian.index_[nnz] = iRow;
      hessian.value_[nnz] = hessian.value_[iEl];
      if (iRow == iCol && nnz > nnz0) {
        // Diagonal entry is not first in column so swap it in
        hessian.index_[nnz] = hessian.index_[nnz0];
        hessian.value_[nnz] = hessian.value_[nnz0];
        hessian.index_[nnz0] = iRow;
        hessian.value_[nnz0] = hessian.value_[iEl];
      }
      nnz++;
    }
    hessian.start_[iCol] = nnz0;
  }
  const HighsInt num_ignored_nz = hessian.start_[dim] - nnz;
  assert(num_ignored_nz >= 0);
  if (num_ignored_nz) {
    if (hessian.format_ == HessianFormat::kTriangular) {
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "Ignored %" HIGHSINT_FORMAT
                   " entries of Hessian in opposite triangle\n",
                   num_ignored_nz);
      return_status = HighsStatus::kWarning;
    }
    hessian.start_[dim] = nnz;
  }
  assert(hessian.start_[dim] == nnz);
  hessian.format_ = HessianFormat::kTriangular;
  return return_status;
}

void triangularToSquareHessian(const HighsHessian& hessian,
                               vector<HighsInt>& start, vector<HighsInt>& index,
                               vector<double>& value) {
  assert(!hessian.isOracle());
  if (hessian.isOracle()) return;
  const HighsInt dim = hessian.dim_;
  if (dim <= 0) {
    start.assign(1, 0);
    return;
  }
  assert(hessian.format_ == HessianFormat::kTriangular);
  // Cannot use the following, since it puts the diagonal entry first
  // in each column and, strangely, the active set QP solver seems to
  // need the entries column-by-column in row order
  /*
  HighsHessian square_hessian = hessian.toSquare();
  std::vector<HighsInt> nw_start = square_hessian.start_;
  std::vector<HighsInt> nw_index = square_hessian.index_;
  std::vector<double> nw_value = square_hessian.value_;
  */
  const HighsInt nnz = hessian.start_[dim];
  const HighsInt square_nnz = nnz + (nnz - dim);
  start.resize(dim + 1);
  index.resize(square_nnz);
  value.resize(square_nnz);
  vector<HighsInt> length;
  length.assign(dim, 0);
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    HighsInt iRow = hessian.index_[hessian.start_[iCol]];
    assert(iRow == iCol);
    length[iCol]++;
    for (HighsInt iEl = hessian.start_[iCol] + 1;
         iEl < hessian.start_[iCol + 1]; iEl++) {
      HighsInt iRow = hessian.index_[iEl];
      assert(iRow > iCol);
      length[iRow]++;
      length[iCol]++;
    }
  }
  start[0] = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++)
    start[iCol + 1] = start[iCol] + length[iCol];
  assert(square_nnz == start[dim]);
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    HighsInt iEl = hessian.start_[iCol];
    HighsInt iRow = hessian.index_[iEl];
    HighsInt toEl = start[iCol];
    index[toEl] = iRow;
    value[toEl] = hessian.value_[iEl];
    start[iCol]++;
    for (HighsInt iEl = hessian.start_[iCol] + 1;
         iEl < hessian.start_[iCol + 1]; iEl++) {
      HighsInt iRow = hessian.index_[iEl];
      HighsInt toEl = start[iRow];
      index[toEl] = iCol;
      value[toEl] = hessian.value_[iEl];
      start[iRow]++;
      toEl = start[iCol];
      index[toEl] = iRow;
      value[toEl] = hessian.value_[iEl];
      start[iCol]++;
    }
  }
  start[0] = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++)
    start[iCol + 1] = start[iCol] + length[iCol];
}

HighsStatus normaliseHessian(const HighsOptions& options,
                             HighsHessian& hessian) {
  assert(!hessian.isOracle());
  if (hessian.isOracle()) return HighsStatus::kError;
  HighsInt dim = hessian.dim_;
  const bool triangular = hessian.format_ == HessianFormat::kTriangular;
  const bool square = !triangular;
  assert((hessian.format_ == HessianFormat::kSquare) == square);
  HighsSparseMatrix upper;
  std::vector<HighsInt> upper_length;
  upper_length.assign(dim, 0);
  HighsInt num_upper_nz = 0;
  HighsInt num_lower_nz = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol + 1];
         iEl++) {
      HighsInt iRow = hessian.index_[iEl];
      if (iRow < iCol) {
        upper_length[iRow]++;
        num_upper_nz++;
      } else if (iRow > iCol) {
        num_lower_nz++;
      }
    }
  }
  if (triangular && num_upper_nz == 0) return HighsStatus::kOk;
  // If square, do quick check for asymmetry
  if (square && num_upper_nz != num_lower_nz) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Hessian has %d / %d lower / upper triangular entries so is "
                 "not symmetric\n",
                 int(num_lower_nz), int(num_upper_nz));
    return HighsStatus::kError;
  }
  if (num_upper_nz > 0) {
    // Form the HighsSparseMatrix of strict upper triangular entries
    num_lower_nz = 0;
    // Define the upper starts, and use upper_length to store the
    // element number for the next nonzero in each row
    for (HighsInt iRow = 0; iRow < dim; iRow++) {
      upper.start_.push_back(upper.start_[iRow] + upper_length[iRow]);
      upper_length[iRow] = upper.start_[iRow];
    }
    upper.index_.resize(num_upper_nz);
    upper.value_.resize(num_upper_nz);
    for (HighsInt iCol = 0; iCol < dim; iCol++) {
      // Loop through the Hessian column, moving its upper triangular
      // nonzeros to upper, so have to shift the lower nonzeros
      HighsInt from_el = hessian.start_[iCol];
      // Define the new start for this lower column
      hessian.start_[iCol] = num_lower_nz;
      for (HighsInt iEl = from_el; iEl < hessian.start_[iCol + 1]; iEl++) {
        HighsInt iRow = hessian.index_[iEl];
        if (iRow < iCol) {
          // Store the nonzero in the upper row
          upper.index_[upper_length[iRow]] = iCol;
          upper.value_[upper_length[iRow]] = hessian.value_[iEl];
          upper_length[iRow]++;
        } else {
          // Shift the nonzero in the lower column
          hessian.index_[num_lower_nz] = iRow;
          hessian.value_[num_lower_nz] = hessian.value_[iEl];
          num_lower_nz++;
        }
      }
    }
    hessian.start_[dim] = num_lower_nz;
    // Check that upper_length has reached the start of the next row
    for (HighsInt iRow = 0; iRow < dim; iRow++)
      assert(upper_length[iRow] == upper.start_[iRow + 1]);
    assert(upper_length[dim - 1] == num_upper_nz);
  } else {
    upper.start_.resize(dim + 1, 0);
  }
  // Have to work from a copy of the Hessian if it is triangular and
  // there are upper triangular entries to insert
  HighsHessian hessian_copy;
  if (triangular) hessian_copy = hessian;
  const HighsHessian& from_hessian = triangular ? hessian_copy : hessian;

  // Determine the lower (upper) off-diagonal entries in each column
  // (row) of the Hessian
  std::vector<double> lower_on_below_diagonal;
  std::vector<double> upper_off_diagonal;
  lower_on_below_diagonal.assign(dim, 0);
  upper_off_diagonal.assign(dim, 0);
  // Should be no duplicates
  HighsInt debug_num_duplicate = 0;
  HighsInt num_summation = 0;
  HighsInt num_upper_triangle = 0;
  HighsInt num_hessian_el = 0;
  HighsInt num_illegal_asymmetry = 0;
  double min_illegal_asymmetry = kHighsInf;
  double max_illegal_asymmetry = 0;
  HighsInt num_ok_asymmetry = 0;
  double min_ok_asymmetry = kHighsInf;
  double max_ok_asymmetry = 0;

  //  const bool expensive_2821_check = true;
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    for (HighsInt iEl = from_hessian.start_[iCol];
         iEl < from_hessian.start_[iCol + 1]; iEl++) {
      HighsInt iRow = from_hessian.index_[iEl];
      assert(iRow >= iCol);
      if (lower_on_below_diagonal[iRow]) debug_num_duplicate++;
      lower_on_below_diagonal[iRow] = from_hessian.value_[iEl];
    }
    // Now look at the corresponding row of the upper triangular
    // matrix - deliberately referring to it as iCol and picking out
    // its indices as iRow
    for (HighsInt iEl = upper.start_[iCol]; iEl < upper.start_[iCol + 1];
         iEl++) {
      HighsInt iRow = upper.index_[iEl];
      assert(iRow > iCol);
      if (upper_off_diagonal[iRow]) debug_num_duplicate++;
      upper_off_diagonal[iRow] = upper.value_[iEl];
      if (square) {
        // When square, ensure that the upper triangular value matches
        // the corresponding lower triangular value to within the
        // toleranace used by JuMP
        const double asymmetry =
            std::fabs(upper_off_diagonal[iRow] - lower_on_below_diagonal[iRow]);
        if (asymmetry > kSquareHessianAsymmetryTolerance) {
          num_illegal_asymmetry++;
          min_illegal_asymmetry = std::min(asymmetry, min_illegal_asymmetry);
          max_illegal_asymmetry = std::max(asymmetry, max_illegal_asymmetry);
        } else if (asymmetry) {
          num_ok_asymmetry++;
          min_ok_asymmetry = std::min(asymmetry, min_ok_asymmetry);
          max_ok_asymmetry = std::max(asymmetry, max_ok_asymmetry);
          double average =
              (upper_off_diagonal[iRow] + lower_on_below_diagonal[iRow]) * 0.5;
          upper_off_diagonal[iRow] = average;
          lower_on_below_diagonal[iRow] = average;
        }
        // Don't zero the upper off diagonal entry, so that it's
        // possible to check that nonzeros in the lower off diagonal
        // match the upper off diagonal entry
      } else {
        // When triangular, add the upper triangular value to the
        // lower triangular value, and zero the upper triangular value
        // so that it doesn't generate a duplicate
        assert(triangular);
        num_upper_triangle++;
        if (lower_on_below_diagonal[iRow]) num_summation++;
        lower_on_below_diagonal[iRow] += upper_off_diagonal[iRow];
        upper_off_diagonal[iRow] = 0;
      }
    }
    // Now gather the nonzeros, zeroing lower_on_below_diagonal and
    // upper_off_diagonal
    HighsInt from_el = from_hessian.start_[iCol];
    hessian.start_[iCol] = num_hessian_el;
    // Assign any nonzero diagonal entry
    HighsInt iRow = iCol;
    if (lower_on_below_diagonal[iRow]) {
      hessian.index_[num_hessian_el] = iRow;
      hessian.value_[num_hessian_el] = lower_on_below_diagonal[iRow];
      num_hessian_el++;
      lower_on_below_diagonal[iRow] = 0;
    }
    // Assign nonzeros below the diagonal
    for (HighsInt iEl = from_el; iEl < from_hessian.start_[iCol + 1]; iEl++) {
      HighsInt iRow = from_hessian.index_[iEl];
      if (lower_on_below_diagonal[iRow]) {
        hessian.index_[num_hessian_el] = iRow;
        hessian.value_[num_hessian_el] = lower_on_below_diagonal[iRow];
        num_hessian_el++;
        lower_on_below_diagonal[iRow] = 0;
      }
    }
    // Assign nonzeros above the diagonal
    for (HighsInt iEl = upper.start_[iCol]; iEl < upper.start_[iCol + 1];
         iEl++) {
      HighsInt iRow = upper.index_[iEl];
      if (triangular) {
        assert(!upper_off_diagonal[iRow]);
      }
      // Look for nonzeros in lower_on_below_diagonal created by
      // adding an entry that was in the upper triangle and added in
      if (lower_on_below_diagonal[iRow]) {
        hessian.index_[num_hessian_el] = iRow;
        hessian.value_[num_hessian_el] = lower_on_below_diagonal[iRow];
        num_hessian_el++;
        lower_on_below_diagonal[iRow] = 0;
      }
      // Any upper entry retained for checking symmetry for square
      // Hessians needs to be zeroed
      upper_off_diagonal[iRow] = 0;
    }
    /*
    if (expensive_2821_check) {
      // Check that lower_on_below_diagonal and upper_off_diagonal
      // have been zeroed
      for (HighsInt iRow = 0; iRow < dim; iRow++) {
        assert(!lower_on_below_diagonal[iRow]);
        assert(!upper_off_diagonal[iRow]);
      }
    }
    */
  }  // Loop iCol = 0; iCol < dim; iCol++
  hessian.start_[dim] = num_hessian_el;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.index_.resize(num_hessian_el);
  hessian.value_.resize(num_hessian_el);

  assert(!debug_num_duplicate);

  bool warning_found = false;
  bool error_found = false;
  if (num_ok_asymmetry) {
    assert(square);
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "Square Hessian contains %d non-symmetr%s in [%.2g, %.2g] "
                 "within tolerance of %.1g\n",
                 int(num_ok_asymmetry), num_ok_asymmetry == 1 ? "y" : "ies",
                 min_ok_asymmetry, max_ok_asymmetry,
                 kSquareHessianAsymmetryTolerance);
  }
  if (num_illegal_asymmetry) {
    assert(square);
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Square Hessian contains %d non-symmetr%s in [%.2g, %.2g] "
                 "exceeding tolerance of %.1g\n",
                 int(num_illegal_asymmetry),
                 num_illegal_asymmetry == 1 ? "y" : "ies",
                 min_illegal_asymmetry, max_illegal_asymmetry,
                 kSquareHessianAsymmetryTolerance);
    error_found = true;
  }
  if (num_upper_triangle) {
    assert(triangular);
    highsLogUser(
        options.log_options, HighsLogType::kWarning,
        "Triangular Hessian contains %d entr%s in upper triangle: added to "
        "lower triangle, requiring %d non-trivial summation%s\n",
        int(num_upper_triangle), num_upper_triangle == 1 ? "y" : "ies",
        int(num_summation), num_summation == 1 ? "" : "s");
    assert(num_summation <= num_upper_triangle);
    warning_found = true;
  }
  HighsStatus return_status = HighsStatus::kOk;
  if (error_found)
    return_status = HighsStatus::kError;
  else if (warning_found)
    return_status = HighsStatus::kWarning;
  return return_status;
}

void completeHessian(const HighsInt full_dim, HighsHessian& hessian) {
  assert(!hessian.isOracle());
  if (hessian.isOracle()) return;
  // Ensure that any non-zero Hessian of dimension less than the
  // number of columns in the model is completed with explicit zero
  // diagonal entries
  assert(hessian.dim_ <= full_dim);
  if (hessian.dim_ == full_dim) return;
  HighsInt nnz = hessian.numNz();
  hessian.exactResize();
  for (HighsInt iCol = hessian.dim_; iCol < full_dim; iCol++) {
    hessian.index_.push_back(iCol);
    hessian.value_.push_back(0);
    nnz++;
    hessian.start_.push_back(nnz);
  }
  hessian.dim_ = full_dim;
  assert(HighsInt(hessian.start_.size()) == hessian.dim_ + 1);
}

void reportHessian(const HighsLogOptions& log_options, const HighsInt dim,
                   const HighsInt num_nz, const HighsInt* start,
                   const HighsInt* index, const double* value) {
  if (dim <= 0) return;
  highsLogUser(log_options, HighsLogType::kInfo,
               "Hessian Index              Value\n");
  for (HighsInt col = 0; col < dim; col++) {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "    %8" HIGHSINT_FORMAT " Start   %10" HIGHSINT_FORMAT "\n",
                 col, start[col]);
    HighsInt to_el = (col < dim - 1 ? start[col + 1] : num_nz);
    for (HighsInt el = start[col]; el < to_el; el++)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "          %8" HIGHSINT_FORMAT " %12g\n", index[el],
                   value[el]);
  }
  highsLogUser(log_options, HighsLogType::kInfo,
               "             Start   %10" HIGHSINT_FORMAT "\n", num_nz);
}

void userScaleHessian(HighsHessian& hessian, HighsUserScaleData& data,
                      const bool apply) {
  assert(!hessian.isOracle());
  if (hessian.isOracle()) return;
  data.num_infinite_hessian_values = 0;
  if (!hessian.dim_) return;
  const HighsInt user_objective_scale = data.user_objective_scale;
  const HighsInt user_bound_scale = data.user_bound_scale;
  if (!user_objective_scale && !user_bound_scale) return;
  // If variable bounds are scaled by bound_scale_value, then linear
  // term in objective is scaled by bound_scale_value, but Hessian
  // term is scaled by bound_scale_value**2, so have to scale down the
  // Hessian values so linear and Hessian terms are both scaled by the
  // same constant value
  double objective_scale_value = std::pow(2, user_objective_scale);
  double bound_scale_value = std::pow(2, -user_bound_scale);
  for (HighsInt iEl = 0; iEl < hessian.start_[hessian.dim_]; iEl++) {
    double value =
        hessian.value_[iEl] * objective_scale_value * bound_scale_value;
    if (std::abs(value) > data.infinite_cost)
      data.num_infinite_hessian_values++;
    if (apply) hessian.value_[iEl] = value;
  }
}

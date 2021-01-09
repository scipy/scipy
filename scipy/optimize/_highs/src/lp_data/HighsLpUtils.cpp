/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsLpUtils.h"

#include <algorithm>
#include <cassert>

#include "HConfig.h"
#include "io/Filereader.h"
#include "io/HMPSIO.h"
#include "io/HighsIO.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsSort.h"
#include "util/HighsTimer.h"

HighsStatus assessLp(HighsLp& lp, const HighsOptions& options) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Assess the LP dimensions and vector sizes, returning on error
  call_status = assessLpDimensions(options, lp);
  return_status =
      interpretCallStatus(call_status, return_status, "assessLpDimensions");
  if (return_status == HighsStatus::Error) return return_status;

  // If the LP has no columns there is nothing left to test
  // NB assessLpDimensions returns HighsStatus::Error if lp.numCol_ < 0
  if (lp.numCol_ == 0) return HighsStatus::OK;

  // From here, any LP has lp.numCol_ > 0 and lp.Astart_[lp.numCol_] exists (as
  // the number of nonzeros)
  assert(lp.numCol_ > 0);

  // Assess the LP column costs
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = 0;
  index_collection.to_ = lp.numCol_ - 1;
  call_status = assessCosts(options, 0, index_collection, lp.colCost_,
                            options.infinite_cost);
  return_status =
      interpretCallStatus(call_status, return_status, "assessCosts");
  if (return_status == HighsStatus::Error) return return_status;
  // Assess the LP column bounds
  call_status = assessBounds(options, "Col", 0, index_collection, lp.colLower_,
                             lp.colUpper_, options.infinite_bound);
  return_status =
      interpretCallStatus(call_status, return_status, "assessBounds");
  if (return_status == HighsStatus::Error) return return_status;
  if (lp.numRow_) {
    // Assess the LP row bounds
    index_collection.dimension_ = lp.numRow_;
    index_collection.is_interval_ = true;
    index_collection.from_ = 0;
    index_collection.to_ = lp.numRow_ - 1;
    call_status =
        assessBounds(options, "Row", 0, index_collection, lp.rowLower_,
                     lp.rowUpper_, options.infinite_bound);
    return_status =
        interpretCallStatus(call_status, return_status, "assessBounds");
    if (return_status == HighsStatus::Error) return return_status;
    // Assess the LP matrix
    //
    // The start of column 0 must be zero. It's possible to make
    // everything consistent with the start being positive but, in
    // particular, this means that the number of nonzeros is no longer
    // lp.Astart_[lp.numCol_], a real banana skin to avoid just for
    // the few wierdos who want lp.Astart_[0] to be positive. Hence
    // it's not permitted!
    if (lp.Astart_[0]) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "LP has nonzero value (%d) for the start of column 0\n",
                      lp.Astart_[0]);
      return HighsStatus::Error;
    }
    call_status = assessMatrix(
        options, lp.numRow_, lp.numCol_, lp.Astart_, lp.Aindex_, lp.Avalue_,
        options.small_matrix_value, options.large_matrix_value);
    return_status =
        interpretCallStatus(call_status, return_status, "assessMatrix");
    if (return_status == HighsStatus::Error) return return_status;
    int lp_num_nz = lp.Astart_[lp.numCol_];
    // If entries have been removed from the matrix, resize the index
    // and value vectors to prevent bug in presolve
    if ((int)lp.Aindex_.size() > lp_num_nz) lp.Aindex_.resize(lp_num_nz);
    if ((int)lp.Avalue_.size() > lp_num_nz) lp.Avalue_.resize(lp_num_nz);
  }
  if (return_status == HighsStatus::Error)
    return_status = HighsStatus::Error;
  else
    return_status = HighsStatus::OK;
#ifdef HiGHSDEV
  HighsLogMessage(options.logfile, HighsMessageType::INFO,
                  "assess_lp returns HighsStatus = %s",
                  HighsStatusToString(return_status).c_str());
#endif
  return return_status;
}

HighsStatus assessLpDimensions(const HighsOptions& options, const HighsLp& lp) {
  HighsStatus return_status = HighsStatus::OK;

  // Use error_found to track whether an error has been found in multiple tests
  bool error_found = false;

  // Don't expect the matrix_start_size to be legal if there are no columns
  bool check_matrix_start_size = lp.numCol_ > 0;

  // Assess column-related dimensions
  bool legal_num_col = lp.numCol_ >= 0;
  if (!legal_num_col) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "LP has illegal number of cols = %d\n", lp.numCol_);
    error_found = true;
  } else {
    // Check the size of the column vectors
    int col_cost_size = lp.colCost_.size();
    int col_lower_size = lp.colLower_.size();
    int col_upper_size = lp.colUpper_.size();
    int matrix_start_size = lp.Astart_.size();
    bool legal_col_cost_size = col_cost_size >= lp.numCol_;
    bool legal_col_lower_size = col_lower_size >= lp.numCol_;
    bool legal_col_upper_size = col_lower_size >= lp.numCol_;

    if (!legal_col_cost_size) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "LP has illegal colCost size = %d < %d\n", col_cost_size,
                      lp.numCol_);
      error_found = true;
    }
    if (!legal_col_lower_size) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "LP has illegal colLower size = %d < %d\n",
                      col_lower_size, lp.numCol_);
      error_found = true;
    }
    if (!legal_col_upper_size) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "LP has illegal colUpper size = %d < %d\n",
                      col_upper_size, lp.numCol_);
      error_found = true;
    }
    if (check_matrix_start_size) {
      bool legal_matrix_start_size = matrix_start_size >= lp.numCol_ + 1;
      if (!legal_matrix_start_size) {
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "LP has illegal Astart size = %d < %d\n",
                        matrix_start_size, lp.numCol_ + 1);
        error_found = true;
      }
    }
  }

  // Assess row-related dimensions
  bool legal_num_row = lp.numRow_ >= 0;
  if (!legal_num_row) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "LP has illegal number of rows = %d\n", lp.numRow_);
    error_found = true;
  } else {
    int row_lower_size = lp.rowLower_.size();
    int row_upper_size = lp.rowUpper_.size();
    bool legal_row_lower_size = row_lower_size >= lp.numRow_;
    bool legal_row_upper_size = row_lower_size >= lp.numRow_;
    if (!legal_row_lower_size) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "LP has illegal rowLower size = %d < %d\n",
                      row_lower_size, lp.numRow_);
      error_found = true;
    }
    if (!legal_row_upper_size) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "LP has illegal rowUpper size = %d < %d\n",
                      row_upper_size, lp.numRow_);
      error_found = true;
    }
  }

  // Assess matrix-related dimensions
  if (check_matrix_start_size) {
    int lp_num_nz = lp.Astart_[lp.numCol_];
    bool legal_num_nz = lp_num_nz >= 0;
    if (!legal_num_nz) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "LP has illegal number of nonzeros = %d\n", lp_num_nz);
      error_found = true;
    } else {
      int matrix_index_size = lp.Aindex_.size();
      int matrix_value_size = lp.Avalue_.size();
      bool legal_matrix_index_size = matrix_index_size >= lp_num_nz;
      bool legal_matrix_value_size = matrix_value_size >= lp_num_nz;
      if (!legal_matrix_index_size) {
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "LP has illegal Aindex size = %d < %d\n",
                        matrix_index_size, lp_num_nz);
        error_found = true;
      }
      if (!legal_matrix_value_size) {
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "LP has illegal Avalue size = %d < %d\n",
                        matrix_value_size, lp_num_nz);
        error_found = true;
      }
    }
  }
  if (error_found)
    return_status = HighsStatus::Error;
  else
    return_status = HighsStatus::OK;

  return return_status;
}

HighsStatus assessCosts(const HighsOptions& options, const int ml_col_os,
                        const HighsIndexCollection& index_collection,
                        vector<double>& cost, const double infinite_cost) {
  HighsStatus return_status = HighsStatus::OK;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return return_status;

  return_status = HighsStatus::OK;
  bool error_found = false;
  // Work through the data to be assessed.
  //
  // Loop is k \in [from_k...to_k) covering the entries in the
  // interval, set or mask to be considered.
  //
  // For an interval or mask, these values of k are the columns to be
  // considered in a local sense, as well as the entries in the
  // cost data to be assessed
  //
  // For a set, these values of k are the indices in the set, from
  // which the columns to be considered in a local sense are
  // drawn. The entries in the cost data to be assessed correspond
  // to the values of k
  //
  // Adding the value of ml_col_os to local_col yields the value of
  // ml_col, being the column in a global (whole-model) sense. This is
  // necessary when assessing the costs of columns being added to a
  // model, since they are specified using an interval
  // [0...num_new_col) which must be offset by the current number of
  // columns in the model.
  //
  int local_col;
  int data_col;
  int ml_col;
  for (int k = from_k; k < to_k + 1; k++) {
    if (index_collection.is_interval_ || index_collection.is_mask_) {
      local_col = k;
      data_col = k;
    } else {
      local_col = index_collection.set_[k];
      data_col = k;
    }
    ml_col = ml_col_os + local_col;
    if (index_collection.is_mask_ && !index_collection.mask_[local_col])
      continue;
    double abs_cost = fabs(cost[data_col]);
    bool legal_cost = abs_cost < infinite_cost;
    if (!legal_cost) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Col  %12d has |cost| of %12g >= %12g", ml_col, abs_cost,
                      infinite_cost);
      error_found = true;
    }
  }
  if (error_found)
    return_status = HighsStatus::Error;
  else
    return_status = HighsStatus::OK;

  return return_status;
}

HighsStatus assessBounds(const HighsOptions& options, const char* type,
                         const int ml_ix_os,
                         const HighsIndexCollection& index_collection,
                         vector<double>& lower, vector<double>& upper,
                         const double infinite_bound) {
  HighsStatus return_status = HighsStatus::OK;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::OK;

  return_status = HighsStatus::OK;
  bool error_found = false;
  bool warning_found = false;
  // Work through the data to be assessed.
  //
  // Loop is k \in [from_k...to_k) covering the entries in the
  // interval, set or mask to be considered.
  //
  // For an interval or mask, these values of k are the row/column
  // indices to be considered in a local sense, as well as the entries
  // in the lower and upper bound data to be assessed
  //
  // For a set, these values of k are the indices in the set, from
  // which the indices to be considered in a local sense are
  // drawn. The entries in the lower and
  // upper bound data to be assessed correspond to the values of
  // k.
  //
  // Adding the value of ml_ix_os to local_ix yields the value of
  // ml_ix, being the index in a global (whole-model) sense. This is
  // necessary when assessing the bounds of rows/columns being added
  // to a model, since they are specified using an interval
  // [0...num_new_row/col) which must be offset by the current number
  // of rows/columns (generically indices) in the model.
  //
  int num_infinite_lower_bound = 0;
  int num_infinite_upper_bound = 0;
  int local_ix;
  int data_ix;
  int ml_ix;
  for (int k = from_k; k < to_k + 1; k++) {
    if (index_collection.is_interval_ || index_collection.is_mask_) {
      local_ix = k;
      data_ix = k;
    } else {
      local_ix = index_collection.set_[k];
      data_ix = k;
    }
    ml_ix = ml_ix_os + local_ix;
    if (index_collection.is_mask_ && !index_collection.mask_[local_ix])
      continue;

    if (!highs_isInfinity(-lower[data_ix])) {
      // Check whether a finite lower bound will be treated as -Infinity
      bool infinite_lower_bound = lower[data_ix] <= -infinite_bound;
      if (infinite_lower_bound) {
        lower[data_ix] = -HIGHS_CONST_INF;
        num_infinite_lower_bound++;
      }
    }
    if (!highs_isInfinity(upper[data_ix])) {
      // Check whether a finite upper bound will be treated as Infinity
      bool infinite_upper_bound = upper[data_ix] >= infinite_bound;
      if (infinite_upper_bound) {
        upper[data_ix] = HIGHS_CONST_INF;
        num_infinite_upper_bound++;
      }
    }
    // Check that the lower bound does not exceed the upper bound
    bool legalLowerUpperBound = lower[data_ix] <= upper[data_ix];
    if (!legalLowerUpperBound) {
      // Leave inconsistent bounds to be used to deduce infeasibility
      HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                      "%3s  %12d has inconsistent bounds [%12g, %12g]", type,
                      ml_ix, lower[data_ix], upper[data_ix]);
      warning_found = true;
    }
    // Check that the lower bound is not as much as +Infinity
    bool legalLowerBound = lower[data_ix] < infinite_bound;
    if (!legalLowerBound) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "%3s  %12d has lower bound of %12g >= %12g", type, ml_ix,
                      lower[data_ix], infinite_bound);
      error_found = true;
    }
    // Check that the upper bound is not as little as -Infinity
    bool legalUpperBound = upper[data_ix] > -infinite_bound;
    if (!legalUpperBound) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "%3s  %12d has upper bound of %12g <= %12g", type, ml_ix,
                      upper[data_ix], -infinite_bound);
      error_found = true;
    }
  }
  if (num_infinite_lower_bound) {
    HighsLogMessage(
        options.logfile, HighsMessageType::INFO,
        "%3ss:%12d lower bounds exceeding %12g are treated as -Infinity", type,
        num_infinite_lower_bound, -infinite_bound);
  }
  if (num_infinite_upper_bound) {
    HighsLogMessage(
        options.logfile, HighsMessageType::INFO,
        "%3ss:%12d upper bounds exceeding %12g are treated as +Infinity", type,
        num_infinite_upper_bound, infinite_bound);
  }

  if (error_found)
    return_status = HighsStatus::Error;
  else if (warning_found)
    return_status = HighsStatus::Warning;
  else
    return_status = HighsStatus::OK;

  return return_status;
}

HighsStatus assessMatrix(const HighsOptions& options, const int vec_dim,
                         const int num_vec, vector<int>& Astart,
                         vector<int>& Aindex, vector<double>& Avalue,
                         const double small_matrix_value,
                         const double large_matrix_value) {
  int num_nz = Astart[num_vec];
  if (num_nz > 0 && vec_dim <= 0) return HighsStatus::Error;
  if (num_nz <= 0) return HighsStatus::OK;

  HighsStatus return_status = HighsStatus::OK;
  bool error_found = false;
  bool warning_found = false;

  // Return a error if the first start is not zero
  if (Astart[0]) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Matrix starts do not begin with 0");
    return HighsStatus::Error;
  }
  // Assess the starts
  // Set up previous_start for a fictitious previous empty packed vector
  int previous_start = Astart[0];
  for (int ix = 0; ix < num_vec; ix++) {
    int this_start = Astart[ix];
    bool this_start_too_small = this_start < previous_start;
    if (this_start_too_small) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Matrix packed vector %d has illegal start of %d < %d = "
                      "previous start",
                      ix, this_start, previous_start);
      return HighsStatus::Error;
    }
    bool this_start_too_big = this_start > num_nz;
    if (this_start_too_big) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Matrix packed vector %d has illegal start of %d > %d = "
                      "number of nonzeros",
                      ix, this_start, num_nz);
      return HighsStatus::Error;
    }
  }

  // Assess the indices and values
  // Count the number of acceptable indices/values
  int num_new_nz = 0;
  int num_small_values = 0;
  double max_small_value = 0;
  double min_small_value = HIGHS_CONST_INF;
  // Set up a zeroed vector to detect duplicate indices
  vector<int> check_vector;
  if (vec_dim > 0) check_vector.assign(vec_dim, 0);
  for (int ix = 0; ix < num_vec; ix++) {
    int from_el = Astart[ix];
    int to_el = Astart[ix + 1];
    // Account for any index-value pairs removed so far
    Astart[ix] = num_new_nz;
    for (int el = from_el; el < to_el; el++) {
      int component = Aindex[el];
      // Check that the index is non-negative
      bool legal_component = component >= 0;
      if (!legal_component) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Matrix packed vector %d, entry %d, is illegal index %d", ix, el,
            component);
        return HighsStatus::Error;
      }
      // Check that the index does not exceed the vector dimension
      legal_component = component < vec_dim;
      if (!legal_component) {
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Matrix packed vector %d, entry %d, is illegal index "
                        "%12d >= %d = vector dimension",
                        ix, el, component, vec_dim);
        return HighsStatus::Error;
      }
      // Check that the index has not already ocurred
      legal_component = check_vector[component] == 0;
      if (!legal_component) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Matrix packed vector %d, entry %d, is duplicate index %d", ix, el,
            component);
        return HighsStatus::Error;
      }
      // Indicate that the index has occurred
      check_vector[component] = 1;
      // Check that the value is not too large
      double abs_value = fabs(Avalue[el]);
      bool large_value = abs_value >= large_matrix_value;
      if (large_value) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Matrix packed vector %d, entry %d, is large value |%g| >= %g", ix,
            el, abs_value, large_matrix_value);
        return HighsStatus::Error;
      }
      bool ok_value = abs_value > small_matrix_value;
      if (!ok_value) {
        if (max_small_value < abs_value) max_small_value = abs_value;
        if (min_small_value > abs_value) min_small_value = abs_value;
        num_small_values++;
      }
      if (ok_value) {
        // Shift the index and value of the OK entry to the new
        // position in the index and value vectors, and increment
        // the new number of nonzeros
        Aindex[num_new_nz] = Aindex[el];
        Avalue[num_new_nz] = Avalue[el];
        num_new_nz++;
      } else {
        // Zero the check_vector entry since the small value
        // _hasn't_ occurred
        check_vector[component] = 0;
      }
    }
    // Zero check_vector
    for (int el = Astart[ix]; el < num_new_nz; el++)
      check_vector[Aindex[el]] = 0;
#ifdef HiGHSDEV
    // NB This is very expensive so shouldn't be true
    const bool check_check_vector = false;
    if (check_check_vector) {
      // Check zeroing of check vector
      for (int component = 0; component < vec_dim; component++) {
        if (check_vector[component]) error_found = true;
      }
      if (error_found)
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "assessMatrix: check_vector not zeroed");
    }
#endif
  }
  if (num_small_values) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Matrix packed vector contains %d |values| in [%g, %g] "
                    "less than %g: ignored",
                    num_small_values, min_small_value, max_small_value,
                    small_matrix_value);
    warning_found = true;
  }
  Astart[num_vec] = num_new_nz;
  if (error_found)
    return_status = HighsStatus::Error;
  else if (warning_found)
    return_status = HighsStatus::Warning;
  else
    return_status = HighsStatus::OK;

  return return_status;
}

HighsStatus cleanBounds(const HighsOptions& options, HighsLp& lp) {
  double max_residual = 0;
  int num_change = 0;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double residual = lp.colLower_[iCol] - lp.colUpper_[iCol];
    if (residual > options.primal_feasibility_tolerance) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Column %d has inconsistent bounds [%g, %g] (residual = "
                      "%g) after presolve ",
                      iCol, lp.colLower_[iCol], lp.colUpper_[iCol], residual);
      return HighsStatus::Error;
    } else if (residual > 0) {
      num_change++;
      max_residual = std::max(residual, max_residual);
      double mid = 0.5 * (lp.colLower_[iCol] + lp.colUpper_[iCol]);
      lp.colLower_[iCol] = mid;
      lp.colUpper_[iCol] = mid;
    }
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    double residual = lp.rowLower_[iRow] - lp.rowUpper_[iRow];
    if (residual > options.primal_feasibility_tolerance) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Row %d has inconsistent bounds [%g, %g] (residual = %g) "
                      "after presolve ",
                      iRow, lp.rowLower_[iRow], lp.rowUpper_[iRow], residual);
      return HighsStatus::Error;
    } else if (residual > 0) {
      num_change++;
      max_residual = std::max(residual, max_residual);
      double mid = 0.5 * (lp.rowLower_[iRow] + lp.rowUpper_[iRow]);
      lp.rowLower_[iRow] = mid;
      lp.rowUpper_[iRow] = mid;
    }
  }
  if (num_change) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Resolved %d inconsistent bounds (maximum residual = "
                    "%9.4g) after presolve ",
                    num_change, max_residual);
    return HighsStatus::Warning;
  }
  return HighsStatus::OK;
}

HighsStatus applyScalingToLp(const HighsOptions& options, HighsLp& lp,
                             const HighsScale& scale) {
  if (!scale.is_scaled_) return HighsStatus::OK;
  if ((int)scale.col_.size() < lp.numCol_) return HighsStatus::Error;
  if ((int)scale.row_.size() < lp.numRow_) return HighsStatus::Error;
  bool scale_error = false;
  // Set up column and row index collections for scaling
  HighsIndexCollection all_cols;
  all_cols.is_interval_ = true;
  all_cols.dimension_ = lp.numCol_;
  all_cols.from_ = 0;
  all_cols.to_ = lp.numCol_ - 1;
  HighsIndexCollection all_rows;
  all_rows.is_interval_ = true;
  all_rows.dimension_ = lp.numRow_;
  all_rows.from_ = 0;
  all_rows.to_ = lp.numRow_ - 1;

  scale_error = applyScalingToLpColCost(options, lp, scale.col_, all_cols) !=
                    HighsStatus::OK ||
                scale_error;
  scale_error = applyScalingToLpColBounds(options, lp, scale.col_, all_cols) !=
                    HighsStatus::OK ||
                scale_error;
  scale_error = applyScalingToLpRowBounds(options, lp, scale.row_, all_rows) !=
                    HighsStatus::OK ||
                scale_error;
  scale_error = applyScalingToLpMatrix(options, lp, &scale.col_[0],
                                       &scale.row_[0], 0, lp.numCol_ - 1, 0,
                                       lp.numRow_ - 1) != HighsStatus::OK ||
                scale_error;
  if (scale_error) return HighsStatus::Error;
  return HighsStatus::OK;
}

HighsStatus applyScalingToLpColCost(
    const HighsOptions& options, HighsLp& lp, const vector<double>& colScale,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::OK;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");

  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::OK;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const int* col_set = index_collection.set_;
  const int* col_mask = index_collection.mask_;

  int local_col;
  int ml_col;
  const int ml_col_os = 0;
  for (int k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      local_col = k;
    } else {
      local_col = col_set[k];
    }
    ml_col = ml_col_os + local_col;
    if (mask && !col_mask[local_col]) continue;
    lp.colCost_[ml_col] *= colScale[ml_col];
  }

  return HighsStatus::OK;
}

HighsStatus applyScalingToLpColBounds(
    const HighsOptions& options, HighsLp& lp, const vector<double>& colScale,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::OK;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");

  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::OK;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const int* col_set = index_collection.set_;
  const int* col_mask = index_collection.mask_;

  int local_col;
  int ml_col;
  const int ml_col_os = 0;
  for (int k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      local_col = k;
    } else {
      local_col = col_set[k];
    }
    ml_col = ml_col_os + local_col;
    if (mask && !col_mask[local_col]) continue;
    if (!highs_isInfinity(-lp.colLower_[ml_col]))
      lp.colLower_[ml_col] /= colScale[ml_col];
    if (!highs_isInfinity(lp.colUpper_[ml_col]))
      lp.colUpper_[ml_col] /= colScale[ml_col];
  }

  return HighsStatus::OK;
}

HighsStatus applyScalingToLpRowBounds(
    const HighsOptions& options, HighsLp& lp, const vector<double>& rowScale,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::OK;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");

  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::OK;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const int* row_set = index_collection.set_;
  const int* row_mask = index_collection.mask_;

  int local_row;
  int ml_row;
  const int ml_row_os = 0;
  for (int k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      local_row = k;
    } else {
      local_row = row_set[k];
    }
    ml_row = ml_row_os + local_row;
    if (mask && !row_mask[local_row]) continue;
    if (!highs_isInfinity(-lp.rowLower_[ml_row]))
      lp.rowLower_[ml_row] *= rowScale[ml_row];
    if (!highs_isInfinity(lp.rowUpper_[ml_row]))
      lp.rowUpper_[ml_row] *= rowScale[ml_row];
  }

  return HighsStatus::OK;
}

HighsStatus applyScalingToLpMatrix(const HighsOptions& options, HighsLp& lp,
                                   const double* colScale,
                                   const double* rowScale, const int from_col,
                                   const int to_col, const int from_row,
                                   const int to_row) {
  if (from_col < 0) return HighsStatus::Error;
  if (to_col >= lp.numCol_) return HighsStatus::Error;
  if (from_row < 0) return HighsStatus::Error;
  if (to_row >= lp.numRow_) return HighsStatus::Error;
  if (colScale != NULL) {
    if (rowScale != NULL) {
      for (int iCol = from_col; iCol <= to_col; iCol++) {
        for (int iEl = lp.Astart_[iCol]; iEl < lp.Astart_[iCol + 1]; iEl++) {
          int iRow = lp.Aindex_[iEl];
          if (iRow < from_row || iRow > to_row) continue;
          lp.Avalue_[iEl] *= (colScale[iCol] * rowScale[iRow]);
        }
      }
    } else {
      // No row scaling
      for (int iCol = from_col; iCol <= to_col; iCol++) {
        for (int iEl = lp.Astart_[iCol]; iEl < lp.Astart_[iCol + 1]; iEl++) {
          int iRow = lp.Aindex_[iEl];
          if (iRow < from_row || iRow > to_row) continue;
          lp.Avalue_[iEl] *= colScale[iCol];
        }
      }
    }
  } else {
    // No column scaling
    if (rowScale != NULL) {
      for (int iCol = from_col; iCol <= to_col; iCol++) {
        for (int iEl = lp.Astart_[iCol]; iEl < lp.Astart_[iCol + 1]; iEl++) {
          int iRow = lp.Aindex_[iEl];
          if (iRow < from_row || iRow > to_row) continue;
          lp.Avalue_[iEl] *= rowScale[iRow];
        }
      }
    }
  }
  return HighsStatus::OK;
}

void applyRowScalingToMatrix(const vector<double>& rowScale, const int numCol,
                             const vector<int>& Astart,
                             const vector<int>& Aindex,
                             vector<double>& Avalue) {
  for (int iCol = 0; iCol < numCol; iCol++) {
    for (int el = Astart[iCol]; el < Astart[iCol + 1]; el++) {
      Avalue[el] *= rowScale[Aindex[el]];
    }
  }
}

void colScaleMatrix(const int max_scale_factor_exponent, double* colScale,
                    const int numCol, const vector<int>& Astart,
                    const vector<int>& Aindex, vector<double>& Avalue) {
  const double log2 = log(2.0);
  const double max_allow_scale = pow(2.0, max_scale_factor_exponent);
  const double min_allow_scale = 1 / max_allow_scale;

  const double min_allow_col_scale = min_allow_scale;
  const double max_allow_col_scale = max_allow_scale;

  for (int iCol = 0; iCol < numCol; iCol++) {
    double col_max_value = 0;
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
      col_max_value = max(fabs(Avalue[k]), col_max_value);
    if (col_max_value) {
      double col_scale_value = 1 / col_max_value;
      // Convert the col scale factor to the nearest power of two, and
      // ensure that it is not excessively large or small
      col_scale_value = pow(2.0, floor(log(col_scale_value) / log2 + 0.5));
      col_scale_value =
          min(max(min_allow_col_scale, col_scale_value), max_allow_col_scale);
      colScale[iCol] = col_scale_value;
      // Scale the column
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
        Avalue[k] *= colScale[iCol];
    } else {
      // Empty column
      colScale[iCol] = 1;
    }
  }
}

HighsStatus appendColsToLpVectors(HighsLp& lp, const int num_new_col,
                                  const vector<double>& colCost,
                                  const vector<double>& colLower,
                                  const vector<double>& colUpper) {
  if (num_new_col < 0) return HighsStatus::Error;
  if (num_new_col == 0) return HighsStatus::OK;
  int new_num_col = lp.numCol_ + num_new_col;
  lp.colCost_.resize(new_num_col);
  lp.colLower_.resize(new_num_col);
  lp.colUpper_.resize(new_num_col);
  bool have_names = lp.col_names_.size();
  if (have_names) lp.col_names_.resize(new_num_col);
  for (int new_col = 0; new_col < num_new_col; new_col++) {
    int iCol = lp.numCol_ + new_col;
    lp.colCost_[iCol] = colCost[new_col];
    lp.colLower_[iCol] = colLower[new_col];
    lp.colUpper_[iCol] = colUpper[new_col];
    // Cannot guarantee to create unique names, so name is blank
    if (have_names) lp.col_names_[iCol] = "";
  }
  return HighsStatus::OK;
}

HighsStatus appendRowsToLpVectors(HighsLp& lp, const int num_new_row,
                                  const vector<double>& rowLower,
                                  const vector<double>& rowUpper) {
  if (num_new_row < 0) return HighsStatus::Error;
  if (num_new_row == 0) return HighsStatus::OK;
  int new_num_row = lp.numRow_ + num_new_row;
  lp.rowLower_.resize(new_num_row);
  lp.rowUpper_.resize(new_num_row);
  bool have_names = lp.row_names_.size();
  if (have_names) lp.row_names_.resize(new_num_row);

  for (int new_row = 0; new_row < num_new_row; new_row++) {
    int iRow = lp.numRow_ + new_row;
    lp.rowLower_[iRow] = rowLower[new_row];
    lp.rowUpper_[iRow] = rowUpper[new_row];
    // Cannot guarantee to create unique names, so name is blank
    if (have_names) lp.row_names_[iRow] = "";
  }
  return HighsStatus::OK;
}

HighsStatus appendColsToLpMatrix(HighsLp& lp, const int num_new_col,
                                 const int num_new_nz, const int* XAstart,
                                 const int* XAindex, const double* XAvalue) {
  if (num_new_col < 0) return HighsStatus::Error;
  if (num_new_col == 0) return HighsStatus::OK;
  // Check that nonzeros aren't being appended to a matrix with no rows
  if (num_new_nz > 0 && lp.numRow_ <= 0) return HighsStatus::Error;
  // Determine the new number of columns in the matrix and resize the
  // starts accordingly.
  int new_num_col = lp.numCol_ + num_new_col;
  lp.Astart_.resize(new_num_col + 1);
  // If adding columns to an empty LP then introduce the start for the
  // fictitious column 0
  if (lp.numCol_ == 0) lp.Astart_[0] = 0;

  // Determine the current number of nonzeros and the new number of nonzeros
  int current_num_nz = lp.Astart_[lp.numCol_];
  int new_num_nz = current_num_nz + num_new_nz;

  // Append the starts of the new columns
  if (num_new_nz) {
    // Nontrivial number of nonzeros being added, so use XAstart
    assert(XAstart != NULL);
    for (int col = 0; col < num_new_col; col++)
      lp.Astart_[lp.numCol_ + col] = current_num_nz + XAstart[col];
  } else {
    // No nonzeros being added, so XAstart may be null, but entries of
    // zero are implied.
    for (int col = 0; col < num_new_col; col++)
      lp.Astart_[lp.numCol_ + col] = current_num_nz;
  }
  lp.Astart_[lp.numCol_ + num_new_col] = new_num_nz;

  // If no nonzeros are being added then there's nothing else to do
  if (num_new_nz <= 0) return HighsStatus::OK;

  // Adding a non-trivial matrix: resize the column-wise matrix arrays
  // accordingly
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);
  // Copy in the new indices and values
  for (int el = 0; el < num_new_nz; el++) {
    lp.Aindex_[current_num_nz + el] = XAindex[el];
    lp.Avalue_[current_num_nz + el] = XAvalue[el];
  }
  return HighsStatus::OK;
}

HighsStatus appendRowsToLpMatrix(HighsLp& lp, const int num_new_row,
                                 const int num_new_nz, const int* XARstart,
                                 const int* XARindex, const double* XARvalue) {
  if (num_new_row < 0) return HighsStatus::Error;
  if (num_new_row == 0) return HighsStatus::OK;
  // Check that nonzeros aren't being appended to a matrix with no columns
  if (num_new_nz > 0 && lp.numCol_ <= 0) return HighsStatus::Error;
  if (num_new_nz == 0) return HighsStatus::OK;
  int current_num_nz = lp.Astart_[lp.numCol_];
  vector<int> Alength;
  Alength.assign(lp.numCol_, 0);
  for (int el = 0; el < num_new_nz; el++) Alength[XARindex[el]]++;
  // Determine the new number of nonzeros and resize the column-wise matrix
  // arrays
  int new_num_nz = current_num_nz + num_new_nz;
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);

  // Append the new rows
  // Shift the existing columns to make space for the new entries
  int new_el = new_num_nz;
  for (int col = lp.numCol_ - 1; col >= 0; col--) {
    int start_col_plus_1 = new_el;
    new_el -= Alength[col];
    for (int el = lp.Astart_[col + 1] - 1; el >= lp.Astart_[col]; el--) {
      new_el--;
      lp.Aindex_[new_el] = lp.Aindex_[el];
      lp.Avalue_[new_el] = lp.Avalue_[el];
    }
    lp.Astart_[col + 1] = start_col_plus_1;
  }
  assert(new_el == 0);

  // Insert the new entries
  for (int row = 0; row < num_new_row; row++) {
    int first_el = XARstart[row];
    int last_el = (row < num_new_row - 1 ? XARstart[row + 1] : num_new_nz);
    for (int el = first_el; el < last_el; el++) {
      int col = XARindex[el];
      new_el = lp.Astart_[col + 1] - Alength[col];
      Alength[col]--;
      lp.Aindex_[new_el] = lp.numRow_ + row;
      lp.Avalue_[new_el] = XARvalue[el];
    }
  }
  return HighsStatus::OK;
}

HighsStatus deleteLpCols(const HighsOptions& options, HighsLp& lp,
                         const HighsIndexCollection& index_collection) {
  int new_num_col;
  HighsStatus call_status;
  call_status =
      deleteColsFromLpVectors(options, lp, new_num_col, index_collection);
  if (call_status != HighsStatus::OK) return call_status;
  call_status = deleteColsFromLpMatrix(options, lp, index_collection);
  if (call_status != HighsStatus::OK) return call_status;
  lp.numCol_ = new_num_col;
  return HighsStatus::OK;
}

HighsStatus deleteColsFromLpVectors(
    const HighsOptions& options, HighsLp& lp, int& new_num_col,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::OK;
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0, lp.numCol_ - 1,
                         true))
      return HighsStatus::Error;
  }
  // Initialise new_num_col in case none is removed due to from_k > to_k
  new_num_col = lp.numCol_;
  if (from_k > to_k) return HighsStatus::OK;

  int delete_from_col;
  int delete_to_col;
  int keep_from_col;
  int keep_to_col = -1;
  int current_set_entry = 0;

  int col_dim = lp.numCol_;
  new_num_col = 0;
  bool have_names = lp.col_names_.size();
  for (int k = from_k; k <= to_k; k++) {
    updateIndexCollectionOutInIndex(index_collection, delete_from_col,
                                    delete_to_col, keep_from_col, keep_to_col,
                                    current_set_entry);
    // Account for the initial columns being kept
    if (k == from_k) new_num_col = delete_from_col;
    if (delete_to_col >= col_dim - 1) break;
    assert(delete_to_col < col_dim);
    for (int col = keep_from_col; col <= keep_to_col; col++) {
      lp.colCost_[new_num_col] = lp.colCost_[col];
      lp.colLower_[new_num_col] = lp.colLower_[col];
      lp.colUpper_[new_num_col] = lp.colUpper_[col];
      if (have_names) lp.col_names_[new_num_col] = lp.col_names_[col];
      new_num_col++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
  lp.colCost_.resize(new_num_col);
  lp.colLower_.resize(new_num_col);
  lp.colUpper_.resize(new_num_col);
  if (have_names) lp.col_names_.resize(new_num_col);
  return HighsStatus::OK;
}

HighsStatus deleteColsFromLpMatrix(
    const HighsOptions& options, HighsLp& lp,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::OK;
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0, lp.numCol_ - 1,
                         true))
      return HighsStatus::Error;
  }
  if (from_k > to_k) return HighsStatus::OK;

  int delete_from_col;
  int delete_to_col;
  int keep_from_col;
  int keep_to_col = -1;
  int current_set_entry = 0;

  int col_dim = lp.numCol_;
  int new_num_col = 0;
  int new_num_nz = 0;
  for (int k = from_k; k <= to_k; k++) {
    updateIndexCollectionOutInIndex(index_collection, delete_from_col,
                                    delete_to_col, keep_from_col, keep_to_col,
                                    current_set_entry);
    if (k == from_k) {
      // Account for the initial columns being kept
      new_num_col = delete_from_col;
      new_num_nz = lp.Astart_[delete_from_col];
    }
    // Ensure that the starts of the deleted columns are zeroed to
    // avoid redundant start information for columns whose indices
    // are't used after the deletion takes place. In particular, if
    // all columns are deleted then something must be done to ensure
    // that the matrix isn't magially recreated by increasing the
    // number of columns from zero when there are no rows in the LP.
    for (int col = delete_from_col; col <= delete_to_col; col++)
      lp.Astart_[col] = 0;
    for (int col = keep_from_col; col <= keep_to_col; col++) {
      lp.Astart_[new_num_col] =
          new_num_nz + lp.Astart_[col] - lp.Astart_[keep_from_col];
      new_num_col++;
    }
    for (int el = lp.Astart_[keep_from_col]; el < lp.Astart_[keep_to_col + 1];
         el++) {
      lp.Aindex_[new_num_nz] = lp.Aindex_[el];
      lp.Avalue_[new_num_nz] = lp.Avalue_[el];
      new_num_nz++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
  // Ensure that the start of the spurious last column is zeroed so
  // that it doesn't give a positive number of matrix entries if the
  // number of columns in the LP is increased when there are no rows
  // in the LP.
  lp.Astart_[lp.numCol_] = 0;
  lp.Astart_[new_num_col] = new_num_nz;
  lp.Astart_.resize(new_num_col + 1);
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);
  return HighsStatus::OK;
}

HighsStatus deleteLpRows(const HighsOptions& options, HighsLp& lp,
                         const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  int new_num_row;
  call_status =
      deleteRowsFromLpVectors(options, lp, new_num_row, index_collection);
  return_status = interpretCallStatus(call_status, return_status,
                                      "deleteRowsFromLpVectors");
  if (return_status == HighsStatus::Error) return return_status;
  call_status = deleteRowsFromLpMatrix(options, lp, index_collection);
  return_status =
      interpretCallStatus(call_status, return_status, "deleteRowsFromLpMatrix");
  if (return_status == HighsStatus::Error) return return_status;
  lp.numRow_ = new_num_row;
  return HighsStatus::OK;
}

HighsStatus deleteRowsFromLpVectors(
    const HighsOptions& options, HighsLp& lp, int& new_num_row,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::OK;
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0, lp.numRow_ - 1,
                         true))
      return HighsStatus::Error;
  }
  // Initialise new_num_row in case none is removed due to from_k > to_k
  new_num_row = lp.numRow_;
  if (from_k > to_k) return HighsStatus::OK;

  int delete_from_row;
  int delete_to_row;
  int keep_from_row;
  int keep_to_row = -1;
  int current_set_entry = 0;

  int row_dim = lp.numRow_;
  new_num_row = 0;
  bool have_names = lp.row_names_.size();
  for (int k = from_k; k <= to_k; k++) {
    updateIndexCollectionOutInIndex(index_collection, delete_from_row,
                                    delete_to_row, keep_from_row, keep_to_row,
                                    current_set_entry);
    if (k == from_k) {
      // Account for the initial rows being kept
      new_num_row = delete_from_row;
    }
    if (delete_to_row >= row_dim - 1) break;
    assert(delete_to_row < row_dim);
    for (int row = keep_from_row; row <= keep_to_row; row++) {
      lp.rowLower_[new_num_row] = lp.rowLower_[row];
      lp.rowUpper_[new_num_row] = lp.rowUpper_[row];
      if (have_names) lp.row_names_[new_num_row] = lp.row_names_[row];
      new_num_row++;
    }
    if (keep_to_row == row_dim) break;
  }
  lp.rowLower_.resize(new_num_row);
  lp.rowUpper_.resize(new_num_row);
  if (have_names) lp.row_names_.resize(new_num_row);
  return HighsStatus::OK;
}

HighsStatus deleteRowsFromLpMatrix(
    const HighsOptions& options, HighsLp& lp,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::OK;
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0, lp.numRow_ - 1,
                         true))
      return HighsStatus::Error;
  }
  if (from_k > to_k) return HighsStatus::OK;

  int delete_from_row;
  int delete_to_row;
  int keep_from_row;
  int row_dim = lp.numRow_;
  int keep_to_row = -1;
  int current_set_entry = 0;

  // Set up a row mask to indicate the new row index of kept rows and
  // -1 for deleted rows so that the kept entries in the column-wise
  // matrix can be identified and have their correct row index.
  int* new_index = (int*)malloc(sizeof(int) * lp.numRow_);
  int new_num_row = 0;
  bool mask = index_collection.is_mask_;
  const int* row_mask = index_collection.mask_;
  if (!mask) {
    keep_to_row = -1;
    current_set_entry = 0;
    for (int k = from_k; k <= to_k; k++) {
      updateIndexCollectionOutInIndex(index_collection, delete_from_row,
                                      delete_to_row, keep_from_row, keep_to_row,
                                      current_set_entry);
      if (k == from_k) {
        // Account for any initial rows being kept
        for (int row = 0; row < delete_from_row; row++) {
          new_index[row] = new_num_row;
          new_num_row++;
        }
      }
      for (int row = delete_from_row; row <= delete_to_row; row++) {
        new_index[row] = -1;
      }
      for (int row = keep_from_row; row <= keep_to_row; row++) {
        new_index[row] = new_num_row;
        new_num_row++;
      }
      if (keep_to_row >= row_dim - 1) break;
    }
  } else {
    for (int row = 0; row < lp.numRow_; row++) {
      if (row_mask[row]) {
        new_index[row] = -1;
      } else {
        new_index[row] = new_num_row;
        new_num_row++;
      }
    }
  }
  int new_num_nz = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    int from_el = lp.Astart_[col];
    lp.Astart_[col] = new_num_nz;
    for (int el = from_el; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      int new_row = new_index[row];
      if (new_row >= 0) {
        lp.Aindex_[new_num_nz] = new_row;
        lp.Avalue_[new_num_nz] = lp.Avalue_[el];
        new_num_nz++;
      }
    }
  }
  lp.Astart_[lp.numCol_] = new_num_nz;
  free(new_index);
  lp.Astart_.resize(lp.numCol_ + 1);
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);
  return HighsStatus::OK;
}

HighsStatus changeLpMatrixCoefficient(HighsLp& lp, const int row, const int col,
                                      const double new_value) {
  if (row < 0 || row > lp.numRow_) return HighsStatus::Error;
  if (col < 0 || col > lp.numCol_) return HighsStatus::Error;
  int changeElement = -1;
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
    // printf("Column %d: Element %d is row %d. Is it %d?\n", col, el,
    // lp.Aindex_[el], row);
    if (lp.Aindex_[el] == row) {
      changeElement = el;
      break;
    }
  }
  if (changeElement < 0) {
    //    printf("changeLpMatrixCoefficient: Cannot find row %d in column %d\n",
    //    row, col);
    changeElement = lp.Astart_[col + 1];
    int new_num_nz = lp.Astart_[lp.numCol_] + 1;
    //    printf("changeLpMatrixCoefficient: Increasing Nnonz from %d to %d\n",
    //    lp.Astart_[lp.numCol_], new_num_nz);
    lp.Aindex_.resize(new_num_nz);
    lp.Avalue_.resize(new_num_nz);
    for (int i = col + 1; i <= lp.numCol_; i++) lp.Astart_[i]++;
    for (int el = new_num_nz - 1; el > changeElement; el--) {
      lp.Aindex_[el] = lp.Aindex_[el - 1];
      lp.Avalue_[el] = lp.Avalue_[el - 1];
    }
  }
  lp.Aindex_[changeElement] = row;
  lp.Avalue_[changeElement] = new_value;

  return HighsStatus::OK;
}

HighsStatus changeLpCosts(const HighsOptions& options, HighsLp& lp,
                          const HighsIndexCollection& index_collection,
                          const vector<double>& new_col_cost) {
  HighsStatus return_status = HighsStatus::OK;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::OK;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const int* col_set = index_collection.set_;
  const int* col_mask = index_collection.mask_;

  // Change the costs to the user-supplied costs, according to the technique
  int usr_col;
  for (int k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      usr_col = k;
    } else {
      usr_col = col_set[k];
    }
    int col = usr_col;
    if (mask && !col_mask[col]) continue;
    lp.colCost_[col] = new_col_cost[k];
  }
  return HighsStatus::OK;
}

HighsStatus changeLpColBounds(const HighsOptions& options, HighsLp& lp,
                              const HighsIndexCollection& index_collection,
                              const vector<double>& new_col_lower,
                              const vector<double>& new_col_upper) {
  return changeBounds(options, lp.colLower_, lp.colUpper_, index_collection,
                      new_col_lower, new_col_upper);
}

HighsStatus changeLpRowBounds(const HighsOptions& options, HighsLp& lp,
                              const HighsIndexCollection& index_collection,
                              const vector<double>& new_row_lower,
                              const vector<double>& new_row_upper) {
  return changeBounds(options, lp.rowLower_, lp.rowUpper_, index_collection,
                      new_row_lower, new_row_upper);
}

HighsStatus changeBounds(const HighsOptions& options, vector<double>& lower,
                         vector<double>& upper,
                         const HighsIndexCollection& index_collection,
                         const vector<double>& new_lower,
                         const vector<double>& new_upper) {
  HighsStatus return_status = HighsStatus::OK;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::OK;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const int* ix_set = index_collection.set_;
  const int* ix_mask = index_collection.mask_;

  // Change the bounds to the user-supplied bounds, according to the technique
  int usr_ix;
  for (int k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      usr_ix = k;
    } else {
      usr_ix = ix_set[k];
    }
    int ix = usr_ix;
    if (mask && !ix_mask[ix]) continue;
    lower[ix] = new_lower[k];
    upper[ix] = new_upper[k];
  }
  return HighsStatus::OK;
}

int getNumInt(const HighsLp& lp) {
  int num_int = 0;
  if (lp.integrality_.size()) {
    for (int iCol = 0; iCol < lp.numCol_; iCol++)
      if (lp.integrality_[iCol]) num_int++;
  }
  return num_int;
}

HighsStatus getLpCosts(const HighsLp& lp, const int from_col, const int to_col,
                       double* XcolCost) {
  if (from_col < 0 || to_col >= lp.numCol_) return HighsStatus::Error;
  if (from_col > to_col) return HighsStatus::OK;
  for (int col = from_col; col < to_col + 1; col++)
    XcolCost[col - from_col] = lp.colCost_[col];
  return HighsStatus::OK;
}

HighsStatus getLpColBounds(const HighsLp& lp, const int from_col,
                           const int to_col, double* XcolLower,
                           double* XcolUpper) {
  if (from_col < 0 || to_col >= lp.numCol_) return HighsStatus::Error;
  if (from_col > to_col) return HighsStatus::OK;
  for (int col = from_col; col < to_col + 1; col++) {
    if (XcolLower != NULL) XcolLower[col - from_col] = lp.colLower_[col];
    if (XcolUpper != NULL) XcolUpper[col - from_col] = lp.colUpper_[col];
  }
  return HighsStatus::OK;
}

HighsStatus getLpRowBounds(const HighsLp& lp, const int from_row,
                           const int to_row, double* XrowLower,
                           double* XrowUpper) {
  if (from_row < 0 || to_row >= lp.numRow_) return HighsStatus::Error;
  if (from_row > to_row) return HighsStatus::OK;
  for (int row = from_row; row < to_row + 1; row++) {
    if (XrowLower != NULL) XrowLower[row - from_row] = lp.rowLower_[row];
    if (XrowUpper != NULL) XrowUpper[row - from_row] = lp.rowUpper_[row];
  }
  return HighsStatus::OK;
}

// Get a single coefficient from the matrix
HighsStatus getLpMatrixCoefficient(const HighsLp& lp, const int Xrow,
                                   const int Xcol, double* val) {
#ifdef HiGHSDEV
  printf("Called getLpMatrixCoefficient(row=%d, col=%d)\n", Xrow, Xcol);
#endif
  if (Xrow < 0 || Xrow >= lp.numRow_) return HighsStatus::Error;
  if (Xcol < 0 || Xcol >= lp.numCol_) return HighsStatus::Error;

  int get_el = -1;
  for (int el = lp.Astart_[Xcol]; el < lp.Astart_[Xcol + 1]; el++) {
    if (lp.Aindex_[el] == Xrow) {
      get_el = el;
      break;
    }
  }
  if (get_el < 0) {
    *val = 0;
  } else {
    *val = lp.Avalue_[get_el];
  }
  return HighsStatus::OK;
}

// Methods for reporting an LP, including its row and column data and matrix
//
// Report the whole LP
void reportLp(const HighsOptions& options, const HighsLp& lp,
              const int report_level) {
  reportLpBrief(options, lp);
  if (report_level >= 1) {
    reportLpColVectors(options, lp);
    reportLpRowVectors(options, lp);
    if (report_level >= 2) reportLpColMatrix(options, lp);
  }
}

// Report the LP briefly
void reportLpBrief(const HighsOptions& options, const HighsLp& lp) {
  reportLpDimensions(options, lp);
  reportLpObjSense(options, lp);
}

// Report the LP dimensions
void reportLpDimensions(const HighsOptions& options, const HighsLp& lp) {
  int lp_num_nz;
  if (lp.numCol_ == 0)
    lp_num_nz = 0;
  else
    lp_num_nz = lp.Astart_[lp.numCol_];
  HighsPrintMessage(options.output, options.message_level, ML_MINIMAL,
                    "LP has %d columns, %d rows", lp.numCol_, lp.numRow_);
  int num_int = getNumInt(lp);
  if (num_int) {
    HighsPrintMessage(options.output, options.message_level, ML_MINIMAL,
                      ", %d nonzeros and %d integer columns\n", lp_num_nz,
                      num_int);
  } else {
    HighsPrintMessage(options.output, options.message_level, ML_MINIMAL,
                      " and %d nonzeros\n", lp_num_nz, num_int);
  }
}

// Report the LP objective sense
void reportLpObjSense(const HighsOptions& options, const HighsLp& lp) {
  if (lp.sense_ == ObjSense::MINIMIZE)
    HighsPrintMessage(options.output, options.message_level, ML_MINIMAL,
                      "Objective sense is minimize\n");
  else if (lp.sense_ == ObjSense::MAXIMIZE)
    HighsPrintMessage(options.output, options.message_level, ML_MINIMAL,
                      "Objective sense is maximize\n");
  else
    HighsPrintMessage(options.output, options.message_level, ML_MINIMAL,
                      "Objective sense is ill-defined as %d\n", lp.sense_);
}

std::string getBoundType(const double lower, const double upper) {
  std::string type;
  if (highs_isInfinity(-lower)) {
    if (highs_isInfinity(upper)) {
      type = "FR";
    } else {
      type = "UB";
    }
  } else {
    if (highs_isInfinity(upper)) {
      type = "LB";
    } else {
      if (lower < upper) {
        type = "BX";
      } else {
        type = "FX";
      }
    }
  }
  return type;
}

// Report the vectors of LP column data
void reportLpColVectors(const HighsOptions& options, const HighsLp& lp) {
  if (lp.numCol_ <= 0) return;
  std::string type;
  int count;
  bool have_integer_columns = getNumInt(lp);
  bool have_col_names = lp.col_names_.size();

  HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                    "  Column        Lower        Upper         Cost       "
                    "Type        Count");
  if (have_integer_columns)
    HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                      "  Discrete");
  if (have_col_names)
    HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                      "  Name");
  HighsPrintMessage(options.output, options.message_level, ML_VERBOSE, "\n");

  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    type = getBoundType(lp.colLower_[iCol], lp.colUpper_[iCol]);
    count = lp.Astart_[iCol + 1] - lp.Astart_[iCol];
    HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                      "%8d %12g %12g %12g         %2s %12d", iCol,
                      lp.colLower_[iCol], lp.colUpper_[iCol], lp.colCost_[iCol],
                      type.c_str(), count);
    if (have_integer_columns) {
      std::string integer_column = "";
      if (lp.integrality_[iCol]) {
        if (lp.colLower_[iCol] == 0 && lp.colUpper_[iCol] == 1) {
          integer_column = "Binary";
        } else {
          integer_column = "Integer";
        }
      }
      HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                        "  %-8s", integer_column.c_str());
    }
    if (have_col_names)
      HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                        "  %-s", lp.col_names_[iCol].c_str());
    HighsPrintMessage(options.output, options.message_level, ML_VERBOSE, "\n");
  }
}

// Report the vectors of LP row data
void reportLpRowVectors(const HighsOptions& options, const HighsLp& lp) {
  if (lp.numRow_ <= 0) return;
  std::string type;
  vector<int> count;
  bool have_row_names = lp.row_names_.size();

  count.resize(lp.numRow_, 0);
  if (lp.numCol_ > 0) {
    for (int el = 0; el < lp.Astart_[lp.numCol_]; el++) count[lp.Aindex_[el]]++;
  }

  HighsPrintMessage(
      options.output, options.message_level, ML_VERBOSE,
      "     Row        Lower        Upper       Type        Count");
  if (have_row_names)
    HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                      "  Name");
  HighsPrintMessage(options.output, options.message_level, ML_VERBOSE, "\n");

  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    type = getBoundType(lp.rowLower_[iRow], lp.rowUpper_[iRow]);
    std::string name = "";
    HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                      "%8d %12g %12g         %2s %12d", iRow,
                      lp.rowLower_[iRow], lp.rowUpper_[iRow], type.c_str(),
                      count[iRow]);
    if (have_row_names)
      HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                        "  %-s", lp.row_names_[iRow].c_str());
    HighsPrintMessage(options.output, options.message_level, ML_VERBOSE, "\n");
  }
}

// Report the LP column-wise matrix
void reportLpColMatrix(const HighsOptions& options, const HighsLp& lp) {
  if (lp.numCol_ <= 0) return;
  if (lp.numRow_) {
    // With postitive number of rows, can assume that there are index and value
    // vectors to pass
    reportMatrix(options, "Column", lp.numCol_, lp.Astart_[lp.numCol_],
                 &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);
  } else {
    // With no rows, can's assume that there are index and value vectors to pass
    reportMatrix(options, "Column", lp.numCol_, lp.Astart_[lp.numCol_],
                 &lp.Astart_[0], NULL, NULL);
  }
}

void reportMatrix(const HighsOptions& options, const std::string message,
                  const int num_col, const int num_nz, const int* start,
                  const int* index, const double* value) {
  if (num_col <= 0) return;
  HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                    "%6s Index              Value\n", message.c_str());
  for (int col = 0; col < num_col; col++) {
    HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                      "    %8d Start   %10d\n", col, start[col]);
    int to_el = (col < num_col - 1 ? start[col + 1] : num_nz);
    for (int el = start[col]; el < to_el; el++)
      HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                        "          %8d %12g\n", index[el], value[el]);
  }
  HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                    "             Start   %10d\n", num_nz);
}

#ifdef HiGHSDEV
void analyseLp(const HighsLp& lp, const std::string message) {
  vector<double> min_colBound;
  vector<double> min_rowBound;
  vector<double> colRange;
  vector<double> rowRange;
  min_colBound.resize(lp.numCol_);
  min_rowBound.resize(lp.numRow_);
  colRange.resize(lp.numCol_);
  rowRange.resize(lp.numRow_);
  for (int col = 0; col < lp.numCol_; col++)
    min_colBound[col] = min(fabs(lp.colLower_[col]), fabs(lp.colUpper_[col]));
  for (int row = 0; row < lp.numRow_; row++)
    min_rowBound[row] = min(fabs(lp.rowLower_[row]), fabs(lp.rowUpper_[row]));
  for (int col = 0; col < lp.numCol_; col++)
    colRange[col] = lp.colUpper_[col] - lp.colLower_[col];
  for (int row = 0; row < lp.numRow_; row++)
    rowRange[row] = lp.rowUpper_[row] - lp.rowLower_[row];

  printf("\n%s model data: Analysis\n", message.c_str());
  analyseVectorValues("Column costs", lp.numCol_, lp.colCost_);
  analyseVectorValues("Column lower bounds", lp.numCol_, lp.colLower_);
  analyseVectorValues("Column upper bounds", lp.numCol_, lp.colUpper_);
  analyseVectorValues("Column min abs bound", lp.numCol_, min_colBound);
  analyseVectorValues("Column range", lp.numCol_, colRange);
  analyseVectorValues("Row lower bounds", lp.numRow_, lp.rowLower_);
  analyseVectorValues("Row upper bounds", lp.numRow_, lp.rowUpper_);
  analyseVectorValues("Row min abs bound", lp.numRow_, min_rowBound);
  analyseVectorValues("Row range", lp.numRow_, rowRange);
  analyseVectorValues("Matrix sparsity", lp.Astart_[lp.numCol_], lp.Avalue_,
                      true, lp.model_name_);
  analyseMatrixSparsity("Constraint matrix", lp.numCol_, lp.numRow_, lp.Astart_,
                        lp.Aindex_);
  analyseModelBounds("Column", lp.numCol_, lp.colLower_, lp.colUpper_);
  analyseModelBounds("Row", lp.numRow_, lp.rowLower_, lp.rowUpper_);
}
#endif

// void writeSolutionToFile(FILE* file, const HighsLp& lp, const HighsBasis&
// basis,
//                          const HighsSolution& solution, const bool pretty) {
//   if (pretty) {
//     reportModelBoundSol(file, true, lp.numCol_, lp.colLower_, lp.colUpper_,
//                         lp.col_names_, solution.col_value, solution.col_dual,
//                         basis.col_status);
//     reportModelBoundSol(file, false, lp.numRow_, lp.rowLower_, lp.rowUpper_,
//                         lp.row_names_, solution.row_value, solution.row_dual,
//                         basis.row_status);
//   } else {
//     fprintf(file,
//             "%d %d : Number of columns and rows for primal and dual solution
//             " "and basis\n", lp.numCol_, lp.numRow_);
//     const bool with_basis = basis.valid_;
//     if (with_basis) {
//       fprintf(file, "T\n");
//     } else {
//       fprintf(file, "F\n");
//     }
//     for (int iCol = 0; iCol < lp.numCol_; iCol++) {
//       fprintf(file, "%g %g", solution.col_value[iCol],
//       solution.col_dual[iCol]); if (with_basis) fprintf(file, " %d",
//       (int)basis.col_status[iCol]); fprintf(file, " \n");
//     }
//     for (int iRow = 0; iRow < lp.numRow_; iRow++) {
//       fprintf(file, "%g %g", solution.row_value[iRow],
//       solution.row_dual[iRow]); if (with_basis) fprintf(file, " %d",
//       (int)basis.row_status[iRow]); fprintf(file, " \n");
//     }
//   }
// }

HighsStatus convertBasis(const HighsLp& lp, const SimplexBasis& basis,
                         HighsBasis& new_basis) {
  new_basis.col_status.clear();
  new_basis.row_status.clear();

  new_basis.col_status.resize(lp.numCol_);
  new_basis.row_status.resize(lp.numRow_);

  for (int col = 0; col < lp.numCol_; col++) {
    if (!basis.nonbasicFlag_[col]) {
      new_basis.col_status[col] = HighsBasisStatus::BASIC;
    } else if (basis.nonbasicMove_[col] == NONBASIC_MOVE_UP) {
      new_basis.col_status[col] = HighsBasisStatus::LOWER;
    } else if (basis.nonbasicMove_[col] == NONBASIC_MOVE_DN) {
      new_basis.col_status[col] = HighsBasisStatus::UPPER;
    } else if (basis.nonbasicMove_[col] == NONBASIC_MOVE_ZE) {
      if (lp.colLower_[col] == lp.colUpper_[col]) {
        new_basis.col_status[col] = HighsBasisStatus::LOWER;
      } else {
        new_basis.col_status[col] = HighsBasisStatus::ZERO;
      }
    } else {
      return HighsStatus::Error;
    }
  }

  for (int row = 0; row < lp.numRow_; row++) {
    int var = lp.numCol_ + row;
    if (!basis.nonbasicFlag_[var]) {
      new_basis.row_status[row] = HighsBasisStatus::BASIC;
    } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN) {
      new_basis.row_status[row] = HighsBasisStatus::LOWER;
    } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
      new_basis.row_status[row] = HighsBasisStatus::UPPER;
    } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
      if (lp.rowLower_[row] == lp.rowUpper_[row]) {
        new_basis.row_status[row] = HighsBasisStatus::LOWER;
      } else {
        new_basis.row_status[row] = HighsBasisStatus::ZERO;
      }
    } else {
      return HighsStatus::Error;
    }
  }

  return HighsStatus::OK;
}

HighsBasis getSimplexBasis(const HighsLp& lp, const SimplexBasis& basis) {
  HighsBasis new_basis;
  HighsStatus result = convertBasis(lp, basis, new_basis);
  if (result != HighsStatus::OK) return HighsBasis();
  // Call Julian's code to translate basis once it's out of
  // SimplexInterface. Until it is out of SimplexInteface use code
  // I just added above which does the same but only returns an
  // error and not which basis index has an illegal value.
  return new_basis;
}

HighsStatus calculateColDuals(const HighsLp& lp, HighsSolution& solution) {
  assert(solution.row_dual.size() > 0);
  if (!isSolutionRightSize(lp, solution)) return HighsStatus::Error;

  solution.col_dual.assign(lp.numCol_, 0);

  for (int col = 0; col < lp.numCol_; col++) {
    for (int i = lp.Astart_[col]; i < lp.Astart_[col + 1]; i++) {
      const int row = lp.Aindex_[i];
      assert(row >= 0);
      assert(row < lp.numRow_);

      solution.col_dual[col] -= solution.row_dual[row] * lp.Avalue_[i];
    }
    solution.col_dual[col] += lp.colCost_[col];
  }

  return HighsStatus::OK;
}

HighsStatus calculateRowValues(const HighsLp& lp, HighsSolution& solution) {
  assert(solution.col_value.size() > 0);
  if (!isSolutionRightSize(lp, solution)) return HighsStatus::Error;

  solution.row_value.clear();
  solution.row_value.assign(lp.numRow_, 0);

  for (int col = 0; col < lp.numCol_; col++) {
    for (int i = lp.Astart_[col]; i < lp.Astart_[col + 1]; i++) {
      const int row = lp.Aindex_[i];
      assert(row >= 0);
      assert(row < lp.numRow_);

      solution.row_value[row] += solution.col_value[col] * lp.Avalue_[i];
    }
  }

  return HighsStatus::OK;
}

double calculateObjective(const HighsLp& lp, HighsSolution& solution) {
  assert(isSolutionRightSize(lp, solution));
  double sum = 0;
  for (int col = 0; col < lp.numCol_; col++)
    sum += lp.colCost_[col] * solution.col_value[col];

  return sum;
}

bool isColDataNull(const HighsOptions& options, const double* usr_col_cost,
                   const double* usr_col_lower, const double* usr_col_upper) {
  bool null_data = false;
  null_data =
      doubleUserDataNotNull(options.logfile, usr_col_cost, "column costs") ||
      null_data;
  null_data = doubleUserDataNotNull(options.logfile, usr_col_lower,
                                    "column lower bounds") ||
              null_data;
  null_data = doubleUserDataNotNull(options.logfile, usr_col_upper,
                                    "column upper bounds") ||
              null_data;
  return null_data;
}

bool isRowDataNull(const HighsOptions& options, const double* usr_row_lower,
                   const double* usr_row_upper) {
  bool null_data = false;
  null_data = doubleUserDataNotNull(options.logfile, usr_row_lower,
                                    "row lower bounds") ||
              null_data;
  null_data = doubleUserDataNotNull(options.logfile, usr_row_upper,
                                    "row upper bounds") ||
              null_data;
  return null_data;
}

bool isMatrixDataNull(const HighsOptions& options, const int* usr_matrix_start,
                      const int* usr_matrix_index,
                      const double* usr_matrix_value) {
  bool null_data = false;
  null_data =
      intUserDataNotNull(options.logfile, usr_matrix_start, "matrix starts") ||
      null_data;
  null_data =
      intUserDataNotNull(options.logfile, usr_matrix_index, "matrix indices") ||
      null_data;
  null_data = doubleUserDataNotNull(options.logfile, usr_matrix_value,
                                    "matrix values") ||
              null_data;
  return null_data;
}

HighsStatus transformIntoEqualityProblem(const HighsLp& lp,
                                         HighsLp& equality_lp) {
  // Copy lp.
  equality_lp = lp;

  // Add slacks for each row with more than one bound.
  std::vector<double> rhs(lp.numRow_, 0);

  for (int row = 0; row < lp.numRow_; row++) {
    assert(equality_lp.Astart_[equality_lp.numCol_] ==
           (int)equality_lp.Avalue_.size());
    assert((int)equality_lp.Aindex_.size() == (int)equality_lp.Avalue_.size());
    const int nnz = equality_lp.Astart_[equality_lp.numCol_];

    if (lp.rowLower_[row] <= -HIGHS_CONST_INF &&
        lp.rowUpper_[row] >= HIGHS_CONST_INF) {
      // free row
      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(1.0);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(-HIGHS_CONST_INF);
      equality_lp.colUpper_.push_back(HIGHS_CONST_INF);
      equality_lp.colCost_.push_back(0);
    } else if (lp.rowLower_[row] > -HIGHS_CONST_INF &&
               lp.rowUpper_[row] >= HIGHS_CONST_INF) {
      // only lower bound
      rhs[row] = lp.rowLower_[row];

      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(-1.0);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(0);
      equality_lp.colUpper_.push_back(HIGHS_CONST_INF);
      equality_lp.colCost_.push_back(0);
    } else if (lp.rowLower_[row] <= -HIGHS_CONST_INF &&
               lp.rowUpper_[row] < HIGHS_CONST_INF) {
      // only upper bound
      rhs[row] = lp.rowUpper_[row];

      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(1.0);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(0);
      equality_lp.colUpper_.push_back(HIGHS_CONST_INF);
      equality_lp.colCost_.push_back(0);
    } else if (lp.rowLower_[row] > -HIGHS_CONST_INF &&
               lp.rowUpper_[row] < HIGHS_CONST_INF &&
               lp.rowLower_[row] != lp.rowUpper_[row]) {
      // both lower and upper bound that are different
      double rhs_value, coefficient;
      double difference = lp.rowUpper_[row] - lp.rowLower_[row];
      if (fabs(lp.rowLower_[row]) < fabs(lp.rowUpper_[row])) {
        rhs_value = lp.rowLower_[row];
        coefficient = -1;
      } else {
        rhs_value = lp.rowUpper_[row];
        coefficient = 1;
      }
      rhs[row] = rhs_value;

      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(coefficient);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(0);
      equality_lp.colUpper_.push_back(difference);
      equality_lp.colCost_.push_back(0);
    } else if (lp.rowLower_[row] == lp.rowUpper_[row]) {
      // equality row
      rhs[row] = lp.rowLower_[row];
    } else {
#ifdef HiGHSDEV
      printf(
          "transformIntoEqualityProblem: Unknown row type when adding slacks");
#endif
      return HighsStatus::Error;
    }
  }
  equality_lp.rowLower_ = rhs;
  equality_lp.rowUpper_ = rhs;

  return HighsStatus::OK;
}

// Given (P) returns (D) for the pair
// (P)
//    min c'x st Ax=b
//     st l <= x <= u
// (D)
//    max b'y + l'zl - u'zu
//     st A'y + zl - zu = c
//        y free, zl >=0, zu >= 0
HighsStatus dualizeEqualityProblem(const HighsLp& lp, HighsLp& dual) {
  std::vector<double> colCost = lp.colCost_;
  if (lp.sense_ != ObjSense::MINIMIZE) {
    for (int col = 0; col < lp.numCol_; col++) colCost[col] = -colCost[col];
  }

  assert(lp.rowLower_ == lp.rowUpper_);

  const int ncols = lp.numRow_;
  const int nrows = lp.numCol_;

  dual.numRow_ = nrows;
  dual.rowLower_ = colCost;
  dual.rowUpper_ = colCost;

  // Add columns (y)
  dual.numCol_ = ncols;
  dual.colLower_.resize(ncols);
  dual.colUpper_.resize(ncols);
  dual.colCost_.resize(ncols);

  for (int col = 0; col < ncols; col++) {
    dual.colLower_[col] = -HIGHS_CONST_INF;
    dual.colUpper_[col] = HIGHS_CONST_INF;
    // cost b'y
    dual.colCost_[col] = lp.rowLower_[col];
  }

  // Get transpose of A
  int i, k;
  vector<int> iwork(lp.numRow_, 0);
  dual.Astart_.resize(lp.numRow_ + 1, 0);
  int AcountX = lp.Aindex_.size();
  dual.Aindex_.resize(AcountX);
  dual.Avalue_.resize(AcountX);
  for (int k = 0; k < AcountX; k++) iwork.at(lp.Aindex_.at(k))++;
  for (i = 1; i <= lp.numRow_; i++)
    dual.Astart_.at(i) = dual.Astart_.at(i - 1) + iwork.at(i - 1);
  for (i = 0; i < lp.numRow_; i++) iwork.at(i) = dual.Astart_.at(i);
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    for (k = lp.Astart_.at(iCol); k < lp.Astart_.at(iCol + 1); k++) {
      int iRow = lp.Aindex_.at(k);
      int iPut = iwork.at(iRow)++;
      dual.Aindex_.at(iPut) = iCol;
      dual.Avalue_.at(iPut) = lp.Avalue_[k];
    }
  }

  // Add columns (zl)
  for (int col = 0; col < lp.numCol_; col++) {
    if (lp.colLower_[col] > -HIGHS_CONST_INF) {
      const int nnz = dual.Astart_[dual.numCol_];

      dual.colLower_.push_back(0);
      dual.colUpper_.push_back(HIGHS_CONST_INF);

      dual.colCost_.push_back(lp.colLower_[col]);

      // Add constaints
      dual.Astart_.push_back(nnz + 1);
      dual.Aindex_.push_back(col);
      dual.Avalue_.push_back(1.0);

      dual.numCol_++;
    }
  }

  // Add columns (zu)
  for (int col = 0; col < lp.numCol_; col++) {
    if (lp.colUpper_[col] < HIGHS_CONST_INF) {
      const int nnz = dual.Astart_[dual.numCol_];

      dual.colLower_.push_back(0);
      dual.colUpper_.push_back(HIGHS_CONST_INF);

      dual.colCost_.push_back(-lp.colUpper_[col]);

      // Add constaints
      dual.Astart_.push_back(nnz + 1);
      dual.Aindex_.push_back(col);
      dual.Avalue_.push_back(-1.0);

      dual.numCol_++;
    }
  }

  dual.sense_ = ObjSense::MINIMIZE;
  for (int col = 0; col < dual.numCol_; col++) {
    dual.colCost_[col] = -dual.colCost_[col];
  }

  dual.model_name_ = lp.model_name_ + "_dualized";

#ifdef HiGHSDEV
  printf("Dualized equality LP\n");
#endif

  return HighsStatus::OK;
}

void reportPresolveReductions(const HighsOptions& options, const HighsLp& lp,
                              const HighsLp& presolve_lp) {
  int num_col_from = lp.numCol_;
  int num_row_from = lp.numRow_;
  int num_els_from = lp.Astart_[num_col_from];
  int num_col_to = presolve_lp.numCol_;
  int num_row_to = presolve_lp.numRow_;
  int num_els_to;
  if (num_col_to) {
    num_els_to = presolve_lp.Astart_[num_col_to];
  } else {
    num_els_to = 0;
  }
  HighsPrintMessage(options.logfile, options.message_level, ML_ALWAYS,
                    "Presolve : Reductions: rows %d(-%d); columns %d(-%d); "
                    "elements %d(-%d)\n",
                    num_row_to, (num_row_from - num_row_to), num_col_to,
                    (num_col_from - num_col_to), num_els_to,
                    (num_els_from - num_els_to));
}

void reportPresolveReductions(const HighsOptions& options, const HighsLp& lp,
                              const bool presolve_to_empty) {
  int num_col_from = lp.numCol_;
  int num_row_from = lp.numRow_;
  int num_els_from = lp.Astart_[num_col_from];
  int num_col_to;
  int num_row_to;
  int num_els_to;
  std::string message;
  if (presolve_to_empty) {
    num_col_to = 0;
    num_row_to = 0;
    num_els_to = 0;
    message = "- Reduced to empty";
  } else {
    num_col_to = num_col_from;
    num_row_to = num_row_from;
    num_els_to = num_els_from;
    message = "- Not reduced";
  }
  HighsPrintMessage(options.logfile, options.message_level, ML_ALWAYS,
                    "Presolve : Reductions: rows %d(-%d); columns %d(-%d); "
                    "elements %d(-%d) %s\n",
                    num_row_to, (num_row_from - num_row_to), num_col_to,
                    (num_col_from - num_col_to), num_els_to,
                    (num_els_from - num_els_to), message.c_str());
}

bool isLessInfeasibleDSECandidate(const HighsOptions& options,
                                  const HighsLp& lp) {
  int max_col_num_en = -1;
  const int max_allowed_col_num_en = 24;
  const int max_assess_col_num_en = std::max(9, max_allowed_col_num_en);
  const int max_average_col_num_en = 6;
  vector<int> col_length_k;
  col_length_k.resize(1 + max_assess_col_num_en, 0);
  bool LiDSE_candidate = true;
  bool all_unit_nonzeros = true;
  for (int col = 0; col < lp.numCol_; col++) {
    // Check limit on number of entries in the column has not been breached
    int col_num_en = lp.Astart_[col + 1] - lp.Astart_[col];
    max_col_num_en = std::max(col_num_en, max_col_num_en);
    if (col_num_en > max_assess_col_num_en) {
#ifdef HiGHSDEV
      if (LiDSE_candidate)
        printf("Column %d has %d > %d entries so LP is not LiDSE candidate\n",
               col, col_num_en, max_allowed_col_num_en);
      LiDSE_candidate = false;
#else
      LiDSE_candidate = false;
      return LiDSE_candidate;
#endif
    } else {
      col_length_k[col_num_en]++;
    }
    for (int en = lp.Astart_[col]; en < lp.Astart_[col + 1]; en++) {
      double value = lp.Avalue_[en];
      // All nonzeros must be +1 or -1
      if (fabs(value) != 1) {
        all_unit_nonzeros = false;
#ifdef HiGHSDEV
        if (LiDSE_candidate)
          printf(
              "Column %d has entry %d with value %g so LP is not LiDSE "
              "candidate\n",
              col, en - lp.Astart_[col], value);
        LiDSE_candidate = false;
#else
        LiDSE_candidate = false;
        return LiDSE_candidate;
#endif
      }
    }
  }
#ifdef HiGHSDEV
  /*
  printf("LP has\n");
  int to_num_en = std::min(max_assess_col_num_en, max_col_num_en);
  for (int col_num_en = 0; col_num_en < to_num_en+1; col_num_en++)
    printf("%7d columns of count %1d\n", col_length_k[col_num_en], col_num_en);
  */
#endif
  double average_col_num_en = lp.Astart_[lp.numCol_];
  average_col_num_en = average_col_num_en / lp.numCol_;
  LiDSE_candidate =
      LiDSE_candidate && average_col_num_en <= max_average_col_num_en;
  std::string logic0 = "has";
  if (!all_unit_nonzeros) logic0 = "does not have";
  std::string logic1 = "is not";
  if (LiDSE_candidate) logic1 = "is";
  HighsLogMessage(
      options.logfile, HighsMessageType::INFO,
      "LP %s %s all |entries|=1; max column count = %d (limit %d); average "
      "column count = %0.2g (limit %d): So %s a candidate for LiDSE",
      lp.model_name_.c_str(), logic0.c_str(), max_col_num_en,
      max_allowed_col_num_en, average_col_num_en, max_average_col_num_en,
      logic1.c_str());
#ifdef HiGHSDEV
  int int_average_col_num_en = average_col_num_en;
  printf("grep_count_distrib,%s,%d,%d,%d\n", lp.model_name_.c_str(),
         max_col_num_en, int_average_col_num_en, LiDSE_candidate);
#endif
  return LiDSE_candidate;
}

void convertToMinimization(HighsLp& lp) {
  if (lp.sense_ != ObjSense::MINIMIZE) {
    for (int col = 0; col < lp.numCol_; col++)
      lp.colCost_[col] = -lp.colCost_[col];
  }
}

bool isEqualityProblem(const HighsLp& lp) {
  for (int row = 0; row < lp.numRow_; row++)
    if (lp.rowLower_[row] != lp.rowUpper_[row]) return false;

  return true;
}

double vectorProduct(const std::vector<double>& v1,
                     const std::vector<double>& v2) {
  assert(v1.size() == v2.size());
  double sum = 0;
  for (int i = 0; i < (int)v1.size(); i++) sum += v1[i] * v2[i];
  return sum;
}

HighsStatus calculateResidual(const HighsLp& lp, HighsSolution& solution,
                              std::vector<double>& residual) {
  HighsStatus status = calculateRowValues(lp, solution);
  if (status != HighsStatus::OK) return status;

  residual.clear();
  residual.resize(lp.numRow_);

  for (int row = 0; row < lp.numRow_; row++) {
    if (solution.row_value[row] < lp.rowLower_[row]) {
      residual[row] = lp.rowLower_[row] - solution.row_value[row];
    } else if (solution.row_value[row] > lp.rowUpper_[row]) {
      residual[row] = solution.row_value[row] - lp.rowUpper_[row];
    }
  }

  return status;
}

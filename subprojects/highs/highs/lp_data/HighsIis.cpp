/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsIis.cpp
 * @brief IIS utilities for HiGHS
 */

#include "Highs.h"

void HighsIis::clearData() {
  this->valid_ = false;
  this->status_ = kIisModelStatusUnknown;
  this->strategy_ = kIisStrategyMin;
  this->col_index_.clear();
  this->row_index_.clear();
  this->col_bound_.clear();
  this->row_bound_.clear();
  this->col_status_.clear();
  this->row_status_.clear();
  this->model_.clear();
}

void HighsIis::clear() {
  this->clearData();
  this->info_.clear();
}

void HighsIis::clearLogInfo() {
  this->info_.iis_last_disptime = -kHighsInf;
  this->info_.iis_num_disp_lines = 0;
}

void HighsIis::invalid(const HighsLp& lp) {
  this->clear();
  this->col_status_.assign(lp.num_col_, kIisStatusMaybeInConflict);
  this->row_status_.assign(lp.num_row_, kIisStatusMaybeInConflict);
}

std::string HighsIis::iisBoundStatusToString(HighsInt bound_status) const {
  if (bound_status == kIisBoundStatusDropped) return "Dropped";
  if (bound_status == kIisBoundStatusNull) return "   Null";
  if (bound_status == kIisBoundStatusFree) return "   Free";
  if (bound_status == kIisBoundStatusLower) return "  Lower";
  if (bound_status == kIisBoundStatusUpper) return "  Upper";
  if (bound_status == kIisBoundStatusBoxed) return "  Boxed";
  return "*****";
}

std::string HighsIis::iisModelStatusToString(HighsInt model_status) const {
  if (model_status == kIisModelStatusFeasible) return "Feasible";
  if (model_status == kIisModelStatusUnknown) return "Unknown";
  if (model_status == kIisModelStatusTimeLimit) return "Time limit reached";
  if (model_status == kIisModelStatusReducible) return "Reducible";
  if (model_status == kIisModelStatusIrreducible) return "Irreducible";
  return "*****";
}

void HighsIis::report(const std::string& message, const HighsLp& lp) const {
  HighsInt num_iis_col = this->col_index_.size();
  HighsInt num_iis_row = this->row_index_.size();
  if (num_iis_col > 10 || num_iis_row > 10) return;
  printf("\nIIS %s\n===\n", message.c_str());
  printf("Column: ");
  for (HighsInt iCol = 0; iCol < num_iis_col; iCol++) printf("%9d ", int(iCol));
  printf("\nStatus: ");
  for (HighsInt iCol = 0; iCol < num_iis_col; iCol++)
    printf("%9s ", iisBoundStatusToString(this->col_bound_[iCol]).c_str());
  printf("\nLower:  ");
  for (HighsInt iCol = 0; iCol < num_iis_col; iCol++)
    printf("%9.2g ", lp.col_lower_[iCol]);
  printf("\nUpper:  ");
  for (HighsInt iCol = 0; iCol < num_iis_col; iCol++)
    printf("%9.2g ", lp.col_upper_[iCol]);
  printf("\n");
  printf("Row:    Status     Lower     Upper\n");
  for (HighsInt iRow = 0; iRow < num_iis_row; iRow++)
    printf("%2d   %9s %9.2g %9.2g\n", int(iRow),
           iisBoundStatusToString(this->row_bound_[iRow]).c_str(),
           lp.row_lower_[iRow], lp.row_upper_[iRow]);
  printf("\n");
}

void HighsIis::reportIteration(const HighsOptions& options, const HighsInt iter,
                               const HighsInt num_rows_remaining,
                               const bool force) {
  const bool output_flag = *options.log_options.output_flag;
  if (!output_flag) return;
  const double min_interval = 5.0;
  const double runtime = info_.sum_simplex_times;
  if (!force && info_.iis_last_disptime > -0.5 * kHighsInf &&
      runtime - info_.iis_last_disptime < min_interval)
    return;
  // Update last time only when we actually print a line
  info_.iis_last_disptime = runtime;
  const int gap = 17;
  const int w_iter = int(strlen("Iteration")) + 2;
  const int w_rows = int(strlen("Rows")) + gap;
  const int w_time = int(strlen("Runtime")) + gap;

  // Print header every 20 lines (and on first line)
  if (info_.iis_num_disp_lines % 20 == 0) {
    highsLogUser(options.log_options, HighsLogType::kInfo, "%*s%*s%*s\n",
                 w_iter, "Iteration", w_rows, "Rows", w_time, "Runtime");
  }
  ++info_.iis_num_disp_lines;
  const std::string time_string = highsFormatToString("%.2fs", runtime);
  highsLogUser(options.log_options, HighsLogType::kInfo, "%*d%*d%*s\n", w_iter,
               int(iter), w_rows, int(num_rows_remaining), w_time,
               time_string.c_str());
}

void HighsIis::reportFinal(const HighsOptions& options) const {
  const bool output_flag = *options.log_options.output_flag;
  if (!output_flag) return;
  const HighsLogOptions& log_options = options.log_options;
  // Align colons by padding labels to the same width
  const int kLabelWidth = 19;
  highsLogUser(log_options, HighsLogType::kInfo, "\n");
  highsLogUser(log_options, HighsLogType::kInfo, "%-*s : %s\n", kLabelWidth,
               "IIS status", iisModelStatusToString(status_).c_str());
  highsLogUser(log_options, HighsLogType::kInfo, "%-*s : %d\n", kLabelWidth,
               "Rows", int(row_index_.size()));
  highsLogUser(log_options, HighsLogType::kInfo, "%-*s : %d\n", kLabelWidth,
               "Columns", int(col_index_.size()));
  highsLogUser(log_options, HighsLogType::kInfo, "%-*s : %.2fs\n", kLabelWidth,
               "HiGHS run time", info_.sum_simplex_times);
}

void HighsIis::addCol(const HighsInt col, const HighsInt status) {
  this->col_index_.push_back(col);
  this->col_bound_.push_back(status);
}

void HighsIis::addRow(const HighsInt row, const HighsInt status) {
  this->row_index_.push_back(row);
  this->row_bound_.push_back(status);
}

void HighsIis::removeCol(const HighsInt col) {
  HighsInt num_col = this->col_index_.size();
  assert(col < num_col);
  this->col_index_[col] = this->col_index_[num_col - 1];
  this->col_index_.resize(num_col - 1);
}

void HighsIis::removeRow(const HighsInt row) {
  HighsInt num_row = this->row_index_.size();
  assert(row < num_row);
  this->row_index_[row] = this->row_index_[num_row - 1];
  this->row_index_.resize(num_row - 1);
}

bool HighsIis::trivial(const HighsLp& lp, const HighsOptions& options) {
  this->clear();
  const bool col_priority = kIisStrategyColPriority & options.iis_strategy;
  for (HighsInt k = 0; k < 2; k++) {
    if ((col_priority && k == 0) || (!col_priority && k == 1)) {
      // Loop over columns first
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        if (lp.col_lower_[iCol] - lp.col_upper_[iCol] >
            2 * options.primal_feasibility_tolerance) {
          this->addCol(iCol, kIisBoundStatusBoxed);
          break;
        }
      }
      if (this->col_index_.size() > 0) break;
    } else {
      // Loop over rows first
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        if (lp.row_lower_[iRow] - lp.row_upper_[iRow] >
            2 * options.primal_feasibility_tolerance) {
          this->addRow(iRow, kIisBoundStatusBoxed);
          break;
        }
      }
      if (this->row_index_.size() > 0) break;
    }
  }
  HighsInt num_iis_col = this->col_index_.size();
  HighsInt num_iis_row = this->row_index_.size();
  // If one is found then we've found an IIS
  if (num_iis_col + num_iis_row > 0) {
    // Should have found exactly 1
    assert(num_iis_col + num_iis_row == 1);
    this->valid_ = true;
    this->status_ = kIisModelStatusIrreducible;
    this->strategy_ = options.iis_strategy;
    return true;
  }
  // Now look for empty rows that cannot have zero activity
  std::vector<HighsInt> count;
  // Get the row counts
  if (lp.a_matrix_.isColwise()) {
    count.assign(lp.num_row_, 0);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
           iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
        count[lp.a_matrix_.index_[iEl]]++;
    }
  } else {
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
      count.push_back(lp.a_matrix_.start_[iRow + 1] -
                      lp.a_matrix_.start_[iRow]);
  }
  assert(this->row_index_.size() == 0);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (count[iRow] > 0) continue;
    if (lp.row_lower_[iRow] > options.primal_feasibility_tolerance) {
      this->addRow(iRow, kIisBoundStatusLower);
    } else if (lp.row_upper_[iRow] < -options.primal_feasibility_tolerance) {
      this->addRow(iRow, kIisBoundStatusUpper);
    }
    if (this->row_index_.size() > 0) {
      // If one is found then we've found an IIS
      this->valid_ = true;
      this->status_ = kIisModelStatusIrreducible;
      this->strategy_ = options.iis_strategy;
      return true;
    }
  }
  return false;
}

bool HighsIis::rowValueBounds(const HighsLp& lp, const HighsOptions& options) {
  // Look for infeasible rows based on row value bounds
  this->clear();
  std::vector<double> lower_value;
  std::vector<double> upper_value;
  if (lp.a_matrix_.isColwise()) {
    lower_value.assign(lp.num_row_, 0);
    upper_value.assign(lp.num_row_, 0);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      const double lower = lp.col_lower_[iCol];
      const double upper = lp.col_upper_[iCol];
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
           iEl < lp.a_matrix_.start_[iCol + 1]; iEl++) {
        HighsInt iRow = lp.a_matrix_.index_[iEl];
        double value = lp.a_matrix_.value_[iEl];
        if (value > 0) {
          lower_value[iRow] += value * lower;
          upper_value[iRow] += value * upper;
        } else {
          lower_value[iRow] += value * upper;
          upper_value[iRow] += value * lower;
        }
      }
    }
  } else {
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
      double lower_row_value = 0;
      double upper_row_value = 0;
      for (HighsInt iEl = lp.a_matrix_.start_[iRow];
           iEl < lp.a_matrix_.start_[iRow + 1]; iEl++) {
        HighsInt iCol = lp.a_matrix_.index_[iEl];
        const double lower = lp.col_lower_[iCol];
        const double upper = lp.col_upper_[iCol];
        double value = lp.a_matrix_.value_[iEl];
        if (value > 0) {
          lower_row_value += value * lower;
          upper_row_value += value * upper;
        } else {
          lower_row_value += value * upper;
          upper_row_value += value * lower;
        }
      }
      lower_value.push_back(lower_row_value);
      upper_value.push_back(upper_row_value);
    }
  }
  bool below_lower = false;
  bool above_upper = false;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    below_lower = upper_value[iRow] <
                  lp.row_lower_[iRow] - options.primal_feasibility_tolerance;
    above_upper = lower_value[iRow] >
                  lp.row_upper_[iRow] + options.primal_feasibility_tolerance;
    if (below_lower || above_upper) {
      this->row_index_.push_back(iRow);
      if (below_lower) {
        this->row_bound_.push_back(kIisBoundStatusLower);
      } else {
        this->row_bound_.push_back(kIisBoundStatusUpper);
      }
      break;
    }
  }
  if (this->row_index_.size() == 0) {
    // Nothing found, but IIS data still valid
    this->clear();
    this->valid_ = true;
    this->status_ = kIisModelStatusUnknown;
    this->strategy_ = options.iis_strategy;
    return false;
  }
  assert(below_lower || above_upper);
  assert(!(below_lower && above_upper));
  // Found an infeasible row
  const HighsInt iRow = this->row_index_[0];
  const std::string row_name_string =
      lp.row_names_.size() > 0 ? "(" + lp.row_names_[iRow] + ")" : "";
  if (below_lower) {
    highsLogUser(
        options.log_options, HighsLogType::kInfo,
        "LP row %d %shas maximum row value of %g, below lower bound of %g\n",
        int(iRow), row_name_string.c_str(), upper_value[iRow],
        lp.row_lower_[iRow]);
  } else {
    highsLogUser(
        options.log_options, HighsLogType::kInfo,
        "LP row %d %shas minimum row value of %g, above upper bound of %g\n",
        int(iRow), row_name_string.c_str(), lower_value[iRow],
        lp.row_upper_[iRow]);
  }
  double value;
  auto setColBound = [&]() {
    if (below_lower) {
      if (value > 0) {
        this->col_bound_.push_back(kIisBoundStatusUpper);
      } else {
        this->col_bound_.push_back(kIisBoundStatusLower);
      }
    } else {
      if (value > 0) {
        this->col_bound_.push_back(kIisBoundStatusLower);
      } else {
        this->col_bound_.push_back(kIisBoundStatusUpper);
      }
    }
  };

  if (lp.a_matrix_.isColwise()) {
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
           iEl < lp.a_matrix_.start_[iCol + 1]; iEl++) {
        value = lp.a_matrix_.value_[iEl];
        if (lp.a_matrix_.index_[iEl] == iRow && value != 0) {
          this->col_index_.push_back(iCol);
          setColBound();
        }
      }
    }
  } else {
    for (HighsInt iEl = lp.a_matrix_.start_[iRow];
         iEl < lp.a_matrix_.start_[iRow + 1]; iEl++) {
      HighsInt iCol = lp.a_matrix_.index_[iEl];
      value = lp.a_matrix_.value_[iEl];
      if (value != 0) {
        this->col_index_.push_back(iCol);
        setColBound();
      }
    }
  }

  // There must be at least one column in the IIS
  assert(this->col_index_.size() > 0);
  assert(this->col_index_.size() == this->col_bound_.size());
  assert(this->row_index_.size() == this->row_bound_.size());
  this->valid_ = true;
  this->status_ = kIisModelStatusIrreducible;
  this->strategy_ = options.iis_strategy;
  return true;
}

HighsStatus HighsIis::deduce(const HighsLp& lp, const HighsOptions& options,
                             const HighsCallback& callback,
                             const HighsBasis& basis) {
  // The number of infeasible rows must be positive
  assert(this->row_index_.size() > 0);
  // Identify the LP corresponding to the set of infeasible rows
  std::vector<HighsInt> from_row = this->row_index_;
  std::vector<HighsInt> from_col;
  std::vector<HighsInt> to_row;
  to_row.assign(lp.num_row_, -1);
  // Only uses this->row_index_ to initialise from_row, so can clear
  this->clearData();
  // ToDo Exploit the known col_index_ and row_bound_ HighsIis
  // information
  //
  // To get the IIS data needs the matrix to be column-wise
  assert(lp.a_matrix_.isColwise());
  // Determine how to detect whether a row is in from_row and (then)
  // gather information about it
  for (HighsInt iX = 0; iX < HighsInt(from_row.size()); iX++)
    to_row[from_row[iX]] = iX;
  // Identify the columns (from_col) with nonzeros in the infeasible
  // rows
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    bool use_col = false;
    for (HighsInt iEl = lp.a_matrix_.start_[iCol];
         iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
      use_col = use_col || to_row[lp.a_matrix_.index_[iEl]] >= 0;
    if (use_col) from_col.push_back(iCol);
  }
  HighsInt to_num_col = from_col.size();
  HighsInt to_num_row = from_row.size();
  HighsLp to_lp;
  to_lp.num_col_ = to_num_col;
  to_lp.num_row_ = to_num_row;
  to_lp.a_matrix_.num_col_ = to_lp.num_col_;
  to_lp.a_matrix_.num_row_ = to_lp.num_row_;
  const bool has_col_names = lp.col_names_.size() > 0;
  for (HighsInt iCol = 0; iCol < to_num_col; iCol++) {
    to_lp.col_cost_.push_back(0);
    to_lp.col_lower_.push_back(lp.col_lower_[from_col[iCol]]);
    to_lp.col_upper_.push_back(lp.col_upper_[from_col[iCol]]);
    if (has_col_names)
      to_lp.col_names_.push_back(lp.col_names_[from_col[iCol]]);
    for (HighsInt iEl = lp.a_matrix_.start_[from_col[iCol]];
         iEl < lp.a_matrix_.start_[from_col[iCol] + 1]; iEl++) {
      HighsInt iRow = lp.a_matrix_.index_[iEl];
      if (to_row[iRow] >= 0) {
        to_lp.a_matrix_.index_.push_back(to_row[iRow]);
        to_lp.a_matrix_.value_.push_back(lp.a_matrix_.value_[iEl]);
      }
    }
    to_lp.a_matrix_.start_.push_back(to_lp.a_matrix_.index_.size());
  }
  const bool has_row_names = lp.row_names_.size() > 0;
  for (HighsInt iRow = 0; iRow < to_num_row; iRow++) {
    to_lp.row_lower_.push_back(lp.row_lower_[from_row[iRow]]);
    to_lp.row_upper_.push_back(lp.row_upper_[from_row[iRow]]);
    if (has_row_names)
      to_lp.row_names_.push_back(lp.row_names_[from_row[iRow]]);
  }
  HighsStatus return_status = this->compute(to_lp, options, callback);
  // Indirect the values into the original LP
  for (HighsInt& colindex : this->col_index_) colindex = from_col[colindex];
  for (HighsInt& rowindex : this->row_index_) rowindex = from_row[rowindex];
  if (kIisDevReport) this->report("On exit", lp);
  return return_status;
}

void HighsIis::setLp(const HighsLp& lp) {
  HighsLp& iis_lp = this->model_.lp_;
  iis_lp.clear();
  HighsInt iis_num_col = this->col_index_.size();
  HighsInt iis_num_row = this->row_index_.size();
  const bool colwise = lp.a_matrix_.isColwise();
  // Scatter the IIS rows (cols) into a full-length vector to identify
  // IIS rows (cols) with LP rows (cols) according to whether the
  // incumbent matrix is col-wise or row-wise
  std::vector<HighsInt> iis_row;
  std::vector<HighsInt> iis_col;
  if (colwise) {
    iis_row.assign(lp.num_row_, -1);
    for (HighsInt iisRow = 0; iisRow < iis_num_row; iisRow++)
      iis_row[this->row_index_[iisRow]] = iisRow;
  } else {
    iis_col.assign(lp.num_col_, -1);
    for (HighsInt iisCol = 0; iisCol < iis_num_col; iisCol++)
      iis_col[this->col_index_[iisCol]] = iisCol;
  }
  double bound;

  const bool has_row_name = lp.row_names_.size() > 0;
  for (HighsInt iisRow = 0; iisRow < iis_num_row; iisRow++) {
    HighsInt iRow = this->row_index_[iisRow];
    if (has_row_name) iis_lp.row_names_.push_back(lp.row_names_[iRow]);
    HighsInt row_bound = this->row_bound_[iisRow];
    assert(row_bound == kIisBoundStatusLower ||
           row_bound == kIisBoundStatusUpper ||
           row_bound == kIisBoundStatusBoxed);
    bound =
        row_bound == kIisBoundStatusLower || row_bound == kIisBoundStatusBoxed
            ? lp.row_lower_[iRow]
            : -kHighsInf;
    iis_lp.row_lower_.push_back(bound);
    bound =
        row_bound == kIisBoundStatusUpper || row_bound == kIisBoundStatusBoxed
            ? lp.row_upper_[iRow]
            : kHighsInf;
    iis_lp.row_upper_.push_back(bound);
    if (!colwise) {
      for (HighsInt iEl = lp.a_matrix_.start_[iRow];
           iEl < lp.a_matrix_.start_[iRow + 1]; iEl++) {
        HighsInt iCol = lp.a_matrix_.index_[iEl];
        HighsInt iisCol = iis_col[iCol];
        if (iisCol >= 0) {
          iis_lp.a_matrix_.index_.push_back(iisCol);
          iis_lp.a_matrix_.value_.push_back(lp.a_matrix_.value_[iEl]);
        }
      }
    }
  }

  const bool has_col_name = lp.col_names_.size() > 0;
  for (HighsInt iisCol = 0; iisCol < iis_num_col; iisCol++) {
    HighsInt iCol = this->col_index_[iisCol];
    // Costs in the IIS LP are zero since they play no role in IIS,
    // and when dropping bounds, optimality is the only valid model
    // status
    iis_lp.col_cost_.push_back(0);
    if (has_col_name) iis_lp.col_names_.push_back(lp.col_names_[iCol]);
    HighsInt col_bound = this->col_bound_[iisCol];
    assert(col_bound == kIisBoundStatusLower ||
           col_bound == kIisBoundStatusUpper ||
           col_bound == kIisBoundStatusBoxed ||
           col_bound == kIisBoundStatusFree);
    bound =
        col_bound == kIisBoundStatusLower || col_bound == kIisBoundStatusBoxed
            ? lp.col_lower_[iCol]
            : -kHighsInf;
    iis_lp.col_lower_.push_back(bound);
    bound =
        col_bound == kIisBoundStatusUpper || col_bound == kIisBoundStatusBoxed
            ? lp.col_upper_[iCol]
            : kHighsInf;
    iis_lp.col_upper_.push_back(bound);
    if (colwise) {
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
           iEl < lp.a_matrix_.start_[iCol + 1]; iEl++) {
        HighsInt iRow = lp.a_matrix_.index_[iEl];
        HighsInt iisRow = iis_row[iRow];
        if (iisRow >= 0) {
          iis_lp.a_matrix_.index_.push_back(iisRow);
          iis_lp.a_matrix_.value_.push_back(lp.a_matrix_.value_[iEl]);
        }
      }
    }
    iis_lp.a_matrix_.start_.push_back(iis_lp.a_matrix_.index_.size());
  }
  iis_lp.num_col_ = iis_lp.col_cost_.size();
  iis_lp.num_row_ = iis_lp.row_lower_.size();
  // The IIS LP matrix will have the same format as the incumbent LP
  iis_lp.a_matrix_.format_ = lp.a_matrix_.format_;
  iis_lp.a_matrix_.num_col_ = iis_lp.num_col_;
  iis_lp.a_matrix_.num_row_ = iis_lp.num_row_;
  iis_lp.model_name_ = lp.model_name_ + "_IIS";
}

HighsInt HighsIis::nonIsStatus() const {
  const bool is_feasible = this->status_ == kIisModelStatusFeasible;
  const bool has_is = this->col_index_.size() || this->row_index_.size();
  // If the model is known to be feasible, then there should be no IS,
  // and all columns and rows are kIisStatusNotInConflict
  if (is_feasible) assert(!has_is);
  if (has_is) assert(this->status_ >= kIisModelStatusTimeLimit);
  // If there is an IS, then all columns and rows not in the IS are
  // kIisStatusNotInConflict
  const HighsInt default_iis_status = is_feasible || has_is
                                          ? kIisStatusNotInConflict
                                          : kIisStatusMaybeInConflict;
  return default_iis_status;
}

void HighsIis::setStatus(const HighsLp& lp) {
  if (!this->valid_) return;
  const HighsInt non_is_status = nonIsStatus();
  const HighsInt in_is_status = this->status_ == kIisModelStatusIrreducible
                                    ? kIisStatusInConflict
                                    : kIisStatusMaybeInConflict;
  this->col_status_.assign(lp.num_col_, non_is_status);
  this->row_status_.assign(lp.num_row_, non_is_status);
  const HighsInt iis_num_col = this->col_index_.size();
  const HighsInt iis_num_row = this->row_index_.size();
  for (HighsInt iisCol = 0; iisCol < iis_num_col; iisCol++)
    this->col_status_[this->col_index_[iisCol]] = in_is_status;
  for (HighsInt iisRow = 0; iisRow < iis_num_row; iisRow++)
    this->row_status_[this->row_index_[iisRow]] = in_is_status;
}

HighsInt HighsIis::determineBoundStatus(const double lower, const double upper,
                                        const bool is_row) const {
  HighsInt iss_bound_status = kIisBoundStatusNull;
  if (lower <= -kHighsInf) {
    if (upper >= kHighsInf) {
      if (is_row) {
        // Free rows can be dropped
        iss_bound_status = kIisBoundStatusDropped;
      } else {
        // Free columns can only be dropped if they are empty
        iss_bound_status = kIisBoundStatusFree;
      }
    } else {
      iss_bound_status = kIisBoundStatusUpper;
    }
  } else {
    if (upper >= kHighsInf) {
      iss_bound_status = kIisBoundStatusLower;
    } else {
      // FX or BX
      iss_bound_status = kIisBoundStatusBoxed;
    }
  }
  assert(iss_bound_status != kIisBoundStatusNull);
  return iss_bound_status;
}

HighsStatus HighsIis::compute(const HighsLp& lp, const HighsOptions& options,
                              const HighsCallback& callback,
                              const HighsBasis* basis) {
  const HighsLogOptions& log_options = options.log_options;
  const bool col_priority = kIisStrategyColPriority & options.iis_strategy;
  const bool row_priority = !col_priority;
  // Initially all columns and rows are candidates for the IIS
  HighsInt num_rows = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    HighsInt col_status = this->determineBoundStatus(
        lp.col_lower_[iCol], lp.col_upper_[iCol], false);
    this->addCol(iCol, col_status);
  }
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    HighsInt row_status = this->determineBoundStatus(lp.row_lower_[iRow],
                                                     lp.row_upper_[iRow], true);
    this->addRow(iRow, row_status);
    if (row_status != kIisBoundStatusDropped) num_rows++;
  }

  Highs highs;
  highs.passOptions(options);
  highs.setOptionValue("output_flag", kIisDevReport);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue(
      "time_limit",
      std::max(options.iis_time_limit - this->info_.sum_simplex_times, 0.0));
  // Handle the callback propagation after setting output_flag false,
  // otherwise deprecation message for setLogCallback is echoed
  if (log_options.user_log_callback || callback.active[kCallbackLogging] ||
      callback.active[kCallbackSimplexInterrupt]) {
    // Setting the logging callbacks currently serves no purpose since
    // output_flag is set to kIisDevReport which is false (unless
    // developing) so that solves with this Highs instance are silent
    /*
    if (log_options.user_log_callback)
      highs.setLogCallback(log_options.user_log_callback);
    */
    highs.setCallback(callback.user_callback, callback.user_callback_data);
    if (callback.active[kCallbackLogging])
      highs.startCallback(kCallbackLogging);
    if (callback.active[kCallbackSimplexInterrupt])
      highs.startCallback(kCallbackSimplexInterrupt);
  }
  const HighsLp& incumbent_lp = highs.getLp();
  const HighsBasis& incumbent_basis = highs.getBasis();
  const HighsSolution& solution = highs.getSolution();
  const HighsInfo& info = highs.getInfo();
  HighsStatus run_status = highs.passModel(lp);
  assert(run_status == HighsStatus::kOk);
  if (basis) highs.setBasis(*basis);
  // Zero the objective
  std::vector<double> cost;
  cost.assign(lp.num_col_, 0);
  run_status = highs.changeColsCost(0, lp.num_col_ - 1, cost.data());
  assert(run_status == HighsStatus::kOk);
  // Solve the LP
  if (basis) highs.setBasis(*basis);
  const bool use_sensitivity_filter = false;
  std::vector<double> primal_phase1_dual;
  bool row_deletion = false;
  HighsInt iX = -1;
  bool drop_lower = false;
  // Lambda for gathering data when solving an LP
  auto solveLp = [&]() -> HighsStatus {
    double simplex_time = -highs.getRunTime();
    HighsInt simplex_iterations = -info.simplex_iteration_count;
    run_status = highs.optimizeModel();
    simplex_time += highs.getRunTime();
    simplex_iterations += info.simplex_iteration_count;
    this->info_.update(simplex_time, simplex_iterations);
    if (run_status != HighsStatus::kOk) return run_status;
    HighsModelStatus model_status = highs.getModelStatus();
    if (use_sensitivity_filter &&
        model_status == HighsModelStatus::kInfeasible) {
      printf("\nHighsIis::compute %s deletion for %d and %s bound\n",
             row_deletion ? "Row" : "Col", int(iX),
             drop_lower ? "Lower" : "Upper");
      bool output_flag;
      highs.getOptionValue("output_flag", output_flag);
      highs.setOptionValue("output_flag", true);
      HighsInt simplex_strategy;
      highs.getOptionValue("simplex_strategy", simplex_strategy);
      highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
      // Solve the LP
      run_status = highs.optimizeModel();
      if (run_status != HighsStatus::kOk) return run_status;
      highs.writeSolution("", kSolutionStylePretty);
      const HighsInt* basic_index = highs.getBasicVariablesArray();
      std::vector<double> rhs;
      rhs.assign(lp.num_row_, 0);
      // Get duals for nonbasic rows, and initialise duals so that basic duals
      // are zero
      assert(101 == 202);

      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        HighsInt iVar = basic_index[iRow];
        const double lower = iVar < lp.num_col_
                                 ? lp.col_lower_[iVar]
                                 : lp.row_lower_[iVar - lp.num_col_];
        const double upper = iVar < lp.num_col_
                                 ? lp.col_upper_[iVar]
                                 : lp.row_upper_[iVar - lp.num_col_];
        const double value = iVar < lp.num_col_
                                 ? solution.col_value[iVar]
                                 : solution.row_value[iVar - lp.num_col_];
        if (value < lower - options.primal_feasibility_tolerance) {
          rhs[iRow] = -1;
        } else if (value > upper + options.primal_feasibility_tolerance) {
          rhs[iRow] = 1;
        }
      }
      HVector pi;
      pi.setup(lp.num_row_);
      highs.getBasisTransposeSolve(rhs.data(), &pi.array[0], NULL, NULL);
      pi.count = lp.num_row_;
      std::vector<double> reduced_costs_value;
      std::vector<HighsInt> reduced_costs_index;
      lp.a_matrix_.productTransposeQuad(reduced_costs_value,
                                        reduced_costs_index, pi);

      primal_phase1_dual = highs.getPrimalPhase1Dual();
      HighsInt num_zero_dual = 0;
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        const HighsBasisStatus status = incumbent_basis.col_status[iCol];
        const double dual = primal_phase1_dual[iCol];
        const double lower = lp.col_lower_[iCol];
        const double upper = lp.col_upper_[iCol];
        const double value = solution.col_value[iCol];
        if (status != HighsBasisStatus::kBasic &&
            std::fabs(dual) < options.dual_feasibility_tolerance) {
          num_zero_dual++;
          // Small dual for nonbasic variable
          printf(
              "HighsIis::compute Column %d [%g, %g, %g] with status %s has "
              "dual %g\n",
              int(iCol), lower, value, upper,
              highs.basisStatusToString(status).c_str(), dual);
          //	  assert(123 == 456);
        }
      }
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        const HighsBasisStatus status = incumbent_basis.row_status[iRow];
        const double dual = primal_phase1_dual[lp.num_col_ + iRow];
        const double lower = lp.row_lower_[iRow];
        const double upper = lp.row_upper_[iRow];
        const double value = solution.row_value[iRow];
        if (status != HighsBasisStatus::kBasic &&
            std::fabs(dual) < options.dual_feasibility_tolerance) {
          num_zero_dual++;
          // Small dual for nonbasic variable
          printf(
              "HighsIis::compute Row    %d [%g, %g, %g] with status %s has "
              "dual %g\n",
              int(iRow), lower, value, upper,
              highs.basisStatusToString(status).c_str(), dual);
          //	  assert(123 == 456);
        }
      }
      highs.setOptionValue("output_flag", output_flag);
      highs.setOptionValue("simplex_strategy", simplex_strategy);
      assert(!num_zero_dual);
    }

    return run_status;
  };

  run_status = solveLp();
  // If we fail to establish infeasibility, return the initial subset
  if (run_status != HighsStatus::kOk) {
    this->valid_ = true;
    this->strategy_ = options.iis_strategy;
    if (highs.getModelStatus() == HighsModelStatus::kTimeLimit) {
      this->status_ = IisModelStatus::kIisModelStatusTimeLimit;
    } else {
      this->status_ = IisModelStatus::kIisModelStatusReducible;
    }
    return HighsStatus::kWarning;
  }

  assert(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  IisModelStatus iis_status = kIisModelStatusIrreducible;
  HighsStatus search_return_status = HighsStatus::kOk;

  // Initial logging
  this->clearLogInfo();
  HighsInt iter = 0;
  highsLogUser(log_options, HighsLogType::kInfo,
               "\nRunning deletion filter to identify an IIS\n");
  this->reportIteration(options, iter, num_rows, true);

  // Pass twice: rows before columns, or columns before rows, according to
  // row_priority
  for (HighsInt k = 0; k < 2; k++) {
    row_deletion = (row_priority && k == 0) || (!row_priority && k == 1);
    std::string type = row_deletion ? "Row" : "Col";
    // Perform deletion pass
    HighsInt num_index = row_deletion ? lp.num_row_ : lp.num_col_;
    for (iX = 0; iX < num_index; iX++) {
      // Get logging info
      iter++;
      const bool force = row_deletion && (iX == (num_index - 1));
      // Get status
      const HighsInt ix_status =
          row_deletion ? this->row_bound_[iX] : this->col_bound_[iX];
      // Skip if status is already free or dropped
      if (ix_status == kIisBoundStatusDropped ||
          ix_status == kIisBoundStatusFree) {
        // Possibly report
        this->reportIteration(options, iter, num_rows, force);
        continue;
      }

      double lower = row_deletion ? lp.row_lower_[iX] : lp.col_lower_[iX];
      double upper = row_deletion ? lp.row_upper_[iX] : lp.col_upper_[iX];
      // Record whether the upper bound has been dropped due to the lower
      // bound being kept
      if (lower > -kHighsInf && iis_status != kIisModelStatusTimeLimit) {
        // Drop the lower bound temporarily
        bool drop_lower = true;
        run_status = row_deletion
                         ? highs.changeRowBounds(iX, -kHighsInf, upper)
                         : highs.changeColBounds(iX, -kHighsInf, upper);
        assert(run_status == HighsStatus::kOk);
        // Solve the LP
        run_status = solveLp();
        this->processBoundRelaxation(highs, row_deletion, true, iX, lower,
                                     upper, iis_status, search_return_status);
      }
      if (upper < kHighsInf && iis_status != kIisModelStatusTimeLimit) {
        // Drop the upper bound temporarily
        run_status = row_deletion ? highs.changeRowBounds(iX, lower, kHighsInf)
                                  : highs.changeColBounds(iX, lower, kHighsInf);
        assert(run_status == HighsStatus::kOk);
        // Solve the LP
        run_status = solveLp();
        this->processBoundRelaxation(highs, row_deletion, false, iX, lower,
                                     upper, iis_status, search_return_status);
      }
      const bool debug_bound_change = true;
      if (debug_bound_change) {
        // Check bounds have been changed correctly
        double check_lower;
        double check_upper;
        double check_cost;
        HighsInt check_num_ix;
        HighsInt check_num_nz;
        run_status =
            row_deletion
                ? highs.getRows(iX, iX, check_num_ix, &check_lower,
                                &check_upper, check_num_nz, nullptr, nullptr,
                                nullptr)
                : highs.getCols(iX, iX, check_num_ix, &check_cost, &check_lower,
                                &check_upper, check_num_nz, nullptr, nullptr,
                                nullptr);
        assert(run_status == HighsStatus::kOk);
        assert(check_lower == lower);
        assert(check_upper == upper);
      }
      HighsInt iss_bound_status =
          this->determineBoundStatus(lower, upper, row_deletion);
      if (row_deletion) {
        this->row_bound_[iX] = iss_bound_status;
        if (iss_bound_status == kIisBoundStatusDropped) num_rows--;
      } else {
        this->col_bound_[iX] = iss_bound_status;
      }
      // Possibly report on iteration
      this->reportIteration(options, iter, num_rows, force);
      if (kIisDevReport) {
        highsLogUser(log_options, HighsLogType::kInfo, "%s %d has status %s\n",
                     type.c_str(), int(iX),
                     iisBoundStatusToString(iss_bound_status).c_str());
      }
    }
    if (k == 1) continue;
    // End of first pass: look to simplify second pass
    if (kIisDevReport) this->report("End of deletion", incumbent_lp);
    if (row_deletion) {
      // Mark empty columns as dropped
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        bool empty_col = true;
        for (HighsInt iEl = lp.a_matrix_.start_[iCol];
             iEl < lp.a_matrix_.start_[iCol + 1]; iEl++) {
          if (this->row_bound_[lp.a_matrix_.index_[iEl]] !=
              kIisBoundStatusDropped) {
            empty_col = false;
            break;
          }
        }
        if (empty_col) {
          if (kIisDevReport) {
            highsLogUser(log_options, HighsLogType::kInfo,
                         "Col %d has status Dropped: Empty\n", int(iCol));
          }
          this->col_bound_[iCol] = kIisBoundStatusDropped;
          run_status = highs.changeColBounds(iCol, -kHighsInf, kHighsInf);
          assert(run_status == HighsStatus::kOk);
        }
      }
    }
    if (kIisDevReport) this->report("End of pass 1", incumbent_lp);
  }
  if (kIisDevReport) this->report("End of pass 2", incumbent_lp);
  HighsInt iss_num_col = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (this->col_bound_[iCol] != kIisBoundStatusDropped) {
      this->col_index_[iss_num_col] = this->col_index_[iCol];
      this->col_bound_[iss_num_col] = this->col_bound_[iCol];
      iss_num_col++;
    }
  }
  HighsInt iss_num_row = 0;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (this->row_bound_[iRow] != kIisBoundStatusDropped) {
      this->row_index_[iss_num_row] = this->row_index_[iRow];
      this->row_bound_[iss_num_row] = this->row_bound_[iRow];
      iss_num_row++;
    }
  }
  // Return final result
  this->valid_ = true;
  this->status_ = iis_status;
  this->strategy_ = options.iis_strategy;
  this->col_index_.resize(iss_num_col);
  this->col_bound_.resize(iss_num_col);
  this->row_index_.resize(iss_num_row);
  this->row_bound_.resize(iss_num_row);
  return search_return_status;
}

void HighsIis::processBoundRelaxation(Highs& highs, const bool row_deletion,
                                      const bool drop_lower, const HighsInt iX,
                                      double& lower, double& upper,
                                      IisModelStatus& iis_status,
                                      HighsStatus& search_return_status) {
  // Declare run_status
  HighsStatus run_status;
  // Get model status
  HighsModelStatus model_status = highs.getModelStatus();
  if (model_status == HighsModelStatus::kOptimal) {
    // Now feasible, so restore the bound
    run_status = row_deletion ? highs.changeRowBounds(iX, lower, upper)
                              : highs.changeColBounds(iX, lower, upper);
    assert(run_status == HighsStatus::kOk);
    // If the lower bound must be kept, then any finite upper bound
    // must be dropped
    const bool apply_reciprocal_rule = drop_lower && (upper < kHighsInf);
    if (apply_reciprocal_rule) {
      // Drop the upper bound permanently
      upper = kHighsInf;
      run_status = row_deletion ? highs.changeRowBounds(iX, lower, upper)
                                : highs.changeColBounds(iX, lower, upper);
      assert(run_status == HighsStatus::kOk);
    }
    if (drop_lower) {
      assert(upper >= kHighsInf);
    }
    // Since upper = kHighsInf, allow the loop to run so that
    // bound status is set as if upper were set to kHighsInf
    // by relaxing it and finding that the LP was still
    // infeasible
  } else if (model_status == HighsModelStatus::kInfeasible) {
    // Bound can be dropped permanently
    if (drop_lower) {
      lower = -kHighsInf;
    } else {
      upper = kHighsInf;
    }
  } else if (model_status == HighsModelStatus::kTimeLimit) {
    // Time limit reached, so restore the bound and set iis status to
    // kIisModelStatusTimeLimit
    run_status = row_deletion ? highs.changeRowBounds(iX, lower, upper)
                              : highs.changeColBounds(iX, lower, upper);
    assert(run_status == HighsStatus::kOk);
    iis_status = kIisModelStatusTimeLimit;
    search_return_status = HighsStatus::kWarning;
  } else {
    // Unknown failure, so restore the bound and set iis status to
    // kIisModelStatusReducible
    run_status = row_deletion ? highs.changeRowBounds(iX, lower, upper)
                              : highs.changeColBounds(iX, lower, upper);
    assert(run_status == HighsStatus::kOk);
    iis_status = kIisModelStatusReducible;
    search_return_status = HighsStatus::kWarning;
  }
}

bool HighsIis::indexStatusOk(const HighsLp& lp) const {
  HighsInt num_col = lp.num_col_;
  HighsInt num_row = lp.num_row_;
  bool col_status_size_ok =
      this->col_status_.size() == static_cast<size_t>(num_col);
  bool row_status_size_ok =
      this->row_status_.size() == static_cast<size_t>(num_row);
  assert(col_status_size_ok);
  assert(row_status_size_ok);
  if (!col_status_size_ok) return indexStatusOkReturn(false);
  if (!row_status_size_ok) return indexStatusOkReturn(false);
  HighsInt num_iis_col = this->col_index_.size();
  HighsInt num_iis_row = this->row_index_.size();
  // Determine whether this is an IIS or just an IS
  bool true_iis = false;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (this->col_status_[iCol] == kIisStatusInConflict) {
      true_iis = true;
      break;
    }
  }
  if (!true_iis) {
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (this->row_status_[iRow] == kIisStatusInConflict) {
        true_iis = true;
        break;
      }
    }
  }
  // Now check that cols and rows in the IIS are kIisStatusInConflict
  // or kIisStatusMaybeInConflict, according to true_iis, and that all
  // other cols and rows are kIisStatusNotConflict
  std::vector<HighsInt> col_status = col_status_;
  std::vector<HighsInt> row_status = row_status_;
  const HighsInt illegal_status = -99;
  for (HighsInt iX = 0; iX < num_iis_col; iX++) {
    HighsInt iCol = this->col_index_[iX];
    if (col_status_[iCol] !=
        (true_iis ? kIisStatusInConflict : kIisStatusMaybeInConflict))
      return indexStatusOkReturn(false);
    col_status[iCol] = illegal_status;
  }
  for (HighsInt iX = 0; iX < num_iis_row; iX++) {
    HighsInt iRow = this->row_index_[iX];
    if (row_status_[iRow] !=
        (true_iis ? kIisStatusInConflict : kIisStatusMaybeInConflict))
      return indexStatusOkReturn(false);
    row_status[iRow] = illegal_status;
  }
  const HighsInt non_is_status = nonIsStatus();
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (col_status[iCol] > illegal_status && col_status[iCol] != non_is_status)
      return indexStatusOkReturn(false);
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (row_status[iRow] > illegal_status && row_status[iRow] != non_is_status)
      return indexStatusOkReturn(false);
  }
  return indexStatusOkReturn(true);
}

bool HighsIis::lpDataOk(const HighsLp& lp, const HighsOptions& options) const {
  const HighsLp& iis_lp = this->model_.lp_;
  HighsInt iis_num_col = this->col_index_.size();
  HighsInt iis_num_row = this->row_index_.size();
  if (!(iis_lp.num_col_ == iis_num_col)) return lpDataOkReturn(false);
  if (!(iis_lp.num_row_ == iis_num_row)) return lpDataOkReturn(false);

  const bool colwise = lp.a_matrix_.isColwise();

  const HighsInt illegal_index = -1;
  const double illegal_value = kHighsInf;
  // iis_row/col give the row/col in the IIS for each row/col in the
  // LP, or an illegal index if the LP row/col isn't in the IIS
  std::vector<HighsInt> iis_row;
  iis_row.assign(lp.num_row_, illegal_index);
  std::vector<HighsInt> iis_col;
  iis_col.assign(lp.num_col_, illegal_index);
  double bound;
  for (HighsInt iisRow = 0; iisRow < iis_num_row; iisRow++) {
    HighsInt iRow = this->row_index_[iisRow];
    iis_row[iRow] = iisRow;
    HighsInt row_bound = this->row_bound_[iisRow];
    bound =
        row_bound == kIisBoundStatusLower || row_bound == kIisBoundStatusBoxed
            ? lp.row_lower_[iRow]
            : -kHighsInf;
    if (iis_lp.row_lower_[iisRow] != bound) return lpDataOkReturn(false);
    bound =
        row_bound == kIisBoundStatusUpper || row_bound == kIisBoundStatusBoxed
            ? lp.row_upper_[iRow]
            : kHighsInf;
    if (iis_lp.row_upper_[iisRow] != bound) return lpDataOkReturn(false);
  }

  // Work through the LP columns checking the zero costs and bounds
  for (HighsInt iisCol = 0; iisCol < iis_num_col; iisCol++) {
    HighsInt iCol = this->col_index_[iisCol];
    iis_col[iCol] = iisCol;
    if (iis_lp.col_cost_[iisCol]) return lpDataOkReturn(false);
    HighsInt col_bound = this->col_bound_[iisCol];
    bound =
        col_bound == kIisBoundStatusLower || col_bound == kIisBoundStatusBoxed
            ? lp.col_lower_[iCol]
            : -kHighsInf;
    if (iis_lp.col_lower_[iisCol] != bound) return lpDataOkReturn(false);
    bound =
        col_bound == kIisBoundStatusUpper || col_bound == kIisBoundStatusBoxed
            ? lp.col_upper_[iCol]
            : kHighsInf;
    if (iis_lp.col_upper_[iisCol] != bound) return lpDataOkReturn(false);
  }
  std::vector<HighsInt> index;
  std::vector<double> value;
  // Work through the LP matrix, checking the matrix index/value
  if (colwise) {
    for (HighsInt iisCol = 0; iisCol < iis_num_col; iisCol++) {
      HighsInt iCol = this->col_index_[iisCol];
      // Use index/value to scatter the IIS matrix column
      index.assign(iis_num_row, illegal_index);
      value.assign(iis_num_row, illegal_value);
      for (HighsInt iEl = iis_lp.a_matrix_.start_[iisCol];
           iEl < iis_lp.a_matrix_.start_[iisCol + 1]; iEl++) {
        HighsInt iisRow = iis_lp.a_matrix_.index_[iEl];
        HighsInt iRow = this->row_index_[iisRow];
        index[iisRow] = iRow;
        value[iisRow] = iis_lp.a_matrix_.value_[iEl];
      }
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
           iEl < lp.a_matrix_.start_[iCol + 1]; iEl++) {
        HighsInt iRow = lp.a_matrix_.index_[iEl];
        HighsInt iisRow = iis_row[iRow];
        if (iisRow >= 0) {
          if (index[iisRow] != iRow) return lpDataOkReturn(false);
          if (value[iisRow] != lp.a_matrix_.value_[iEl])
            return lpDataOkReturn(false);
          index[iisRow] = illegal_index;
          value[iisRow] = illegal_value;
        }
      }
    }
  } else {
    for (HighsInt iisRow = 0; iisRow < iis_num_row; iisRow++) {
      HighsInt iRow = this->row_index_[iisRow];
      // Use index/value to scatter the IIS matrix row
      index.assign(iis_num_col, illegal_index);
      value.assign(iis_num_col, illegal_value);
      for (HighsInt iEl = iis_lp.a_matrix_.start_[iisRow];
           iEl < iis_lp.a_matrix_.start_[iisRow + 1]; iEl++) {
        HighsInt iisCol = iis_lp.a_matrix_.index_[iEl];
        HighsInt iCol = this->col_index_[iisCol];
        index[iisCol] = iCol;
        value[iisCol] = iis_lp.a_matrix_.value_[iEl];
      }
      for (HighsInt iEl = lp.a_matrix_.start_[iRow];
           iEl < lp.a_matrix_.start_[iRow + 1]; iEl++) {
        HighsInt iCol = lp.a_matrix_.index_[iEl];
        HighsInt iisCol = iis_col[iCol];
        if (iisCol >= 0) {
          if (index[iisCol] != iCol) return lpDataOkReturn(false);
          if (value[iisCol] != lp.a_matrix_.value_[iEl])
            return lpDataOkReturn(false);
          index[iisCol] = illegal_index;
          value[iisCol] = illegal_value;
        }
      }
    }
  }
  // Work through the IIS LP matrix, making sure that the index/value
  // are correct
  if (colwise) {
    for (HighsInt iisCol = 0; iisCol < iis_num_col; iisCol++) {
      HighsInt iCol = this->col_index_[iisCol];
      // Use index/value to scatter the LP matrix column
      index.assign(lp.num_row_, illegal_index);
      value.assign(lp.num_row_, illegal_value);
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
           iEl < lp.a_matrix_.start_[iCol + 1]; iEl++) {
        HighsInt iRow = lp.a_matrix_.index_[iEl];
        HighsInt iisRow = iis_row[iRow];
        index[iRow] = iisRow;
        value[iRow] = lp.a_matrix_.value_[iEl];
      }
      for (HighsInt iEl = iis_lp.a_matrix_.start_[iisCol];
           iEl < iis_lp.a_matrix_.start_[iisCol + 1]; iEl++) {
        HighsInt iisRow = iis_lp.a_matrix_.index_[iEl];
        HighsInt iRow = this->row_index_[iisRow];
        if (index[iRow] != iisRow) return lpDataOkReturn(false);
        if (value[iRow] != iis_lp.a_matrix_.value_[iEl])
          return lpDataOkReturn(false);
      }
    }
  } else {
    for (HighsInt iisRow = 0; iisRow < iis_num_row; iisRow++) {
      HighsInt iRow = this->row_index_[iisRow];
      // Use index/value to scatter the LP matrix row
      index.assign(lp.num_col_, illegal_index);
      value.assign(lp.num_col_, illegal_value);
      for (HighsInt iEl = lp.a_matrix_.start_[iRow];
           iEl < lp.a_matrix_.start_[iRow + 1]; iEl++) {
        HighsInt iCol = lp.a_matrix_.index_[iEl];
        HighsInt iisCol = iis_col[iCol];
        index[iCol] = iisCol;
        value[iCol] = lp.a_matrix_.value_[iEl];
      }
      for (HighsInt iEl = iis_lp.a_matrix_.start_[iisRow];
           iEl < iis_lp.a_matrix_.start_[iisRow + 1]; iEl++) {
        HighsInt iisCol = iis_lp.a_matrix_.index_[iEl];
        HighsInt iCol = this->col_index_[iisCol];
        if (index[iCol] != iisCol) return lpDataOkReturn(false);
        if (value[iCol] != iis_lp.a_matrix_.value_[iEl])
          return lpDataOkReturn(false);
      }
    }
  }
  return lpDataOkReturn(true);
}

bool HighsIis::lpOk(const HighsOptions& options) const {
  // Check that the IIS LP is OK (infeasible and optimal if
  // any bound is relaxed)
  if (!this->valid_) return lpOkReturn(true);
  HighsInt num_iis_col = this->col_index_.size();
  HighsInt num_iis_row = this->row_index_.size();
  // If an LP has a row with inconsistent bounds, or an empty row with
  // a positive lower bound or negative upper bound, then it is
  // infeasible, but the IIS contains no columns
  if (num_iis_col == 0) return true;
  const HighsLogOptions& log_options = options.log_options;
  const HighsLp& iis_lp = this->model_.lp_;
  assert(iis_lp.num_col_ == num_iis_col);
  assert(iis_lp.num_row_ == num_iis_row);
  Highs h;
  h.passOptions(options);
  h.setOptionValue("output_flag", false);
  h.passModel(iis_lp);
  h.writeModel("");
  HighsStatus status = h.optimizeModel();
  if (status != HighsStatus::kOk ||
      h.getModelStatus() != HighsModelStatus::kInfeasible) {
    highsLogUser(log_options, HighsLogType::kWarning,
                 "HighsIis: Failed to prove infeasibility for IIS LP\n");
    return lpOkReturn(false);
  }
  if (!(this->status_ == kIisModelStatusIrreducible)) return lpOkReturn(true);
  auto optimal = [&]() -> bool {
    if (options.log_dev_level > 0) h.writeModel("");
    h.optimizeModel();
    return h.getModelStatus() == HighsModelStatus::kOptimal;
  };
  for (HighsInt iisCol = 0; iisCol < num_iis_col; iisCol++) {
    HighsInt iCol = this->col_index_[iisCol];
    if (this->col_bound_[iisCol] == kIisBoundStatusLower) {
      h.changeColBounds(iisCol, -kHighsInf, iis_lp.col_upper_[iisCol]);
      if (!optimal()) {
        highsLogUser(log_options, HighsLogType::kWarning,
                     "HighsIis: IIS column %d (LP column %d): relaxing lower "
                     "bound of %g yield IIS LP with status %s\n",
                     int(iisCol), int(iCol), iis_lp.col_lower_[iisCol],
                     h.modelStatusToString(h.getModelStatus()).c_str());
        return lpOkReturn(false);
      }
      h.changeColBounds(iisCol, iis_lp.col_lower_[iisCol],
                        iis_lp.col_upper_[iisCol]);
    }
    if (this->col_bound_[iisCol] == kIisBoundStatusUpper) {
      h.changeColBounds(iisCol, iis_lp.col_lower_[iisCol], kHighsInf);
      if (!optimal()) {
        highsLogUser(log_options, HighsLogType::kWarning,
                     "HighsIis: IIS column %d (LP column %d): relaxing upper "
                     "bound of %g yield IIS LP with status %s\n",
                     int(iisCol), int(iCol), iis_lp.col_upper_[iisCol],
                     h.modelStatusToString(h.getModelStatus()).c_str());
        return lpOkReturn(false);
      }
      h.changeColBounds(iisCol, iis_lp.col_lower_[iisCol],
                        iis_lp.col_upper_[iisCol]);
    }
  }
  for (HighsInt iisRow = 0; iisRow < num_iis_row; iisRow++) {
    HighsInt iRow = this->row_index_[iisRow];
    if (this->row_bound_[iisRow] == kIisBoundStatusLower) {
      h.changeRowBounds(iisRow, -kHighsInf, iis_lp.row_upper_[iisRow]);
      if (!optimal()) {
        highsLogUser(log_options, HighsLogType::kError,
                     "HighsIis: IIS row %d (LP row %d): relaxing lower bound "
                     "of %g yield IIS LP with status %s\n",
                     int(iisRow), int(iRow), iis_lp.row_lower_[iisRow],
                     h.modelStatusToString(h.getModelStatus()).c_str());
        return lpOkReturn(false);
      }
      h.changeRowBounds(iisRow, iis_lp.row_lower_[iisRow],
                        iis_lp.row_upper_[iisRow]);
    }
    if (this->row_bound_[iisRow] == kIisBoundStatusUpper) {
      h.changeRowBounds(iisRow, iis_lp.row_lower_[iisRow], kHighsInf);
      if (!optimal()) {
        highsLogUser(log_options, HighsLogType::kError,
                     "HighsIis: IIS row %d (LP row %d): relaxing upper "
                     "bound of %g yield IIS LP with status %s\n",
                     int(iisRow), int(iRow), iis_lp.row_upper_[iisRow],
                     h.modelStatusToString(h.getModelStatus()).c_str());
        return lpOkReturn(false);
      }
      h.changeRowBounds(iisRow, iis_lp.row_lower_[iisRow],
                        iis_lp.row_upper_[iisRow]);
    }
  }
  return lpOkReturn(true);
}

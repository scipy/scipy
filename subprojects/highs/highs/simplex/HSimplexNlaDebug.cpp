/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSimplexNlaDebug.cpp
 *
 * @brief Debugging code for simplex NLA
 */
#include "simplex/HSimplexNla.h"
#include "util/HighsRandom.h"

// #include <stdio.h>

const double kResidualLargeError = 1e-8;
const double kResidualExcessiveError = sqrt(kResidualLargeError);

const double kSolveLargeError = 1e-8;
const double kSolveExcessiveError = sqrt(kSolveLargeError);

HighsDebugStatus HSimplexNla::debugCheckInvert(
    const std::string& message, const HighsInt alt_debug_level) const {
  // Sometimes a value other than highs_debug_level is passed as
  // alt_debug_level, either to force debugging, or to limit
  // debugging. If no value is passed, then alt_debug_level = -1, and
  // highs_debug_level is used
  const HighsInt use_debug_level = alt_debug_level >= 0
                                       ? alt_debug_level
                                       : this->options_->highs_debug_level;
  if (use_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  // If highs_debug_level isn't being used, then indicate that it's
  // being forced, and also force reporting of OK errors
  const bool force = alt_debug_level > this->options_->highs_debug_level;
  if (force)
    highsLogDev(this->options_->log_options, HighsLogType::kInfo,
                "CheckNlaINVERT:   Forcing debug\n");

  HighsDebugStatus return_status = HighsDebugStatus::kNotChecked;
  return_status = HighsDebugStatus::kOk;

  const HighsInt num_row = this->lp_->num_row_;
  const HighsInt num_col = this->lp_->num_col_;
  const vector<HighsInt>& a_matrix_start = this->lp_->a_matrix_.start_;
  const vector<HighsInt>& a_matrix_index = this->lp_->a_matrix_.index_;
  const vector<double>& a_matrix_value = this->lp_->a_matrix_.value_;
  const HighsInt* basic_index = this->basic_index_;
  const HighsOptions* options = this->options_;
  // Make sure that this isn't called between the matrix and LP resizing
  assert(num_row == this->lp_->a_matrix_.num_row_);
  assert(num_col == this->lp_->a_matrix_.num_col_);
  const bool report = (options->log_dev_level != 0);

  highsLogDev(options->log_options, HighsLogType::kInfo, "\nCheckINVERT: %s\n",
              message.c_str());

  HVector column;
  HVector rhs;
  column.setup(num_row);
  rhs.setup(num_row);
  HVector residual;
  double expected_density = 1;

  const bool random_ftran_test = true;
  const bool report_basis = options->log_dev_level > 1 && num_row < 20;
  if (random_ftran_test) {
    // Solve for a random solution
    HighsRandom random(1);
    // Solve Bx=b
    column.clear();
    rhs.clear();
    column.count = -1;
    if (report_basis)
      highsLogDev(options_->log_options, HighsLogType::kInfo, "Basis:");
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      rhs.index[rhs.count++] = iRow;
      double value = random.fraction();
      column.array[iRow] = value;
      HighsInt iCol = basic_index[iRow];
      if (report_basis)
        highsLogDev(options_->log_options, HighsLogType::kInfo, " %1d",
                    (int)iCol);
      if (iCol < num_col) {
        for (HighsInt iEl = a_matrix_start[iCol];
             iEl < a_matrix_start[iCol + 1]; iEl++) {
          HighsInt index = a_matrix_index[iEl];
          rhs.array[index] += value * a_matrix_value[iEl];
        }
      } else {
        HighsInt index = iCol - num_col;
        assert(index < num_row);
        rhs.array[index] += value;
      }
    }
    if (report_basis)
      highsLogDev(options_->log_options, HighsLogType::kInfo, "\n");
    residual = rhs;
    //  if (options->log_dev_level) {
    //    reportArray("Random solution before FTRAN", &column, true);
    //    reportArray("Random RHS      before FTRAN", &rhs, true);
    //  }
    this->ftran(rhs, expected_density);

    return_status =
        debugReportInvertSolutionError(false, column, rhs, residual, force);
  }

  const bool random_btran_test = true;
  if (random_btran_test) {
    // Solve B^Tx=b
    rhs.clear();
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      rhs.index[rhs.count++] = iRow;
      HighsInt iCol = basic_index[iRow];
      if (iCol < num_col) {
        for (HighsInt iEl = a_matrix_start[iCol];
             iEl < a_matrix_start[iCol + 1]; iEl++) {
          rhs.array[iRow] +=
              column.array[a_matrix_index[iEl]] * a_matrix_value[iEl];
        }
      } else {
        rhs.array[iRow] += column.array[iCol - num_col];
      }
    }
    residual = rhs;
    this->btran(rhs, expected_density);

    return_status =
        debugReportInvertSolutionError(true, column, rhs, residual, force);
  }

  if (use_debug_level < kHighsDebugLevelExpensive) return return_status;

  std::string value_adjective;
  expected_density = 0;
  double inverse_error_norm = 0;
  double residual_error_norm = 0;
  const bool test_inverse_ftran = true;
  if (test_inverse_ftran) {
    // Solve BX=B
    HighsInt check_col = -1;
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      HighsInt iCol = basic_index[iRow];
      column.clear();
      column.packFlag = true;
      if (iCol < num_col) {
        for (HighsInt k = a_matrix_start[iCol]; k < a_matrix_start[iCol + 1];
             k++) {
          HighsInt index = a_matrix_index[k];
          column.array[index] = a_matrix_value[k];
          column.index[column.count++] = index;
        }
      } else {
        HighsInt index = iCol - num_col;
        column.array[index] = 1.0;
        column.index[column.count++] = index;
      }
      const bool report_col = report && iRow == check_col;
      if (report_col) {
        reportArray("Check col before FTRAN", &column, true);
        //      factor_.reportLu(kReportLuBoth, true);
      }
      HVector residual = column;
      this->ftran(column, expected_density);
      if (report_col) reportArray("Check col after  FTRAN", &column, true);
      double inverse_col_error_norm = 0;
      for (HighsInt lc_iRow = 0; lc_iRow < num_row; lc_iRow++) {
        double ckValue = lc_iRow == iRow ? 1 : 0;
        double inverse_error = fabs(column.array[lc_iRow] - ckValue);
        inverse_col_error_norm =
            std::max(inverse_error, inverse_col_error_norm);
      }
      // Extra printing of intermediate errors
      if (report_col)
        highsLogDev(
            options->log_options, HighsLogType::kInfo,
            "CheckINVERT: Basic column %2d = %2d has inverse error %11.4g\n",
            (int)iRow, int(iCol), inverse_col_error_norm);
      inverse_error_norm = std::max(inverse_col_error_norm, inverse_error_norm);
      double residual_col_error_norm =
          debugInvertResidualError(false, column, residual);
      residual_error_norm =
          std::max(residual_col_error_norm, residual_error_norm);
    }
    return_status = debugReportInvertSolutionError(
        "inverse", false, inverse_error_norm, residual_error_norm, force);
  }

  const bool test_inverse_btran = true;
  if (test_inverse_btran) {
    // Solve B^TX=B^T
    HighsInt check_row = -2;
    inverse_error_norm = 0;
    residual_error_norm = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      column.clear();
      column.packFlag = true;
      for (HighsInt iCol = 0; iCol < num_row; iCol++) {
        HighsInt iVar = basic_index[iCol];
        if (iVar < num_col) {
          for (HighsInt k = a_matrix_start[iVar]; k < a_matrix_start[iVar + 1];
               k++) {
            if (a_matrix_index[k] == iRow) {
              column.array[iCol] = a_matrix_value[k];
              break;
            }
          }
        } else {
          if (iVar == num_col + iRow) column.array[iCol] = 1.0;
        }
      }
      for (HighsInt iCol = 0; iCol < num_row; iCol++) {
        if (column.array[iCol]) column.index[column.count++] = iCol;
      }
      const bool report_row = report && iRow == check_row;
      if (report_row) {
        reportArray("Check col before BTRAN", &column, true);
        //      factor_.reportLu(kReportLuBoth, true);
      }
      HVector residual = column;
      this->btran(column, expected_density);
      if (report_row) reportArray("Check col after  BTRAN", &column, true);
      double inverse_col_error_norm = 0;
      for (HighsInt lc_iRow = 0; lc_iRow < num_row; lc_iRow++) {
        double value = column.array[lc_iRow];
        double ckValue = lc_iRow == iRow ? 1 : 0;
        double inverse_error = fabs(value - ckValue);
        inverse_col_error_norm =
            std::max(inverse_error, inverse_col_error_norm);
      }
      // Extra printing of intermediate errors
      if (report_row)
        highsLogDev(
            options->log_options, HighsLogType::kInfo,
            "CheckINVERT: Basis matrix row %2d has inverse error %11.4g\n",
            (int)iRow, inverse_col_error_norm);
      inverse_error_norm = std::max(inverse_col_error_norm, inverse_error_norm);
      double residual_col_error_norm =
          debugInvertResidualError(true, column, residual);
      residual_error_norm =
          std::max(residual_col_error_norm, residual_error_norm);
    }
    return_status = debugReportInvertSolutionError(
        "inverse", true, inverse_error_norm, residual_error_norm, force);
  }

  return return_status;
}

double HSimplexNla::debugInvertResidualError(const bool transposed,
                                             const HVector& solution,
                                             HVector& residual) const {
  const HighsInt num_row = this->lp_->num_row_;
  const HighsInt num_col = this->lp_->num_col_;
  const vector<HighsInt>& a_matrix_start = this->lp_->a_matrix_.start_;
  const vector<HighsInt>& a_matrix_index = this->lp_->a_matrix_.index_;
  const vector<double>& a_matrix_value = this->lp_->a_matrix_.value_;
  const HighsInt* basic_index = this->basic_index_;

  if (transposed) {
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      HighsInt iCol = basic_index[iRow];
      if (iCol < num_col) {
        for (HighsInt iEl = a_matrix_start[iCol];
             iEl < a_matrix_start[iCol + 1]; iEl++) {
          HighsInt index = a_matrix_index[iEl];
          double value = solution.array[index];
          residual.array[iRow] -= value * a_matrix_value[iEl];
        }
      } else {
        HighsInt index = iCol - num_col;
        double value = solution.array[index];
        residual.array[iRow] -= value;
      }
    }
  } else {
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      double value = solution.array[iRow];
      HighsInt iCol = basic_index[iRow];
      if (iCol < num_col) {
        for (HighsInt iEl = a_matrix_start[iCol];
             iEl < a_matrix_start[iCol + 1]; iEl++) {
          HighsInt index = a_matrix_index[iEl];
          residual.array[index] -= value * a_matrix_value[iEl];
        }
      } else {
        HighsInt index = iCol - num_col;
        residual.array[index] -= value;
      }
    }
  }

  double residual_error_norm = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    double residual_error = fabs(residual.array[iRow]);
    residual_error_norm = std::max(residual_error, residual_error_norm);
  }
  return residual_error_norm;
}

HighsDebugStatus HSimplexNla::debugReportInvertSolutionError(
    const bool transposed, const HVector& true_solution,
    const HVector& solution, HVector& residual, const bool force) const {
  const HighsInt num_row = this->lp_->num_row_;
  double solve_error_norm = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    double solve_error = fabs(solution.array[iRow] - true_solution.array[iRow]);
    solve_error_norm = std::max(solve_error, solve_error_norm);
  }
  double residual_error_norm =
      debugInvertResidualError(transposed, solution, residual);

  return debugReportInvertSolutionError("random solution", transposed,
                                        solve_error_norm, residual_error_norm,
                                        force);
}

HighsDebugStatus HSimplexNla::debugReportInvertSolutionError(
    const std::string& source, const bool transposed,
    const double solve_error_norm, const double residual_error_norm,
    const bool force) const {
  const HighsOptions* options = this->options_;
  std::string value_adjective;
  HighsLogType report_level;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  std::string type = "";
  if (transposed) type = "transposed ";
  if (solve_error_norm) {
    if (solve_error_norm > kSolveExcessiveError) {
      value_adjective = "Excessive";
      report_level = HighsLogType::kError;
    } else if (solve_error_norm > kSolveLargeError) {
      value_adjective = "Large";
      report_level = HighsLogType::kWarning;
    } else {
      value_adjective = "Small";
      report_level = HighsLogType::kInfo;
    }
    if (force) report_level = HighsLogType::kInfo;
    //    printf("%s\n", value_adjective.c_str());
    //    printf("%g\n", solve_error_norm);
    //    printf("%s\n", type.c_str());
    //    printf("%s\n", source.c_str());
    highsLogDev(options->log_options, report_level,
                "CheckINVERT:   %-9s (%9.4g) norm for %s%s solve error\n",
                value_adjective.c_str(), solve_error_norm, type.c_str(),
                source.c_str());
  }

  if (residual_error_norm) {
    if (residual_error_norm > kResidualExcessiveError) {
      value_adjective = "Excessive";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kError;
    } else if (residual_error_norm > kResidualLargeError) {
      value_adjective = "Large";
      report_level = HighsLogType::kWarning;
      return_status = HighsDebugStatus::kWarning;
    } else {
      value_adjective = "Small";
      report_level = HighsLogType::kInfo;
    }
    if (force) report_level = HighsLogType::kInfo;
    highsLogDev(options->log_options, report_level,
                "CheckINVERT:   %-9s (%9.4g) norm for %s%s "
                "residual error\n",
                value_adjective.c_str(), residual_error_norm, type.c_str(),
                source.c_str());
  }
  return return_status;
}

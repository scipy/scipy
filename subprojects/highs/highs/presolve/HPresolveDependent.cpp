/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"

namespace presolve {

HPresolve::Result HPresolve::removeDependentEquations(
    HighsPostsolveStack& postsolve_stack) {
  assert(this->allow_rule_[kPresolveRuleDependentEquations]);
  const bool logging_on = analysis_.logging_on_;
  if (equations.empty()) return Result::kOk;

  auto returnOk = [&]() {
    analysis_.logging_on_ = logging_on;
    if (logging_on)
      analysis_.stopPresolveRuleLog(kPresolveRuleDependentEquations);
    return Result::kOk;
  };

  if (logging_on)
    analysis_.startPresolveRuleLog(kPresolveRuleDependentEquations);

  HighsSparseMatrix matrix;
  HighsInt num_equations = equations.size();
  matrix.num_col_ = num_equations;
  matrix.num_row_ = model->num_col_ + 1;
  matrix.start_.resize(num_equations + 1);
  matrix.start_[0] = 0;
  const HighsInt maxCapacity = numNonzeros() + num_equations;
  matrix.value_.reserve(maxCapacity);
  matrix.index_.reserve(maxCapacity);

  std::vector<HighsInt> eqSet(num_equations);
  std::vector<HighsInt> row_count;
  row_count.assign(model->num_col_, 0);

  HighsInt i = 0;
  for (const std::pair<HighsInt, HighsInt>& p : equations) {
    HighsInt eq = p.second;
    eqSet[i++] = eq;

    // add entries of equation
    for (const HighsSliceNonzero& nonz : getRowVector(eq)) {
      HighsInt iCol = nonz.index();
      row_count[iCol]++;
      matrix.value_.push_back(nonz.value());
      matrix.index_.push_back(iCol);
    }

    // add entry for artificial rhs column
    if (model->row_lower_[eq] != 0.0) {
      matrix.value_.push_back(model->row_lower_[eq]);
      matrix.index_.push_back(model->num_col_);
    }

    matrix.start_[i] = matrix.value_.size();
  }
  // Find the number of (true) variables in the system of equations as
  // the number of columns with entries in at least one equation
  HighsInt num_variables = 0;
  for (HighsInt iCol = 0; iCol < model->num_col_; iCol++)
    if (row_count[iCol]) num_variables++;
  HighsInt num_nz = matrix.numNz();
  const bool silent = silentLog();
  if (!silent)
    highsLogUser(options->log_options, HighsLogType::kInfo,
                 "Considering dependency of %d equation%s in %d variable%s "
                 "with %d nonzero%s\n",
                 int(num_equations), highsIntToPlural(num_equations).c_str(),
                 int(num_variables), highsIntToPlural(num_variables).c_str(),
                 int(num_nz), highsIntToPlural(num_nz).c_str());
  // Identify any dependent equations
  std::vector<HighsInt> colSet(num_equations);
  std::iota(colSet.begin(), colSet.end(), 0);
  HFactor factor;
  factor.setup(matrix, colSet);
  // Set up a time limit to prevent the redundant rows factorization
  // taking forever.
  //
  // Allow no more than 1% of the time limit to be spent on removing
  // dependent equations, but ensure that there is some limit since
  // options->time_limit is infinity by default
  //
  // ToDo: This is strictly non-deterministic, but so conservative
  // that it'll only reap the cases when factor.build never finishes
  const double kMaxDependentEquationsTime = 100;
  const double time_limit = std::max(
      1.0, std::min(0.01 * options->time_limit, kMaxDependentEquationsTime));
  factor.setTimeLimit(time_limit);
  // Determine rank deficiency of the equations
  if (!silent)
    highsLogDev(options->log_options, HighsLogType::kInfo,
                "   Dependent equations search running with time "
                "limit of %.2fs\n",
                time_limit);
  double time_taken = -this->timer->read();
  HighsInt build_return = factor.build();
  time_taken += this->timer->read();
  // Analyse what's been removed
  HighsInt num_removed_row = 0;
  HighsInt num_removed_nz = 0;
  HighsInt num_fictitious_rows_skipped = 0;
  if (build_return == kBuildKernelReturnTimeout) {
    // HFactor::build has timed out, so just return
    if (!silent) {
      if (options->log_dev_level > 0)
        highsLogUser(
            options->log_options, HighsLogType::kInfo,
            "GrepDependentEq,%s,%d,%d,%d,%d,%d,%d,%g,Terminated\n",
            model->model_name_.c_str(), static_cast<int>(num_equations),
            static_cast<int>(num_variables), static_cast<int>(model->num_col_),
            static_cast<int>(num_nz), static_cast<int>(num_removed_row),
            static_cast<int>(num_removed_nz), time_taken);
      highsLogUser(
          options->log_options, HighsLogType::kInfo,
          "   Dependent equations search terminated after %.3gs due to "
          "expected time exceeding limit\n",
          time_taken);
    }
    return returnOk();
  } else {
    double pct_off_timeout =
        1e2 * std::fabs(time_taken - time_limit) / time_limit;
    if (!silent && pct_off_timeout < 1.0)
      highsLogUser(
          options->log_options, HighsLogType::kWarning,
          "   Dependent equations search finished within %.2f%% of limit "
          "of %.2fs: "
          "risk of non-deterministic behaviour if solve is repeated\n",
          pct_off_timeout, time_limit);
  }
  // build_return as rank_deficiency must be valid
  assert(build_return >= 0);
  const HighsInt rank_deficiency = build_return;
  for (HighsInt k = 0; k < rank_deficiency; k++) {
    if (factor.var_with_no_pivot[k] >= 0) {
      HighsInt redundant_row = eqSet[factor.var_with_no_pivot[k]];
      num_removed_row++;
      num_removed_nz += rowsize[redundant_row];
      postsolve_stack.redundantRow(redundant_row);
      removeRow(redundant_row);
    } else {
      num_fictitious_rows_skipped++;
    }
  }
  if (!silent) {
    std::stringstream ss;
    ss.str(std::string());
    ss << highsFormatToString("Dependent equations search");
    if (options->log_dev_level > 0)
      ss << highsFormatToString(
          " with %d / %d variable%s and %d nonzero%s",
          static_cast<int>(num_variables), static_cast<int>(model->num_col_),
          highsIntToPlural(num_variables).c_str(), static_cast<int>(num_nz),
          highsIntToPlural(num_nz).c_str());
    ss << highsFormatToString(" removed %d equation%s and %d nonzero%s",
                              static_cast<int>(num_removed_row),
                              highsIntToPlural(num_removed_row).c_str(),
                              static_cast<int>(num_removed_nz),
                              highsIntToPlural(num_removed_nz).c_str());
    if (options->log_dev_level > 0) {
      ss << highsFormatToString(" in %.2fs with bounds in (%.2f, %.2f)s",
                                time_taken, factor.min_time_bound_,
                                factor.max_time_bound_);
      if (num_fictitious_rows_skipped)
        ss << highsFormatToString(
            ", avoiding %d fictitious row%s",
            static_cast<int>(num_fictitious_rows_skipped),
            highsIntToPlural(num_fictitious_rows_skipped).c_str());
      highsLogDev(options->log_options, HighsLogType::kInfo, "%s\n",
                  ss.str().c_str());
    } else {
      highsLogUser(options->log_options, HighsLogType::kInfo, "%s\n",
                   ss.str().c_str());
    }
  }
  return returnOk();
}

HPresolve::Result HPresolve::removeDependentFreeCols(
    HighsPostsolveStack& postsolve_stack) {
  return Result::kOk;

  // Commented out unreachable code
  //  assert(this->allow_rule_[kPresolveRuleDependentFreeCols]);
  //  const bool logging_on = analysis_.logging_on_;
  //  if (logging_on)
  //    analysis_.startPresolveRuleLog(kPresolveRuleDependentFreeCols);
  //
  //  // todo the postsolve step does not work properly
  //  std::vector<HighsInt> freeCols;
  //  freeCols.reserve(model->num_col_);
  //
  //  for (HighsInt i = 0; i < model->num_col_; ++i) {
  //    if (colDeleted[i]) continue;
  //    if (model->col_lower_[i] == -kHighsInf && model->col_upper_[i] ==
  //    kHighsInf)
  //      freeCols.push_back(i);
  //  }
  //
  //  if (freeCols.empty()) return Result::kOk;
  //
  //  HighsSparseMatrix matrix;
  //  matrix.num_col_ = freeCols.size();
  //  highsLogDev(options->log_options, HighsLogType::kInfo,
  //              "HPresolve::removeDependentFreeCols Got %d free cols, checking
  //              " "for dependent free cols\n", (int)matrix.num_col_);
  //  matrix.num_row_ = model->num_row_ + 1;
  //  matrix.start_.resize(matrix.num_col_ + 1);
  //  matrix.start_[0] = 0;
  //  const HighsInt maxCapacity = numNonzeros() + matrix.num_col_;
  //  matrix.value_.reserve(maxCapacity);
  //  matrix.index_.reserve(maxCapacity);
  //
  //  for (HighsInt i = 0; i < matrix.num_col_; ++i) {
  //    HighsInt col = freeCols[i];
  //    // add entries of free column
  //    for (const HighsSliceNonzero& nonz : getColumnVector(col)) {
  //      matrix.value_.push_back(nonz.value());
  //      matrix.index_.push_back(nonz.index());
  //    }
  //
  //    // add entry for artificial cost row
  //    if (model->col_cost_[col] != 0.0) {
  //      matrix.value_.push_back(model->col_cost_[col]);
  //      matrix.index_.push_back(model->num_row_);
  //    }
  //
  //    matrix.start_[i + 1] = matrix.value_.size();
  //  }
  //  std::vector<HighsInt> colSet(matrix.num_col_);
  //  std::iota(colSet.begin(), colSet.end(), 0);
  //  HFactor factor;
  //  factor.setup(matrix, colSet);
  //  HighsInt rank_deficiency = factor.build();
  //  // Must not have timed out
  //  assert(rank_deficiency >= 0);
  //  highsLogDev(options->log_options, HighsLogType::kInfo,
  //              "HPresolve::removeDependentFreeCols Got %d free cols, checking
  //              " "for dependent free cols\n", (int)matrix.num_col_);
  //  // Analyse what's been removed
  //  HighsInt num_removed_row = 0;
  //  HighsInt num_removed_nz = 0;
  //  HighsInt num_fictitious_cols_skipped = 0;
  //  for (HighsInt k = 0; k < rank_deficiency; k++) {
  //    if (factor.var_with_no_pivot[k] >= 0) {
  //      HighsInt redundant_col = freeCols[factor.var_with_no_pivot[k]];
  //      num_removed_nz += colsize[redundant_col];
  //      fixColToZero(postsolve_stack, redundant_col);
  //    } else {
  //      num_fictitious_cols_skipped++;
  //    }
  //  }
  //  highsLogDev(
  //      options->log_options, HighsLogType::kInfo,
  //      "HPresolve::removeDependentFreeCols Removed %d rows and %d nonzeros",
  //      (int)num_removed_row, (int)num_removed_nz);
  //  if (num_fictitious_cols_skipped)
  //    highsLogDev(options->log_options, HighsLogType::kInfo,
  //                ", avoiding %d fictitious rows",
  //                (int)num_fictitious_cols_skipped);
  //  highsLogDev(options->log_options, HighsLogType::kInfo, "\n");
  //
  //  analysis_.logging_on_ = logging_on;
  //  if (logging_on)
  //  analysis_.stopPresolveRuleLog(kPresolveRuleDependentFreeCols);
  //
  //  return Result::kOk;
}

}  // namespace presolve

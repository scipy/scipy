/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 */
#include "lp_data/HighsLpUtils.h"

#include <algorithm>
#include <cassert>

#include "HConfig.h"
#include "io/Filereader.h"
#include "io/HMPSIO.h"
#include "io/HighsIO.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsCDouble.h"
#include "util/HighsMatrixUtils.h"
#include "util/HighsSort.h"
#include "util/stringutil.h"

using std::fabs;
using std::max;
using std::min;

const HighsInt kMaxLineLength = 80;

HighsStatus assessLp(HighsLp& lp, const HighsOptions& options) {
  const HighsLogOptions& log_options = options.log_options;
  if (!lp.legalDimensions("assessLp", log_options)) return HighsStatus::kError;
  HighsStatus return_status = lp.trimDimensions(log_options);

  HighsStatus call_status = HighsStatus::kOk;
  if (lp.num_col_) {
    // Assess the LP column costs
    HighsIndexCollection index_collection;
    index_collection.dimension_ = lp.num_col_;
    index_collection.is_interval_ = true;
    index_collection.from_ = 0;
    index_collection.to_ = lp.num_col_ - 1;
    call_status = assessCosts(options, 0, index_collection, lp.col_cost_,
                              lp.has_infinite_cost_, options.infinite_cost);
    return_status = interpretCallStatus(log_options, call_status, return_status,
                                        "assessCosts");
    if (return_status == HighsStatus::kError) return return_status;
    // Assess the LP column bounds
    call_status = assessBounds(
        options, "Col", 0, index_collection, lp.col_lower_, lp.col_upper_,
        options.infinite_bound, lp.isMip() ? lp.integrality_.data() : nullptr,
        lp.col_names_.data());
    return_status = interpretCallStatus(log_options, call_status, return_status,
                                        "assessBounds");
    if (return_status == HighsStatus::kError) return return_status;
  }
  if (lp.num_row_) {
    // Assess the LP row bounds
    HighsIndexCollection index_collection;
    index_collection.dimension_ = lp.num_row_;
    index_collection.is_interval_ = true;
    index_collection.from_ = 0;
    index_collection.to_ = lp.num_row_ - 1;
    call_status = assessBounds(
        options, "Row", 0, index_collection, lp.row_lower_, lp.row_upper_,
        options.infinite_bound, nullptr, lp.row_names_.data());
    return_status = interpretCallStatus(log_options, call_status, return_status,
                                        "assessBounds");
    if (return_status == HighsStatus::kError) return return_status;
  }
  // If the LP has no columns the matrix must be empty and there is
  // nothing left to test
  if (lp.num_col_ == 0) {
    assert(!lp.numNz());
    return return_status;
  }
  // From here, any LP has lp.num_col_ > 0 and lp.a_matrix_.start_[lp.num_col_]
  // exists (as the number of nonzeros)
  assert(lp.num_col_ > 0);

  // Assess the LP matrix - even if there are no rows!
  const bool sum_duplicates = false;
  call_status = lp.a_matrix_.assess(log_options, options.small_matrix_value,
                                    options.large_matrix_value, sum_duplicates,
                                    lp.col_names_.data(), lp.row_names_.data());
  return_status = interpretCallStatus(log_options, call_status, return_status,
                                      "assessMatrix");
  if (return_status == HighsStatus::kError) return return_status;

  // Matrix should have been trimmed in lp.a_matrix_.assess
  assert(lp.a_matrix_.trimDimensions(log_options) == HighsStatus::kOk);
  if (return_status != HighsStatus::kOk)
    highsLogDev(log_options, HighsLogType::kInfo,
                "assessLp returns HighsStatus = %s\n",
                highsStatusToString(return_status).c_str());
  return return_status;
}

HighsStatus assessCosts(const HighsOptions& options, const HighsInt ml_col_os,
                        const HighsIndexCollection& index_collection,
                        vector<double>& cost, bool& has_infinite_cost,
                        const double infinite_cost) {
  HighsStatus return_status = HighsStatus::kOk;
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return return_status;

  return_status = HighsStatus::kOk;
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

  HighsInt local_col;
  HighsInt usr_col = -1;
  HighsInt num_infinite_cost = 0;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (index_collection.is_interval_ || index_collection.is_mask_) {
      local_col = k;
    } else {
      local_col = index_collection.set_[k];
    }
    if (index_collection.is_interval_) {
      usr_col++;
    } else {
      usr_col = k;
    }
    if (index_collection.is_mask_ && !index_collection.mask_[local_col])
      continue;
    if (cost[usr_col] >= infinite_cost) {
      num_infinite_cost++;
      cost[usr_col] = kHighsInf;
    } else if (cost[usr_col] <= -infinite_cost) {
      num_infinite_cost++;
      cost[usr_col] = -kHighsInf;
    }
  }
  if (num_infinite_cost > 0) {
    has_infinite_cost = true;
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "%" HIGHSINT_FORMAT
                 " |cost| value%s greater than or equal to %12g %s treated as "
                 "Infinity\n",
                 num_infinite_cost, highsIntToPlural(num_infinite_cost).c_str(),
                 infinite_cost, highsIntToBe(num_infinite_cost).c_str());
  }
  assert(return_status == HighsStatus::kOk);
  return return_status;
}

HighsStatus assessBounds(const HighsOptions& options, const std::string& type,
                         const HighsInt ml_ix_os,
                         const HighsIndexCollection& index_collection,
                         vector<double>& lower, vector<double>& upper,
                         const double infinite_bound,
                         const HighsVarType* integrality,
                         const std::string* names) {
  HighsStatus return_status = HighsStatus::kOk;
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return HighsStatus::kOk;

  return_status = HighsStatus::kOk;
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
  auto possible_name = [&](const HighsInt ix) {
    std::string name = "";
    if (names) name = " (" + type + " " + names[ix] + ")";
    return name;
  };

  HighsInt num_infinite_lower_bound = 0;
  HighsInt num_infinite_upper_bound = 0;
  HighsInt local_ix;
  HighsInt ml_ix;
  HighsInt usr_ix = -1;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (index_collection.is_interval_ || index_collection.is_mask_) {
      local_ix = k;
    } else {
      local_ix = index_collection.set_[k];
    }
    if (index_collection.is_interval_) {
      usr_ix++;
    } else {
      usr_ix = k;
    }
    ml_ix = ml_ix_os + local_ix;
    if (index_collection.is_mask_ && !index_collection.mask_[local_ix])
      continue;

    if (!highs_isInfinity(-lower[usr_ix])) {
      // Check whether a finite lower bound will be treated as -Infinity
      bool infinite_lower_bound = lower[usr_ix] <= -infinite_bound;
      if (infinite_lower_bound) {
        lower[usr_ix] = -kHighsInf;
        num_infinite_lower_bound++;
      }
    }
    if (!highs_isInfinity(upper[usr_ix])) {
      // Check whether a finite upper bound will be treated as Infinity
      bool infinite_upper_bound = upper[usr_ix] >= infinite_bound;
      if (infinite_upper_bound) {
        upper[usr_ix] = kHighsInf;
        num_infinite_upper_bound++;
      }
    }
    // Check that the lower bound does not exceed the upper bound
    bool legalLowerUpperBound = lower[usr_ix] <= upper[usr_ix];
    if (integrality) {
      // Legal for semi-variables to have inconsistent bounds
      if (integrality[usr_ix] == HighsVarType::kSemiContinuous ||
          integrality[usr_ix] == HighsVarType::kSemiInteger)
        legalLowerUpperBound = true;
    }
    if (!legalLowerUpperBound) {
      // Leave inconsistent bounds to be used to deduce
      // infeasibility. Cannot lissue a warning, as they are not
      // removed (since the infeasibility is trivially discovered when
      // there is an attempt to solve the model) and the call to
      // Highs::model_.assess in Highs::run cannot yield a warning.
      highsLogUser(options.log_options, HighsLogType::kInfo,
                   "%3s  %12" HIGHSINT_FORMAT
                   "%s has inconsistent bounds [%12g, %12g]\n",
                   type.c_str(), ml_ix, possible_name(ml_ix).c_str(),
                   lower[usr_ix], upper[usr_ix]);
    }
    // Check that the lower bound is not as much as +Infinity
    bool legalLowerBound = lower[usr_ix] < infinite_bound;
    if (!legalLowerBound) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "%3s  %12" HIGHSINT_FORMAT
                   "%s  has lower bound of %12g >= %12g\n",
                   type.c_str(), ml_ix, possible_name(ml_ix).c_str(),
                   lower[usr_ix], infinite_bound);
      error_found = true;
    }
    // Check that the upper bound is not as little as -Infinity
    bool legalUpperBound = upper[usr_ix] > -infinite_bound;
    if (!legalUpperBound) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "%3s  %12" HIGHSINT_FORMAT
                   "%s  has upper bound of %12g <= %12g\n",
                   type.c_str(), ml_ix, possible_name(ml_ix).c_str(),
                   upper[usr_ix], -infinite_bound);
      error_found = true;
    }
  }
  if (num_infinite_lower_bound) {
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "%3ss:%12" HIGHSINT_FORMAT
                 " lower bound%s    less than or equal to %12g are treated as "
                 "-Infinity\n",
                 type.c_str(), num_infinite_lower_bound,
                 highsIntToPlural(num_infinite_lower_bound).c_str(),
                 -infinite_bound);
  }
  if (num_infinite_upper_bound) {
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "%3ss:%12" HIGHSINT_FORMAT
                 " upper bound%s greater than or equal to %12g are treated as "
                 "+Infinity\n",
                 type.c_str(), num_infinite_upper_bound,
                 highsIntToPlural(num_infinite_upper_bound).c_str(),
                 infinite_bound);
  }

  assert(!warning_found);
  if (error_found)
    return_status = HighsStatus::kError;
  else if (warning_found)
    return_status = HighsStatus::kWarning;
  else
    return_status = HighsStatus::kOk;

  return return_status;
}

HighsStatus assessSemiVariables(HighsLp& lp, const HighsOptions& options,
                                bool& made_semi_variable_mods) {
  made_semi_variable_mods = false;
  HighsStatus return_status = HighsStatus::kOk;
  if (!lp.integrality_.size()) return return_status;
  assert((HighsInt)lp.integrality_.size() == lp.num_col_);
  HighsInt num_illegal_lower = 0;
  HighsInt num_illegal_upper = 0;
  HighsInt num_tightened_upper = 0;
  HighsInt num_inconsistent_semi = 0;
  HighsInt num_non_semi = 0;
  HighsInt num_non_continuous_variables = 0;
  const double kLowerBoundMu = 10.0;

  std::vector<HighsInt>& save_non_semi_variable_index =
      lp.mods_.save_non_semi_variable_index;
  std::vector<HighsInt>& inconsistent_semi_variable_index =
      lp.mods_.save_inconsistent_semi_variable_index;
  std::vector<double>& inconsistent_semi_variable_lower_bound_value =
      lp.mods_.save_inconsistent_semi_variable_lower_bound_value;
  std::vector<double>& inconsistent_semi_variable_upper_bound_value =
      lp.mods_.save_inconsistent_semi_variable_upper_bound_value;
  std::vector<HighsVarType>& inconsistent_semi_variable_type =
      lp.mods_.save_inconsistent_semi_variable_type;

  std::vector<HighsInt>& tightened_semi_variable_upper_bound_index =
      lp.mods_.save_tightened_semi_variable_upper_bound_index;
  std::vector<double>& tightened_semi_variable_upper_bound_value =
      lp.mods_.save_tightened_semi_variable_upper_bound_value;

  assert(int(lp.mods_.save_inconsistent_semi_variable_index.size()) == 0);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
        lp.integrality_[iCol] == HighsVarType::kSemiInteger) {
      if (lp.col_lower_[iCol] > lp.col_upper_[iCol]) {
        // Semi-variables with inconsistent bounds become continuous
        // and fixed at zero
        num_inconsistent_semi++;
        inconsistent_semi_variable_index.push_back(iCol);
        inconsistent_semi_variable_lower_bound_value.push_back(
            lp.col_lower_[iCol]);
        inconsistent_semi_variable_upper_bound_value.push_back(
            lp.col_upper_[iCol]);
        inconsistent_semi_variable_type.push_back(lp.integrality_[iCol]);
        continue;
      }
      // Semi-variables with zero lower bound are not semi
      if (lp.col_lower_[iCol] == 0) {
        num_non_semi++;
        save_non_semi_variable_index.push_back(iCol);
        // Semi-integer become integer so still have a non-continuous variable
        if (lp.integrality_[iCol] == HighsVarType::kSemiInteger)
          num_non_continuous_variables++;
        continue;
      }
      if (lp.col_lower_[iCol] < 0) {
        // Semi-variables must have a positive lower bound
        num_illegal_lower++;
      } else if (lp.col_upper_[iCol] > kMaxSemiVariableUpper) {
        // Semi-variables must have upper bound that's not too large,
        // so see whether the limiting value is sufficiently larger than the
        // lower bound
        if (kLowerBoundMu * lp.col_lower_[iCol] > kMaxSemiVariableUpper) {
          num_illegal_upper++;
        } else {
          // Record the upper bound change to be made later
          tightened_semi_variable_upper_bound_index.push_back(iCol);
          tightened_semi_variable_upper_bound_value.push_back(
              kMaxSemiVariableUpper);
          num_tightened_upper++;
        }
      }
      num_non_continuous_variables++;
    } else if (lp.integrality_[iCol] == HighsVarType::kInteger) {
      num_non_continuous_variables++;
    }
  }
  if (num_inconsistent_semi) {
    highsLogUser(
        options.log_options, HighsLogType::kWarning,
        "%" HIGHSINT_FORMAT
        " semi-continuous/integer variable(s) have inconsistent bounds "
        "so are fixed at zero\n",
        num_inconsistent_semi);
    return_status = HighsStatus::kWarning;
  }
  if (num_non_semi) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "%" HIGHSINT_FORMAT
                 " semi-continuous/integer variable(s) have zero lower bound "
                 "so are continuous/integer\n",
                 num_non_semi);
    return_status = HighsStatus::kWarning;
  }
  if (!num_non_continuous_variables) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "No semi-integer/integer variables in model with non-empty "
                 "integrality\n");
    return_status = HighsStatus::kWarning;
  }
  const bool has_illegal_bounds = num_illegal_lower || num_illegal_upper;
  if (num_tightened_upper) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "%" HIGHSINT_FORMAT
                 " semi-continuous/integer variable(s) have upper bounds "
                 "exceeding %g that can be tightened to %g > %g*lower)\n",
                 num_tightened_upper, kMaxSemiVariableUpper,
                 kMaxSemiVariableUpper, kLowerBoundMu);
    return_status = HighsStatus::kWarning;
    if (has_illegal_bounds) {
      // Don't apply upper bound tightenings if there are illegal bounds
      assert(num_illegal_lower || num_illegal_upper);
      tightened_semi_variable_upper_bound_index.clear();
      tightened_semi_variable_upper_bound_value.clear();
      num_tightened_upper = 0;
    } else {
      // Apply the upper bound tightenings, saving the over-written
      // values
      for (HighsInt k = 0; k < num_tightened_upper; k++) {
        const double use_upper_bound =
            tightened_semi_variable_upper_bound_value[k];
        const HighsInt iCol = tightened_semi_variable_upper_bound_index[k];
        tightened_semi_variable_upper_bound_value[k] = lp.col_upper_[iCol];
        lp.col_upper_[iCol] = use_upper_bound;
      }
    }
  }
  if (num_inconsistent_semi) {
    if (has_illegal_bounds) {
      // Don't apply bound changes if there are illegal bounds
      inconsistent_semi_variable_index.clear();
      inconsistent_semi_variable_lower_bound_value.clear();
      inconsistent_semi_variable_upper_bound_value.clear();
      inconsistent_semi_variable_type.clear();
      num_inconsistent_semi = 0;
    } else {
      for (HighsInt k = 0; k < num_inconsistent_semi; k++) {
        const HighsInt iCol = inconsistent_semi_variable_index[k];
        lp.col_lower_[iCol] = 0;
        lp.col_upper_[iCol] = 0;
        lp.integrality_[iCol] = HighsVarType::kContinuous;
      }
    }
  }
  if (num_non_semi) {
    if (has_illegal_bounds) {
      // Don't apply type changes if there are illegal bounds
      save_non_semi_variable_index.clear();
    } else {
      for (HighsInt k = 0; k < num_non_semi; k++) {
        const HighsInt iCol = save_non_semi_variable_index[k];
        if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous) {
          // Semi-continuous become continuous
          lp.integrality_[iCol] = HighsVarType::kContinuous;
        } else {
          // Semi-integer become integer
          lp.integrality_[iCol] = HighsVarType::kInteger;
        }
      }
    }
  }
  if (num_illegal_lower) {
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "%" HIGHSINT_FORMAT
        " semi-continuous/integer variable(s) have negative lower bounds\n",
        num_illegal_lower);
    return_status = HighsStatus::kError;
  }
  if (num_illegal_upper) {
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "%" HIGHSINT_FORMAT
        " semi-continuous/integer variables have upper bounds "
        "exceeding %g that cannot be modified due to large lower bounds\n",
        num_illegal_upper, kMaxSemiVariableUpper);
    return_status = HighsStatus::kError;
  }
  made_semi_variable_mods =
      num_non_semi > 0 || num_inconsistent_semi > 0 || num_tightened_upper > 0;
  assert(static_cast<size_t>(num_non_semi) <=
         save_non_semi_variable_index.size());
  assert(static_cast<size_t>(num_inconsistent_semi) <=
         inconsistent_semi_variable_index.size());
  assert(static_cast<size_t>(num_tightened_upper) <=
         tightened_semi_variable_upper_bound_index.size());
  //      save_non_semi_variable_index.size() > 0 ||
  //      inconsistent_semi_variable_index.size() > 0 ||
  //      tightened_semi_variable_upper_bound_index.size() > 0;
  return return_status;
}

void relaxSemiVariables(HighsLp& lp, bool& made_semi_variable_mods) {
  // When solving relaxation, semi-variables are continuous between 0
  // and their upper bound, so have to modify the lower bound to be
  // zero
  made_semi_variable_mods = false;
  if (!lp.integrality_.size()) return;
  assert((HighsInt)lp.integrality_.size() == lp.num_col_);
  std::vector<HighsInt>& relaxed_semi_variable_lower_index =
      lp.mods_.save_relaxed_semi_variable_lower_bound_index;
  std::vector<double>& relaxed_semi_variable_lower_value =
      lp.mods_.save_relaxed_semi_variable_lower_bound_value;
  assert(relaxed_semi_variable_lower_index.size() == 0);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
        lp.integrality_[iCol] == HighsVarType::kSemiInteger) {
      relaxed_semi_variable_lower_index.push_back(iCol);
      relaxed_semi_variable_lower_value.push_back(lp.col_lower_[iCol]);
      lp.col_lower_[iCol] = 0;
    }
  }
  made_semi_variable_mods = relaxed_semi_variable_lower_index.size() > 0;
}

bool activeModifiedUpperBounds(const HighsOptions& options, const HighsLp& lp,
                               const std::vector<double>& col_value) {
  const std::vector<HighsInt>& tightened_semi_variable_upper_bound_index =
      lp.mods_.save_tightened_semi_variable_upper_bound_index;
  const HighsInt num_tightened_upper =
      tightened_semi_variable_upper_bound_index.size();
  HighsInt num_active_modified_upper = 0;
  double min_semi_variable_margin = kHighsInf;
  for (HighsInt k = 0; k < num_tightened_upper; k++) {
    const double value =
        col_value[tightened_semi_variable_upper_bound_index[k]];
    const double upper =
        lp.col_upper_[tightened_semi_variable_upper_bound_index[k]];
    double semi_variable_margin = upper - value;
    if (value > upper - options.primal_feasibility_tolerance) {
      min_semi_variable_margin = 0;
      num_active_modified_upper++;
    } else {
      min_semi_variable_margin =
          std::min(semi_variable_margin, min_semi_variable_margin);
    }
  }
  if (num_active_modified_upper) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "%" HIGHSINT_FORMAT
                 " semi-variables are active at modified upper bounds\n",
                 num_active_modified_upper);
  } else if (num_tightened_upper) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "No semi-variables are active at modified upper bounds:"
                 " a large minimum margin (%g) suggests optimality,"
                 " but there is no guarantee\n",
                 min_semi_variable_margin);
  }
  return (num_active_modified_upper != 0);
}

HighsStatus cleanBounds(const HighsOptions& options, HighsLp& lp) {
  double max_residual = 0;
  HighsInt num_change = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    double residual = lp.col_lower_[iCol] - lp.col_upper_[iCol];
    if (residual > options.primal_feasibility_tolerance) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Column %" HIGHSINT_FORMAT
                   " has inconsistent bounds [%g, %g] (residual = "
                   "%g) after presolve\n",
                   iCol, lp.col_lower_[iCol], lp.col_upper_[iCol], residual);
      return HighsStatus::kError;
    } else if (residual > 0) {
      num_change++;
      max_residual = std::max(residual, max_residual);
      double mid = 0.5 * (lp.col_lower_[iCol] + lp.col_upper_[iCol]);
      lp.col_lower_[iCol] = mid;
      lp.col_upper_[iCol] = mid;
    }
  }
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    double residual = lp.row_lower_[iRow] - lp.row_upper_[iRow];
    if (residual > options.primal_feasibility_tolerance) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Row %" HIGHSINT_FORMAT
                   " has inconsistent bounds [%g, %g] (residual = %g) "
                   "after presolve\n",
                   iRow, lp.row_lower_[iRow], lp.row_upper_[iRow], residual);
      return HighsStatus::kError;
    } else if (residual > 0) {
      num_change++;
      max_residual = std::max(residual, max_residual);
      double mid = 0.5 * (lp.row_lower_[iRow] + lp.row_upper_[iRow]);
      lp.row_lower_[iRow] = mid;
      lp.row_upper_[iRow] = mid;
    }
  }
  if (num_change) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Resolved %" HIGHSINT_FORMAT
                 " inconsistent bounds (maximum residual = "
                 "%9.4g) after presolve\n",
                 num_change, max_residual);
    return HighsStatus::kWarning;
  }
  return HighsStatus::kOk;
}

HighsStatus userScaleLp(HighsLp& lp, HighsUserScaleData& data,
                        const HighsLogOptions& log_options) {
  userScaleLp(lp, data, false);
  HighsStatus return_status = userScaleStatus(log_options, data);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  userScaleLp(lp, data);
  return return_status;
}

void userScaleLp(HighsLp& lp, HighsUserScaleData& data, const bool apply) {
  userScaleCosts(lp.integrality_, lp.offset_, lp.col_cost_, data, apply);
  userScaleColBounds(lp.integrality_, lp.col_lower_, lp.col_upper_, data,
                     apply);
  userScaleMatrix(lp.integrality_, lp.a_matrix_, data, apply);
  userScaleRowBounds(lp.row_lower_, lp.row_upper_, data, apply);
}

void userScaleCosts(const vector<HighsVarType>& integrality, double& offset,
                    vector<double>& cost, HighsUserScaleData& data,
                    const bool apply) {
  data.num_infinite_costs = 0;
  const HighsInt user_bound_scale = data.user_bound_scale;
  const HighsInt user_objective_scale = data.user_objective_scale;
  if (!user_bound_scale && !user_objective_scale) return;
  const HighsInt num_col = cost.size();
  if (num_col <= 0) return;
  const HighsInt integrality_size = HighsInt(integrality.size());
  const bool has_integrality = integrality_size > 0;
  double bound_scale_value = std::pow(2, user_bound_scale);
  double objective_scale_value = std::pow(2, user_objective_scale);
  assert(!has_integrality || integrality_size >= num_col);
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    double value = cost[iCol];
    if (has_integrality && integrality[iCol] != HighsVarType::kContinuous)
      value *= bound_scale_value;
    value *= objective_scale_value;
    if (std::abs(value) > data.infinite_cost) data.num_infinite_costs++;
    if (apply) cost[iCol] = value;
  }
  if (apply) offset *= (bound_scale_value * objective_scale_value);
}

void userScaleColBounds(const vector<HighsVarType>& integrality,
                        vector<double>& lower, vector<double>& upper,
                        HighsUserScaleData& data, const bool apply) {
  data.num_infinite_col_bounds = 0;
  const HighsInt user_bound_scale = data.user_bound_scale;
  if (!user_bound_scale) return;
  const HighsInt num_col = lower.size();
  if (num_col <= 0) return;
  const HighsInt integrality_size = HighsInt(integrality.size());
  const bool has_integrality = integrality_size > 0;
  assert(!has_integrality || integrality_size >= num_col);
  double bound_scale_value = std::pow(2, user_bound_scale);
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (!has_integrality || integrality[iCol] == HighsVarType::kContinuous) {
      if (lower[iCol] > -kHighsInf) {
        double value = lower[iCol] * bound_scale_value;
        if (std::abs(value) > data.infinite_bound)
          data.num_infinite_col_bounds++;
        if (apply) lower[iCol] = value;
      }
      if (upper[iCol] < kHighsInf) {
        double value = upper[iCol] * bound_scale_value;
        if (std::abs(value) > data.infinite_bound)
          data.num_infinite_col_bounds++;
        if (apply) upper[iCol] = value;
      }
    }
  }
}

void userScaleRowBounds(vector<double>& lower, vector<double>& upper,
                        HighsUserScaleData& data, const bool apply) {
  data.num_infinite_row_bounds = 0;
  const HighsInt user_bound_scale = data.user_bound_scale;
  if (!user_bound_scale) return;
  const HighsInt num_row = lower.size();
  if (num_row <= 0) return;
  double bound_scale_value = std::pow(2, user_bound_scale);
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (lower[iRow] > -kHighsInf) {
      double value = lower[iRow] * bound_scale_value;
      if (std::abs(value) > data.infinite_bound) data.num_infinite_row_bounds++;
      if (apply) lower[iRow] = value;
    }
    if (upper[iRow] < kHighsInf) {
      double value = upper[iRow] * bound_scale_value;
      if (std::abs(value) > data.infinite_bound) data.num_infinite_row_bounds++;
      if (apply) upper[iRow] = value;
    }
  }
}

void userScaleMatrix(const vector<HighsVarType>& integrality,
                     HighsSparseMatrix& matrix, HighsUserScaleData& data,
                     const bool apply) {
  data.num_small_matrix_values = 0;
  data.num_large_matrix_values = 0;
  const HighsInt user_bound_scale = data.user_bound_scale;
  if (!user_bound_scale) return;
  if (!integrality.size()) return;
  const HighsInt num_col = matrix.num_col_;
  if (num_col <= 0) return;
  const HighsInt num_row = matrix.num_row_;
  if (num_row <= 0) return;
  assert(HighsInt(integrality.size()) >= num_col);
  double bound_scale_value = std::pow(2, user_bound_scale);
  if (matrix.isColwise()) {
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      if (integrality[iCol] == HighsVarType::kContinuous) continue;
      for (HighsInt iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol + 1];
           iEl++) {
        double value = matrix.value_[iEl] * bound_scale_value;
        double abs_value = std::fabs(value);
        if (abs_value <= data.small_matrix_value)
          data.num_small_matrix_values++;
        else if (abs_value >= data.large_matrix_value)
          data.num_large_matrix_values++;
        if (apply) matrix.value_[iEl] = value;
      }
    }
  } else {
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      for (HighsInt iEl = matrix.start_[iRow]; iEl < matrix.start_[iRow + 1];
           iEl++) {
        HighsInt iCol = matrix.index_[iEl];
        if (integrality[iCol] == HighsVarType::kContinuous) continue;
        double value = matrix.value_[iEl] * bound_scale_value;
        double abs_value = std::fabs(value);
        if (abs_value <= data.small_matrix_value)
          data.num_small_matrix_values++;
        else if (abs_value >= data.large_matrix_value)
          data.num_large_matrix_values++;
        if (apply) matrix.value_[iEl] = value;
      }
    }
  }
}

HighsStatus userScaleStatus(const HighsLogOptions& log_options,
                            const HighsUserScaleData& data) {
  HighsStatus return_status = HighsStatus::kOk;
  std::string message;
  if (data.scaleWarning(message)) {
    highsLogUser(log_options, HighsLogType::kWarning, "%s\n", message.c_str());
    return_status = HighsStatus::kWarning;
  }
  if (data.scaleError(message)) {
    highsLogUser(log_options, HighsLogType::kError, "%s\n", message.c_str());
    return_status = HighsStatus::kError;
  }
  return return_status;
}

bool considerScaling(const HighsOptions& options, HighsLp& lp) {
  // Indicate whether new scaling has been determined in the return value.
  bool new_scaling = false;
  // Consider scaling the LP - either by finding new factors or by
  // applying any existing factors
  const bool allow_scaling =
      lp.num_col_ > 0 &&
      options.simplex_scale_strategy != kSimplexScaleStrategyOff;
  if (lp.scale_.has_scaling && !allow_scaling) {
    // LP had scaling before, but now it is not permitted, clear any
    // scaling. Return true as scaling position has changed
    lp.clearScale();
    return true;
  }
  const bool scaling_not_tried = lp.scale_.strategy == kSimplexScaleStrategyOff;
  const bool new_scaling_strategy =
      options.simplex_scale_strategy != lp.scale_.strategy &&
      options.simplex_scale_strategy != kSimplexScaleStrategyChoose;
  const bool try_scaling =
      allow_scaling && (scaling_not_tried || new_scaling_strategy);
  if (try_scaling) {
    // Scaling will be tried, so ensure that any previous scaling is not applied
    lp.unapplyScale();
    const bool analyse_lp_data =
        kHighsAnalysisLevelModelData & options.highs_analysis_level;
    if (analyse_lp_data) analyseLp(options.log_options, lp);
    scaleLp(options, lp);
    // If the LP is now scaled, then the scaling is new
    new_scaling = lp.is_scaled_;
    if (analyse_lp_data && lp.is_scaled_) analyseLp(options.log_options, lp);
  } else if (lp.scale_.has_scaling) {
    // Scaling factors are known, so ensure that they are applied
    lp.applyScale();
  }
  // Ensure that either the LP has scale factors and is scaled, or
  // it doesn't have scale factors and isn't scaled
  assert(lp.scale_.has_scaling == lp.is_scaled_);
  return new_scaling;
}

void scaleLp(const HighsOptions& options, HighsLp& lp,
             const bool force_scaling) {
  lp.clearScaling();
  HighsInt numCol = lp.num_col_;
  HighsInt numRow = lp.num_row_;
  // Scaling not well defined for models with no columns
  assert(numCol > 0);
  vector<double>& colCost = lp.col_cost_;
  vector<double>& colLower = lp.col_lower_;
  vector<double>& colUpper = lp.col_upper_;
  vector<double>& rowLower = lp.row_lower_;
  vector<double>& rowUpper = lp.row_upper_;

  // Save the simplex_scale_strategy so that the option can be
  // modified for the course of this method
  HighsInt simplex_scale_strategy = options.simplex_scale_strategy;
  // Determine the actual strategy to use
  HighsInt use_scale_strategy = simplex_scale_strategy;
  if (use_scale_strategy == kSimplexScaleStrategyChoose) {
    // HiGHS is left to choose: currently use forced equilibration, but maybe do
    // something more intelligent
    use_scale_strategy = kSimplexScaleStrategyForcedEquilibration;
  }
  // Find out range of matrix values and skip matrix scaling if all
  // |values| are in [0.2, 5]
  const double no_scaling_original_matrix_min_value = 0.2;
  const double no_scaling_original_matrix_max_value = 5.0;
  double original_matrix_min_value = kHighsInf;
  double original_matrix_max_value = 0;
  lp.a_matrix_.range(original_matrix_min_value, original_matrix_max_value);
  // Possibly force scaling, otherwise base the decision on the range
  // of values in the matrix, values that will be used later for
  // reporting
  const bool no_scaling = force_scaling
                              ? false
                              : (original_matrix_min_value >=
                                 no_scaling_original_matrix_min_value) &&
                                    (original_matrix_max_value <=
                                     no_scaling_original_matrix_max_value);
  bool scaled_matrix = false;
  if (no_scaling) {
    // No matrix scaling, but possible cost scaling
    if (options.highs_analysis_level)
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Matrix has [min, max] values of [%g, %g] within "
                  "[%g, %g] so no scaling performed\n",
                  original_matrix_min_value, original_matrix_max_value,
                  no_scaling_original_matrix_min_value,
                  no_scaling_original_matrix_max_value);
  } else {
    // Try scaling, so assign unit factors - partly because initial
    // factors may be assumed by the scaling method, but also because
    // scaling factors may not be computed for empty rows/columns
    HighsScale& scale = lp.scale_;
    scale.col.assign(numCol, 1);
    scale.row.assign(numRow, 1);
    const bool equilibration_scaling =
        use_scale_strategy == kSimplexScaleStrategyEquilibration ||
        use_scale_strategy == kSimplexScaleStrategyForcedEquilibration;
    // Try scaling. Value of scaled_matrix indicates whether scaling
    // was considered valuable (and performed). If it's not valuable
    // then the matrix remains unscaled
    if (equilibration_scaling) {
      scaled_matrix = equilibrationScaleMatrix(options, lp, use_scale_strategy);
    } else {
      scaled_matrix = maxValueScaleMatrix(options, lp, use_scale_strategy);
    }
    if (scaled_matrix) {
      // Matrix is scaled, so scale the bounds and costs
      for (HighsInt iCol = 0; iCol < numCol; iCol++) {
        colLower[iCol] /= scale.col[iCol];
        colUpper[iCol] /= scale.col[iCol];
        colCost[iCol] *= scale.col[iCol];
      }
      for (HighsInt iRow = 0; iRow < numRow; iRow++) {
        rowLower[iRow] *= scale.row[iRow];
        rowUpper[iRow] *= scale.row[iRow];
      }
      scale.has_scaling = true;
      scale.num_col = numCol;
      scale.num_row = numRow;
      scale.cost = 1.0;
      lp.is_scaled_ = true;
    } else {
      // Matrix is not scaled, so clear the scaling
      lp.clearScaling();
    }
  }
  // Record the scaling strategy used
  lp.scale_.strategy = use_scale_strategy;
  // Possibly scale the costs
  //  if (allow_cost_scaling) scaleSimplexCost(options, lp, scale.cost);

  // If matrix is unscaled, then LP is only scaled if there is a cost scaling
  // factor
  //  if (!scaled_matrix) lp.is_scaled_ = scale.cost != 1;
}

bool equilibrationScaleMatrix(const HighsOptions& options, HighsLp& lp,
                              const HighsInt use_scale_strategy) {
  HighsInt numCol = lp.num_col_;
  HighsInt numRow = lp.num_row_;
  HighsScale& scale = lp.scale_;
  vector<double>& colScale = scale.col;
  vector<double>& rowScale = scale.row;
  vector<HighsInt>& Astart = lp.a_matrix_.start_;
  vector<HighsInt>& Aindex = lp.a_matrix_.index_;
  vector<double>& Avalue = lp.a_matrix_.value_;
  vector<double>& colCost = lp.col_cost_;

  HighsInt simplex_scale_strategy = use_scale_strategy;

  double original_matrix_min_value = kHighsInf;
  double original_matrix_max_value = 0;
  for (HighsInt k = 0, AnX = Astart[numCol]; k < AnX; k++) {
    double value = fabs(Avalue[k]);
    original_matrix_min_value = min(original_matrix_min_value, value);
    original_matrix_max_value = max(original_matrix_max_value, value);
  }

  // Include cost in scaling if minimum nonzero cost is less than 0.1
  double min_nonzero_cost = kHighsInf;
  for (HighsInt i = 0; i < numCol; i++) {
    if (colCost[i]) min_nonzero_cost = min(fabs(colCost[i]), min_nonzero_cost);
  }
  bool include_cost_in_scaling = false;
  include_cost_in_scaling = min_nonzero_cost < 0.1;

  // Limits on scaling factors
  double max_allow_scale;
  double min_allow_scale;
  // Now that kHighsInf =
  // std::numeric_limits<double>::infinity(), this Qi-trick doesn't
  // work so, in recognition, use the old value of kHighsInf
  const double finite_infinity = 1e200;
  max_allow_scale = pow(2.0, options.allowed_matrix_scale_factor);
  min_allow_scale = 1 / max_allow_scale;

  double min_allow_col_scale = min_allow_scale;
  double max_allow_col_scale = max_allow_scale;
  double min_allow_row_scale = min_allow_scale;
  double max_allow_row_scale = max_allow_scale;

  // Search up to 6 times
  vector<double> row_min_value(numRow, finite_infinity);
  vector<double> row_max_value(numRow, 1 / finite_infinity);
  for (HighsInt search_count = 0; search_count < 6; search_count++) {
    // Find column scale, prepare row data
    for (HighsInt iCol = 0; iCol < numCol; iCol++) {
      // For column scale (find)
      double col_min_value = finite_infinity;
      double col_max_value = 1 / finite_infinity;
      double abs_col_cost = fabs(colCost[iCol]);
      if (include_cost_in_scaling && abs_col_cost != 0) {
        col_min_value = min(col_min_value, abs_col_cost);
        col_max_value = max(col_max_value, abs_col_cost);
      }
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        double value = fabs(Avalue[k]) * rowScale[Aindex[k]];
        col_min_value = min(col_min_value, value);
        col_max_value = max(col_max_value, value);
      }
      double col_equilibration = 1 / sqrt(col_min_value * col_max_value);
      // Ensure that column scale factor is not excessively large or small
      colScale[iCol] =
          min(max(min_allow_col_scale, col_equilibration), max_allow_col_scale);
      // For row scale (only collect)
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        double value = fabs(Avalue[k]) * colScale[iCol];
        row_min_value[iRow] = min(row_min_value[iRow], value);
        row_max_value[iRow] = max(row_max_value[iRow], value);
      }
    }
    // For row scale (find)
    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
      double row_equilibration =
          1 / sqrt(row_min_value[iRow] * row_max_value[iRow]);
      // Ensure that row scale factor is not excessively large or small
      rowScale[iRow] =
          min(max(min_allow_row_scale, row_equilibration), max_allow_row_scale);
    }
    row_min_value.assign(numRow, finite_infinity);
    row_max_value.assign(numRow, 1 / finite_infinity);
  }
  // Make it numerically better
  // Also determine the max and min row and column scaling factors
  double min_col_scale = finite_infinity;
  double max_col_scale = 1 / finite_infinity;
  double min_row_scale = finite_infinity;
  double max_row_scale = 1 / finite_infinity;
  const double log2 = log(2.0);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / log2 + 0.5));
    min_col_scale = min(colScale[iCol], min_col_scale);
    max_col_scale = max(colScale[iCol], max_col_scale);
  }
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    rowScale[iRow] = pow(2.0, floor(log(rowScale[iRow]) / log2 + 0.5));
    min_row_scale = min(rowScale[iRow], min_row_scale);
    max_row_scale = max(rowScale[iRow], max_row_scale);
  }
  // Apply scaling to matrix and bounds
  double matrix_min_value = finite_infinity;
  double matrix_max_value = 0;
  double min_original_col_equilibration = finite_infinity;
  double sum_original_log_col_equilibration = 0;
  double max_original_col_equilibration = 0;
  double min_original_row_equilibration = finite_infinity;
  double sum_original_log_row_equilibration = 0;
  double max_original_row_equilibration = 0;
  double min_col_equilibration = finite_infinity;
  double sum_log_col_equilibration = 0;
  double max_col_equilibration = 0;
  double min_row_equilibration = finite_infinity;
  double sum_log_row_equilibration = 0;
  double max_row_equilibration = 0;
  vector<double> original_row_min_value(numRow, finite_infinity);
  vector<double> original_row_max_value(numRow, 1 / finite_infinity);
  row_min_value.assign(numRow, finite_infinity);
  row_max_value.assign(numRow, 1 / finite_infinity);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    double original_col_min_value = finite_infinity;
    double original_col_max_value = 1 / finite_infinity;
    double col_min_value = finite_infinity;
    double col_max_value = 1 / finite_infinity;
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      HighsInt iRow = Aindex[k];
      const double original_value = fabs(Avalue[k]);
      original_col_min_value = min(original_value, original_col_min_value);
      original_col_max_value = max(original_value, original_col_max_value);
      original_row_min_value[iRow] =
          min(original_row_min_value[iRow], original_value);
      original_row_max_value[iRow] =
          max(original_row_max_value[iRow], original_value);
      Avalue[k] *= (colScale[iCol] * rowScale[iRow]);
      const double value = fabs(Avalue[k]);
      col_min_value = min(value, col_min_value);
      col_max_value = max(value, col_max_value);
      row_min_value[iRow] = min(row_min_value[iRow], value);
      row_max_value[iRow] = max(row_max_value[iRow], value);
    }
    matrix_min_value = min(matrix_min_value, col_min_value);
    matrix_max_value = max(matrix_max_value, col_max_value);

    const double original_col_equilibration =
        1 / sqrt(original_col_min_value * original_col_max_value);
    min_original_col_equilibration =
        min(original_col_equilibration, min_original_col_equilibration);
    sum_original_log_col_equilibration += log(original_col_equilibration);
    max_original_col_equilibration =
        max(original_col_equilibration, max_original_col_equilibration);
    const double col_equilibration = 1 / sqrt(col_min_value * col_max_value);
    min_col_equilibration = min(col_equilibration, min_col_equilibration);
    sum_log_col_equilibration += log(col_equilibration);
    max_col_equilibration = max(col_equilibration, max_col_equilibration);
  }

  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    const double original_row_equilibration =
        1 / sqrt(original_row_min_value[iRow] * original_row_max_value[iRow]);
    min_original_row_equilibration =
        min(original_row_equilibration, min_original_row_equilibration);
    sum_original_log_row_equilibration += log(original_row_equilibration);
    max_original_row_equilibration =
        max(original_row_equilibration, max_original_row_equilibration);
    const double row_equilibration =
        1 / sqrt(row_min_value[iRow] * row_max_value[iRow]);
    min_row_equilibration = min(row_equilibration, min_row_equilibration);
    sum_log_row_equilibration += log(row_equilibration);
    max_row_equilibration = max(row_equilibration, max_row_equilibration);
  }
  const double geomean_original_col_equilibration =
      exp(sum_original_log_col_equilibration / numCol);
  const double geomean_original_row_equilibration =
      exp(sum_original_log_row_equilibration / numRow);
  const double geomean_col_equilibration =
      exp(sum_log_col_equilibration / numCol);
  const double geomean_row_equilibration =
      exp(sum_log_row_equilibration / numRow);
  if (options.log_dev_level) {
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Original equilibration: min/mean/max %11.4g/%11.4g/%11.4g "
        "(cols); min/mean/max %11.4g/%11.4g/%11.4g (rows)\n",
        min_original_col_equilibration, geomean_original_col_equilibration,
        max_original_col_equilibration, min_original_row_equilibration,
        geomean_original_row_equilibration, max_original_row_equilibration);
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Final    equilibration: min/mean/max %11.4g/%11.4g/%11.4g "
        "(cols); min/mean/max %11.4g/%11.4g/%11.4g (rows)\n",
        min_col_equilibration, geomean_col_equilibration, max_col_equilibration,
        min_row_equilibration, geomean_row_equilibration,
        max_row_equilibration);
  }

  // Compute the mean equilibration improvement
  const double geomean_original_col =
      max(geomean_original_col_equilibration,
          1 / geomean_original_col_equilibration);
  const double geomean_original_row =
      max(geomean_original_row_equilibration,
          1 / geomean_original_row_equilibration);
  const double geomean_col =
      max(geomean_col_equilibration, 1 / geomean_col_equilibration);
  const double geomean_row =
      max(geomean_row_equilibration, 1 / geomean_row_equilibration);
  const double mean_equilibration_improvement =
      sqrt((geomean_original_col * geomean_original_row) /
           (geomean_col * geomean_row));
  // Compute the extreme equilibration improvement
  const double original_col_ratio =
      max_original_col_equilibration / min_original_col_equilibration;
  const double original_row_ratio =
      max_original_row_equilibration / min_original_row_equilibration;
  const double col_ratio = max_col_equilibration / min_col_equilibration;
  const double row_ratio = max_row_equilibration / min_row_equilibration;
  const double extreme_equilibration_improvement =
      (original_col_ratio + original_row_ratio) / (col_ratio + row_ratio);
  // Compute the max/min matrix value improvement
  const double matrix_value_ratio = matrix_max_value / matrix_min_value;
  const double original_matrix_value_ratio =
      original_matrix_max_value / original_matrix_min_value;
  const double matrix_value_ratio_improvement =
      original_matrix_value_ratio / matrix_value_ratio;
  if (options.log_dev_level) {
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Extreme equilibration improvement =      ( %11.4g + "
        "%11.4g) / ( %11.4g + %11.4g)  =      %11.4g / %11.4g  = %11.4g\n",
        original_col_ratio, original_row_ratio, col_ratio, row_ratio,
        (original_col_ratio + original_row_ratio), (col_ratio + row_ratio),
        extreme_equilibration_improvement);
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Mean    equilibration improvement = sqrt(( %11.4g * "
        "%11.4g) / ( %11.4g * %11.4g)) = sqrt(%11.4g / %11.4g) = %11.4g\n",
        geomean_original_col, geomean_original_row, geomean_col, geomean_row,
        (geomean_original_col * geomean_original_row),
        (geomean_col * geomean_row), mean_equilibration_improvement);
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Yields [min, max, ratio] matrix values of [%0.4g, %0.4g, "
        "%0.4g]; Originally [%0.4g, %0.4g, %0.4g]: Improvement of %0.4g\n",
        matrix_min_value, matrix_max_value, matrix_value_ratio,
        original_matrix_min_value, original_matrix_max_value,
        original_matrix_value_ratio, matrix_value_ratio_improvement);
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "Scaling: Improves    mean equilibration by a factor %0.4g\n",
                mean_equilibration_improvement);
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "Scaling: Improves extreme equilibration by a factor %0.4g\n",
                extreme_equilibration_improvement);
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "Scaling: Improves max/min matrix values by a factor %0.4g\n",
                matrix_value_ratio_improvement);
  }
  const bool possibly_abandon_scaling =
      simplex_scale_strategy != kSimplexScaleStrategyForcedEquilibration;
  const double improvement_factor = extreme_equilibration_improvement *
                                    mean_equilibration_improvement *
                                    matrix_value_ratio_improvement;

  const double improvement_factor_required = 1.0;
  const bool poor_improvement =
      improvement_factor < improvement_factor_required;

  // Possibly abandon scaling if it's not improved equilibration significantly
  if (possibly_abandon_scaling && poor_improvement) {
    // Unscale the matrix
    for (HighsInt iCol = 0; iCol < numCol; iCol++) {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        Avalue[k] /= (colScale[iCol] * rowScale[iRow]);
      }
    }
    if (options.log_dev_level)
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Improvement factor %0.4g < %0.4g required, so no "
                  "scaling applied\n",
                  improvement_factor, improvement_factor_required);
    return false;
  } else {
    if (options.log_dev_level) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Factors are in [%0.4g, %0.4g] for columns and in "
                  "[%0.4g, %0.4g] for rows\n",
                  min_col_scale, max_col_scale, min_row_scale, max_row_scale);
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Improvement factor is %0.4g >= %0.4g so scale LP\n",
                  improvement_factor, improvement_factor_required);
      if (extreme_equilibration_improvement < 1.0) {
        highsLogDev(
            options.log_options, HighsLogType::kWarning,
            "Scaling: Applying scaling with extreme improvement of %0.4g\n",
            extreme_equilibration_improvement);
      }
      if (mean_equilibration_improvement < 1.0) {
        highsLogDev(
            options.log_options, HighsLogType::kWarning,
            "Scaling: Applying scaling with mean improvement of %0.4g\n",
            mean_equilibration_improvement);
      }
      if (matrix_value_ratio_improvement < 1.0) {
        highsLogDev(options.log_options, HighsLogType::kWarning,
                    "Scaling: Applying scaling with matrix value ratio "
                    "improvement of %0.4g\n",
                    matrix_value_ratio_improvement);
      }
      if (improvement_factor < 10 * improvement_factor_required) {
        highsLogDev(options.log_options, HighsLogType::kWarning,
                    "Scaling: Applying scaling with improvement factor %0.4g "
                    "< 10*(%0.4g) improvement\n",
                    improvement_factor, improvement_factor_required);
      }
    }
  }
  return true;
}

bool maxValueScaleMatrix(const HighsOptions& options, HighsLp& lp,
                         const HighsInt use_scale_strategy) {
  HighsInt numCol = lp.num_col_;
  HighsInt numRow = lp.num_row_;
  HighsScale& scale = lp.scale_;
  vector<double>& colScale = scale.col;
  vector<double>& rowScale = scale.row;
  vector<HighsInt>& Astart = lp.a_matrix_.start_;
  vector<HighsInt>& Aindex = lp.a_matrix_.index_;
  vector<double>& Avalue = lp.a_matrix_.value_;

  assert(options.simplex_scale_strategy == kSimplexScaleStrategyMaxValue);
  assert(kSimplexScaleStrategyMaxValue015 == kSimplexScaleStrategyMaxValue);
  assert(kSimplexScaleStrategyMaxValue0157 == kSimplexScaleStrategyMaxValue);

  // The 015(7) values refer to bit settings in FICO's scaling options.
  // Specifically
  //
  // 0: Row scaling
  //
  // 1: Column scaling
  //
  // 5: Scale by maximum element
  //
  // 7: Scale objective function for the simplex method
  //
  // Note that 7 is not yet implemented, so
  // kSimplexScaleStrategyMaxValue015 and
  // kSimplexScaleStrategyMaxValue0157 are equivalent. However, cost
  // scaling could be well worth adding, now that the unscaled problem
  // can be solved using scaled NLA

  const double log2 = log(2.0);
  const double max_allow_scale = pow(2.0, options.allowed_matrix_scale_factor);
  const double min_allow_scale = 1 / max_allow_scale;

  const double min_allow_col_scale = min_allow_scale;
  const double max_allow_col_scale = max_allow_scale;
  const double min_allow_row_scale = min_allow_scale;
  const double max_allow_row_scale = max_allow_scale;

  double min_row_scale = kHighsInf;
  double max_row_scale = 0;
  double original_matrix_min_value = kHighsInf;
  double original_matrix_max_value = 0;
  // Determine the row scaling. Also determine the max/min row scaling
  // factors, and max/min original matrix values
  vector<double> row_max_value(numRow, 0);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      const HighsInt iRow = Aindex[k];
      const double value = fabs(Avalue[k]);
      row_max_value[iRow] = max(row_max_value[iRow], value);
      original_matrix_min_value = min(original_matrix_min_value, value);
      original_matrix_max_value = max(original_matrix_max_value, value);
    }
  }
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    if (row_max_value[iRow]) {
      double row_scale_value = 1 / row_max_value[iRow];
      // Convert the row scale factor to the nearest power of two, and
      // ensure that it is not excessively large or small
      row_scale_value = pow(2.0, floor(log(row_scale_value) / log2 + 0.5));
      row_scale_value =
          min(max(min_allow_row_scale, row_scale_value), max_allow_row_scale);
      min_row_scale = min(row_scale_value, min_row_scale);
      max_row_scale = max(row_scale_value, max_row_scale);
      rowScale[iRow] = row_scale_value;
    }
  }
  // Determine the column scaling, whilst applying the row scaling
  // Also determine the max/min column scaling factors, and max/min
  // matrix values
  double min_col_scale = kHighsInf;
  double max_col_scale = 0;
  double matrix_min_value = kHighsInf;
  double matrix_max_value = 0;
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    double col_max_value = 0;
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      const HighsInt iRow = Aindex[k];
      Avalue[k] *= rowScale[iRow];
      const double value = fabs(Avalue[k]);
      col_max_value = max(col_max_value, value);
    }
    if (col_max_value) {
      double col_scale_value = 1 / col_max_value;
      // Convert the col scale factor to the nearest power of two, and
      // ensure that it is not excessively large or small
      col_scale_value = pow(2.0, floor(log(col_scale_value) / log2 + 0.5));
      col_scale_value =
          min(max(min_allow_col_scale, col_scale_value), max_allow_col_scale);
      min_col_scale = min(col_scale_value, min_col_scale);
      max_col_scale = max(col_scale_value, max_col_scale);
      colScale[iCol] = col_scale_value;
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        Avalue[k] *= colScale[iCol];
        const double value = fabs(Avalue[k]);
        matrix_min_value = min(matrix_min_value, value);
        matrix_max_value = max(matrix_max_value, value);
      }
    }
  }
  const double matrix_value_ratio = matrix_max_value / matrix_min_value;
  const double original_matrix_value_ratio =
      original_matrix_max_value / original_matrix_min_value;
  const double matrix_value_ratio_improvement =
      original_matrix_value_ratio / matrix_value_ratio;

  const double improvement_factor = matrix_value_ratio_improvement;

  const double improvement_factor_required = 1.0;
  const bool poor_improvement =
      improvement_factor <= improvement_factor_required;

  if (poor_improvement) {
    // Unscale the matrix
    for (HighsInt iCol = 0; iCol < numCol; iCol++) {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        Avalue[k] /= (colScale[iCol] * rowScale[iRow]);
      }
    }
    if (options.log_dev_level)
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Improvement factor %0.4g < %0.4g required, so no "
                  "scaling applied\n",
                  improvement_factor, improvement_factor_required);
    return false;
  } else {
    if (options.log_dev_level) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Factors are in [%0.4g, %0.4g] for columns and in "
                  "[%0.4g, %0.4g] for rows\n",
                  min_col_scale, max_col_scale, min_row_scale, max_row_scale);
      highsLogDev(
          options.log_options, HighsLogType::kInfo,
          "Scaling: Yields [min, max, ratio] matrix values of [%0.4g, %0.4g, "
          "%0.4g]; Originally [%0.4g, %0.4g, %0.4g]: Improvement of %0.4g\n",
          matrix_min_value, matrix_max_value, matrix_value_ratio,
          original_matrix_min_value, original_matrix_max_value,
          original_matrix_value_ratio, matrix_value_ratio_improvement);
    }
    return true;
  }
}

HighsStatus applyScalingToLpCol(HighsLp& lp, const HighsInt col,
                                const double colScale) {
  if (col < 0) return HighsStatus::kError;
  if (col >= lp.num_col_) return HighsStatus::kError;
  if (!colScale) return HighsStatus::kError;

  lp.a_matrix_.scaleCol(col, colScale);
  lp.col_cost_[col] *= colScale;
  if (colScale > 0) {
    lp.col_lower_[col] /= colScale;
    lp.col_upper_[col] /= colScale;
  } else {
    const double new_upper = lp.col_lower_[col] / colScale;
    lp.col_lower_[col] = lp.col_upper_[col] / colScale;
    lp.col_upper_[col] = new_upper;
  }
  return HighsStatus::kOk;
}

HighsStatus applyScalingToLpRow(HighsLp& lp, const HighsInt row,
                                const double rowScale) {
  if (row < 0) return HighsStatus::kError;
  if (row >= lp.num_row_) return HighsStatus::kError;
  if (!rowScale) return HighsStatus::kError;

  lp.a_matrix_.scaleRow(row, rowScale);
  if (rowScale > 0) {
    lp.row_lower_[row] *= rowScale;
    lp.row_upper_[row] *= rowScale;
  } else {
    const double new_upper = lp.row_lower_[row] * rowScale;
    lp.row_lower_[row] = lp.row_upper_[row] * rowScale;
    lp.row_upper_[row] = new_upper;
  }
  return HighsStatus::kOk;
}

void unscaleSolution(HighsSolution& solution, const HighsScale& scale) {
  assert(scale.has_scaling);
  assert(solution.col_value.size() == static_cast<size_t>(scale.num_col));
  assert(solution.col_dual.size() == static_cast<size_t>(scale.num_col));
  assert(solution.row_value.size() == static_cast<size_t>(scale.num_row));
  assert(solution.row_dual.size() == static_cast<size_t>(scale.num_row));

  for (HighsInt iCol = 0; iCol < scale.num_col; iCol++) {
    solution.col_value[iCol] *= scale.col[iCol];
    solution.col_dual[iCol] /= (scale.col[iCol] / scale.cost);
  }
  for (HighsInt iRow = 0; iRow < scale.num_row; iRow++) {
    solution.row_value[iRow] /= scale.row[iRow];
    solution.row_dual[iRow] *= (scale.row[iRow] * scale.cost);
  }
}

void appendColsToLpVectors(HighsLp& lp, const HighsInt num_new_col,
                           const vector<double>& colCost,
                           const vector<double>& colLower,
                           const vector<double>& colUpper) {
  assert(num_new_col >= 0);
  if (num_new_col == 0) return;
  HighsInt new_num_col = lp.num_col_ + num_new_col;
  lp.col_cost_.resize(new_num_col);
  lp.col_lower_.resize(new_num_col);
  lp.col_upper_.resize(new_num_col);
  const bool have_integrality = (lp.integrality_.size() != 0);
  if (have_integrality) {
    assert(HighsInt(lp.integrality_.size()) == lp.num_col_);
    lp.integrality_.resize(new_num_col);
  }
  bool have_names = (lp.col_names_.size() != 0);
  if (have_names) lp.col_names_.resize(new_num_col);
  for (HighsInt new_col = 0; new_col < num_new_col; new_col++) {
    HighsInt iCol = lp.num_col_ + new_col;
    lp.col_cost_[iCol] = colCost[new_col];
    lp.col_lower_[iCol] = colLower[new_col];
    lp.col_upper_[iCol] = colUpper[new_col];
    if (have_names) {
      // Cannot guarantee to create unique names, so name is blank
      lp.col_hash_.addName(iCol, "");
      lp.col_names_[iCol] = "";
    }
    if (have_integrality) lp.integrality_[iCol] = HighsVarType::kContinuous;
  }
}

void appendRowsToLpVectors(HighsLp& lp, const HighsInt num_new_row,
                           const vector<double>& rowLower,
                           const vector<double>& rowUpper) {
  assert(num_new_row >= 0);
  if (num_new_row == 0) return;
  HighsInt new_num_row = lp.num_row_ + num_new_row;
  lp.row_lower_.resize(new_num_row);
  lp.row_upper_.resize(new_num_row);
  bool have_names = (lp.row_names_.size() != 0);
  if (have_names) lp.row_names_.resize(new_num_row);

  for (HighsInt new_row = 0; new_row < num_new_row; new_row++) {
    HighsInt iRow = lp.num_row_ + new_row;
    lp.row_lower_[iRow] = rowLower[new_row];
    lp.row_upper_[iRow] = rowUpper[new_row];
    if (have_names) {
      // Cannot guarantee to create unique names, so name is blank
      lp.row_hash_.addName(iRow, "");
      lp.row_names_[iRow] = "";
    }
  }
}

void deleteScale(vector<double>& scale,
                 const HighsIndexCollection& index_collection) {
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return;

  HighsInt delete_from_col;
  HighsInt delete_to_col;
  HighsInt keep_from_col;
  HighsInt keep_to_col = -1;
  HighsInt current_set_entry = 0;

  HighsInt col_dim = index_collection.dimension_;
  HighsInt new_num_col = 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_col, delete_to_col,
                     keep_from_col, keep_to_col, current_set_entry);
    // Account for the initial columns being kept
    if (k == from_k) new_num_col = delete_from_col;
    if (delete_to_col >= col_dim - 1) break;
    assert(delete_to_col < col_dim);
    for (HighsInt col = keep_from_col; col <= keep_to_col; col++) {
      assert((HighsInt)scale.size() > new_num_col);
      scale[new_num_col] = scale[col];
      new_num_col++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
}

HighsStatus changeLpMatrixCoeff(const HighsLogOptions& log_options, HighsLp& lp,
                                const HighsInt row, const HighsInt col,
                                const double new_value,
                                const double small_matrix_value) {
  assert(0 <= row && row < lp.num_row_);
  assert(0 <= col && col < lp.num_col_);

  const double abs_new_value = std::fabs(new_value);
  const bool zero_new_value = abs_new_value <= small_matrix_value;
  // Determine whether the coefficient corresponds to an existing
  // nonzero
  HighsInt change_el = -1;
  for (HighsInt el = lp.a_matrix_.start_[col];
       el < lp.a_matrix_.start_[col + 1]; el++) {
    if (lp.a_matrix_.index_[el] == row) {
      change_el = el;
      break;
    }
  }
  HighsInt new_num_nz = lp.numNz();
  if (change_el < 0) {
    // Coefficient doesn't correspond to an existing nonzero
    //
    // If coefficient is small, then just ignore it
    if (zero_new_value) {
      highsLogUser(
          log_options, HighsLogType::kWarning,
          "|Value| of %g supplied to Highs::changeCoeff is less than %g and "
          "does not correspond to an existing nonzero: ignored\n",
          abs_new_value, small_matrix_value);
      return HighsStatus::kWarning;
    }
    // New nonzero goes at the end of column "col", so have to
    // increment all subsequent start values by 1, and shift all
    // subsequent index and value entries forward by 1 to accommodate
    // it
    change_el = lp.a_matrix_.start_[col + 1];
    new_num_nz++;
    lp.a_matrix_.index_.resize(new_num_nz);
    lp.a_matrix_.value_.resize(new_num_nz);
    for (HighsInt i = col + 1; i <= lp.num_col_; i++) lp.a_matrix_.start_[i]++;
    for (HighsInt el = new_num_nz - 1; el > change_el; el--) {
      lp.a_matrix_.index_[el] = lp.a_matrix_.index_[el - 1];
      lp.a_matrix_.value_[el] = lp.a_matrix_.value_[el - 1];
    }
  } else if (zero_new_value) {
    // Coefficient zeroes an existing nonzero
    HighsStatus return_status = HighsStatus::kOk;
    if (abs_new_value > 0) {
      // Warn if the new value was small but nonzero, so treated as
      // zero. If the new value is zero, it's OK, as it's the way of
      // zeroing nonzero entries.
      highsLogUser(log_options, HighsLogType::kWarning,
                   "|Value| of %g supplied to Highs::changeCoeff is less than "
                   "%g and zeroes (%d, %d) entry of %g\n",
                   abs_new_value, small_matrix_value, int(row), int(col),
                   lp.a_matrix_.value_[change_el]);
      return_status = HighsStatus::kWarning;
    }
    // Zeroing an existing nonzero, so have to decrement all
    // subsequent start values by 1, and shift all subsequent index
    // and value entries backward by 1 to eliminate it
    new_num_nz--;
    for (HighsInt i = col + 1; i <= lp.num_col_; i++) lp.a_matrix_.start_[i]--;
    for (HighsInt el = change_el; el < new_num_nz; el++) {
      lp.a_matrix_.index_[el] = lp.a_matrix_.index_[el + 1];
      lp.a_matrix_.value_[el] = lp.a_matrix_.value_[el + 1];
    }
    lp.a_matrix_.index_.resize(new_num_nz);
    lp.a_matrix_.value_.resize(new_num_nz);
    return return_status;
  }
  assert(change_el >= 0 && !zero_new_value);
  lp.a_matrix_.index_[change_el] = row;
  lp.a_matrix_.value_[change_el] = new_value;
  return HighsStatus::kOk;
}

HighsStatus changeLpIntegrality(HighsLp& lp,
                                const HighsIndexCollection& index_collection,
                                const vector<HighsVarType>& new_integrality,
                                const HighsOptions& options) {
  HighsStatus return_status = HighsStatus::kOk;
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return return_status;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const vector<HighsInt>& col_set = index_collection.set_;
  const vector<HighsInt>& col_mask = index_collection.mask_;

  auto passThroughNewIntegrality = [&](const bool find_non_continuous) {
    HighsInt lp_col;
    HighsInt usr_col = -1;
    for (HighsInt k = from_k; k < to_k + 1; k++) {
      if (interval || mask) {
        lp_col = k;
      } else {
        lp_col = col_set[k];
      }
      HighsInt col = lp_col;
      if (interval) {
        usr_col++;
      } else {
        usr_col = k;
      }
      if (mask && !col_mask[col]) continue;
      if (find_non_continuous) {
        if (new_integrality[usr_col] != HighsVarType::kContinuous) return true;
      } else {
        lp.integrality_[col] = new_integrality[usr_col];
      }
    }
    return false;
  };

  auto hasNonContinuous = [&]() { return passThroughNewIntegrality(true); };
  auto assignNewIntegrality = [&]() { passThroughNewIntegrality(false); };

  // If changing integrality for a problem without an integrality
  // vector (ie an LP), check whether it is needed before creating it
  // for the incumbent columns - which are naturally continuous -
  // otherwise (#3476) n repeated calls setting single columns to
  // HighsVarType::kContinuous costs O(n^2)
  if (lp.integrality_.size() == 0) {
    if (!hasNonContinuous()) return return_status;
    lp.integrality_.assign(lp.num_col_, HighsVarType::kContinuous);
  }
  assert(HighsInt(lp.integrality_.size()) == lp.num_col_);
  // Change the integrality to the user-supplied integrality, according to the
  // technique
  assignNewIntegrality();

  // If integrality_ now contains only HighsVarType::kContinuous then
  // clear it
  if (!lp.isMip()) lp.integrality_.clear();
  return return_status;
}

void changeLpCosts(HighsLp& lp, const HighsIndexCollection& index_collection,
                   const vector<double>& new_col_cost,
                   const double infinite_cost) {
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const vector<HighsInt>& col_set = index_collection.set_;
  const vector<HighsInt>& col_mask = index_collection.mask_;

  // Change the costs to the user-supplied costs, according to the technique
  HighsInt lp_col;
  HighsInt usr_col = -1;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      lp_col = k;
    } else {
      lp_col = col_set[k];
    }
    HighsInt col = lp_col;
    if (interval) {
      usr_col++;
    } else {
      usr_col = k;
    }
    if (mask && !col_mask[col]) continue;
    lp.col_cost_[col] = new_col_cost[usr_col];
  }
  // Check whether the LP still has an infinite cost
  if (lp.has_infinite_cost_)
    lp.has_infinite_cost_ = lp.hasInfiniteCost(infinite_cost);
}

void changeLpColBounds(HighsLp& lp,
                       const HighsIndexCollection& index_collection,
                       const vector<double>& new_col_lower,
                       const vector<double>& new_col_upper) {
  changeBounds(lp.col_lower_, lp.col_upper_, index_collection, new_col_lower,
               new_col_upper);
}

void changeLpRowBounds(HighsLp& lp,
                       const HighsIndexCollection& index_collection,
                       const vector<double>& new_row_lower,
                       const vector<double>& new_row_upper) {
  changeBounds(lp.row_lower_, lp.row_upper_, index_collection, new_row_lower,
               new_row_upper);
}

void changeBounds(vector<double>& lower, vector<double>& upper,
                  const HighsIndexCollection& index_collection,
                  const vector<double>& new_lower,
                  const vector<double>& new_upper) {
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const vector<HighsInt>& ix_set = index_collection.set_;
  const vector<HighsInt>& ix_mask = index_collection.mask_;

  // Change the bounds to the user-supplied bounds, according to the technique
  HighsInt lp_ix;
  HighsInt usr_ix = -1;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      lp_ix = k;
    } else {
      lp_ix = ix_set[k];
    }
    HighsInt ix = lp_ix;
    if (interval) {
      usr_ix++;
    } else {
      usr_ix = k;
    }
    if (mask && !ix_mask[ix]) continue;
    lower[ix] = new_lower[usr_ix];
    upper[ix] = new_upper[usr_ix];
  }
}

HighsInt getNumInt(const HighsLp& lp) {
  HighsInt num_int = 0;
  if (lp.integrality_.size()) {
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
      if (lp.integrality_[iCol] == HighsVarType::kInteger) num_int++;
  }
  return num_int;
}

void getLpCosts(const HighsLp& lp, const HighsInt from_col,
                const HighsInt to_col, double* XcolCost) {
  assert(0 <= from_col && from_col < lp.num_col_);
  assert(0 <= to_col && to_col < lp.num_col_);
  if (from_col > to_col) return;
  for (HighsInt col = from_col; col < to_col + 1; col++)
    XcolCost[col - from_col] = lp.col_cost_[col];
}

void getLpColBounds(const HighsLp& lp, const HighsInt from_col,
                    const HighsInt to_col, double* XcolLower,
                    double* XcolUpper) {
  assert(0 <= from_col && from_col < lp.num_col_);
  assert(0 <= to_col && to_col < lp.num_col_);
  if (from_col > to_col) return;
  for (HighsInt col = from_col; col < to_col + 1; col++) {
    if (XcolLower != nullptr) XcolLower[col - from_col] = lp.col_lower_[col];
    if (XcolUpper != nullptr) XcolUpper[col - from_col] = lp.col_upper_[col];
  }
}

void getLpRowBounds(const HighsLp& lp, const HighsInt from_row,
                    const HighsInt to_row, double* XrowLower,
                    double* XrowUpper) {
  assert(0 <= from_row && from_row < lp.num_row_);
  assert(0 <= to_row && to_row < lp.num_row_);
  if (from_row > to_row) return;
  for (HighsInt row = from_row; row < to_row + 1; row++) {
    if (XrowLower != nullptr) XrowLower[row - from_row] = lp.row_lower_[row];
    if (XrowUpper != nullptr) XrowUpper[row - from_row] = lp.row_upper_[row];
  }
}

// Get a single coefficient from the matrix
void getLpMatrixCoefficient(const HighsLp& lp, const HighsInt Xrow,
                            const HighsInt Xcol, double* val) {
  assert(0 <= Xrow && Xrow < lp.num_row_);
  assert(0 <= Xcol && Xcol < lp.num_col_);

  HighsInt get_el = -1;
  for (HighsInt el = lp.a_matrix_.start_[Xcol];
       el < lp.a_matrix_.start_[Xcol + 1]; el++) {
    if (lp.a_matrix_.index_[el] == Xrow) {
      get_el = el;
      break;
    }
  }
  if (get_el < 0) {
    *val = 0;
  } else {
    *val = lp.a_matrix_.value_[get_el];
  }
}

// Methods for reporting an LP, including its row and column data and matrix
//
// Report the whole LP
void reportLp(const HighsLogOptions& log_options, const HighsLp& lp,
              const HighsLogType report_level) {
  reportLpBrief(log_options, lp);
  if ((HighsInt)report_level >= (HighsInt)HighsLogType::kDetailed) {
    reportLpColVectors(log_options, lp);
    reportLpRowVectors(log_options, lp);
    if ((HighsInt)report_level >= (HighsInt)HighsLogType::kVerbose) {
      assert(lp.a_matrix_.isColwise());
      reportLpColMatrix(log_options, lp);
    }
  }
}

// Report the LP briefly
void reportLpBrief(const HighsLogOptions& log_options, const HighsLp& lp) {
  reportLpDimensions(log_options, lp);
  reportLpObjSense(log_options, lp);
}

// Report the LP dimensions
void reportLpDimensions(const HighsLogOptions& log_options, const HighsLp& lp) {
  HighsInt lp_num_nz;
  if (lp.num_col_ == 0)
    lp_num_nz = 0;
  else
    lp_num_nz = lp.a_matrix_.start_[lp.num_col_];
  highsLogUser(log_options, HighsLogType::kInfo,
               "LP has %" HIGHSINT_FORMAT " row%s, %" HIGHSINT_FORMAT
               " column%s",
               lp.num_row_, lp.num_row_ == 1 ? "" : "s", lp.num_col_,
               lp.num_col_ == 1 ? "" : "s");
  HighsInt num_int = getNumInt(lp);
  if (num_int) {
    highsLogUser(log_options, HighsLogType::kInfo,
                 ", %" HIGHSINT_FORMAT " nonzero%s and %" HIGHSINT_FORMAT
                 " integer column%s\n",
                 lp_num_nz, lp_num_nz == 1 ? "" : "s", num_int,
                 num_int == 1 ? "" : "s");
  } else {
    highsLogUser(log_options, HighsLogType::kInfo,
                 " and %" HIGHSINT_FORMAT " nonzero%s\n", lp_num_nz,
                 lp_num_nz == 1 ? "" : "s");
  }
}

// Report the LP objective sense
void reportLpObjSense(const HighsLogOptions& log_options, const HighsLp& lp) {
  if (lp.sense_ == ObjSense::kMinimize)
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Objective sense is minimize\n");
  else if (lp.sense_ == ObjSense::kMaximize)
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Objective sense is maximize\n");
  else
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Objective sense is ill-defined as %" HIGHSINT_FORMAT "\n",
                 lp.sense_);
}

static std::string getBoundType(const double lower, const double upper) {
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
void reportLpColVectors(const HighsLogOptions& log_options, const HighsLp& lp) {
  if (lp.num_col_ <= 0) return;
  std::string type;
  HighsInt count;
  bool have_integer_columns = (getNumInt(lp) != 0);
  bool have_col_names =
      lp.col_names_.size() == static_cast<size_t>(lp.num_col_);

  highsLogUser(log_options, HighsLogType::kInfo,
               "  Column        Lower        Upper         Cost       "
               "Type        Count");
  if (have_integer_columns)
    highsLogUser(log_options, HighsLogType::kInfo, "  Discrete");
  if (have_col_names) highsLogUser(log_options, HighsLogType::kInfo, "  Name");
  highsLogUser(log_options, HighsLogType::kInfo, "\n");

  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    type = getBoundType(lp.col_lower_[iCol], lp.col_upper_[iCol]);
    count = lp.a_matrix_.start_[iCol + 1] - lp.a_matrix_.start_[iCol];
    highsLogUser(log_options, HighsLogType::kInfo,
                 "%8" HIGHSINT_FORMAT
                 " %12g %12g %12g         %2s %12" HIGHSINT_FORMAT "",
                 iCol, lp.col_lower_[iCol], lp.col_upper_[iCol],
                 lp.col_cost_[iCol], type.c_str(), count);
    if (have_integer_columns) {
      std::string integer_column = "";
      if (lp.integrality_[iCol] == HighsVarType::kInteger) {
        if (lp.col_lower_[iCol] == 0 && lp.col_upper_[iCol] == 1) {
          integer_column = "Binary";
        } else {
          integer_column = "Integer";
        }
      }
      highsLogUser(log_options, HighsLogType::kInfo, "  %-8s",
                   integer_column.c_str());
    }
    if (have_col_names)
      highsLogUser(log_options, HighsLogType::kInfo, "  %-s",
                   lp.col_names_[iCol].c_str());
    highsLogUser(log_options, HighsLogType::kInfo, "\n");
  }
}

// Report the vectors of LP row data
void reportLpRowVectors(const HighsLogOptions& log_options, const HighsLp& lp) {
  if (lp.num_row_ <= 0) return;
  std::string type;
  vector<HighsInt> count;
  bool have_row_names = (lp.row_names_.size() != 0);

  count.resize(lp.num_row_, 0);
  if (lp.num_col_ > 0) {
    for (HighsInt el = 0; el < lp.a_matrix_.start_[lp.num_col_]; el++)
      count[lp.a_matrix_.index_[el]]++;
  }

  highsLogUser(log_options, HighsLogType::kInfo,
               "     Row        Lower        Upper       Type        Count");
  if (have_row_names) highsLogUser(log_options, HighsLogType::kInfo, "  Name");
  highsLogUser(log_options, HighsLogType::kInfo, "\n");

  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    type = getBoundType(lp.row_lower_[iRow], lp.row_upper_[iRow]);
    std::string name = "";
    highsLogUser(log_options, HighsLogType::kInfo,
                 "%8" HIGHSINT_FORMAT
                 " %12g %12g         %2s %12" HIGHSINT_FORMAT "",
                 iRow, lp.row_lower_[iRow], lp.row_upper_[iRow], type.c_str(),
                 count[iRow]);
    if (have_row_names)
      highsLogUser(log_options, HighsLogType::kInfo, "  %-s",
                   lp.row_names_[iRow].c_str());
    highsLogUser(log_options, HighsLogType::kInfo, "\n");
  }
}

// Report the LP column-wise matrix
void reportLpColMatrix(const HighsLogOptions& log_options, const HighsLp& lp) {
  if (lp.num_col_ <= 0) return;
  if (lp.num_row_) {
    // With positive number of rows, can assume that there are index and value
    // vectors to pass
    reportMatrix(log_options, "Column", lp.num_col_,
                 lp.a_matrix_.start_[lp.num_col_], lp.a_matrix_.start_.data(),
                 lp.a_matrix_.index_.data(), lp.a_matrix_.value_.data());
  } else {
    // With no rows, can's assume that there are index and value vectors to pass
    reportMatrix(log_options, "Column", lp.num_col_,
                 lp.a_matrix_.start_[lp.num_col_], lp.a_matrix_.start_.data(),
                 NULL, NULL);
  }
}

void reportMatrix(const HighsLogOptions& log_options,
                  const std::string& message, const HighsInt num_col,
                  const HighsInt num_nz, const HighsInt* start,
                  const HighsInt* index, const double* value) {
  if (num_col <= 0) return;
  highsLogUser(log_options, HighsLogType::kInfo,
               "%-7s Index              Value\n", message.c_str());
  for (HighsInt col = 0; col < num_col; col++) {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "    %8" HIGHSINT_FORMAT " Start   %10" HIGHSINT_FORMAT "\n",
                 col, start[col]);
    HighsInt to_el = (col < num_col - 1 ? start[col + 1] : num_nz);
    for (HighsInt el = start[col]; el < to_el; el++)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "          %8" HIGHSINT_FORMAT " %12g\n", index[el],
                   value[el]);
  }
  highsLogUser(log_options, HighsLogType::kInfo,
               "             Start   %10" HIGHSINT_FORMAT "\n", num_nz);
}

void analyseLp(const HighsLogOptions& log_options, const HighsLp& lp) {
  /*
  vector<double> min_colBound;
  vector<double> min_rowBound;
  vector<double> colRange;
  vector<double> rowRange;
  min_colBound.resize(lp.num_col_);
  min_rowBound.resize(lp.num_row_);
  colRange.resize(lp.num_col_);
  rowRange.resize(lp.num_row_);
  for (HighsInt col = 0; col < lp.num_col_; col++)
    min_colBound[col] = min(fabs(lp.col_lower_[col]), fabs(lp.col_upper_[col]));
  for (HighsInt row = 0; row < lp.num_row_; row++)
    min_rowBound[row] = min(fabs(lp.row_lower_[row]), fabs(lp.row_upper_[row]));
  for (HighsInt col = 0; col < lp.num_col_; col++)
    colRange[col] = lp.col_upper_[col] - lp.col_lower_[col];
  for (HighsInt row = 0; row < lp.num_row_; row++)
    rowRange[row] = lp.row_upper_[row] - lp.row_lower_[row];
  */
  std::string message;
  if (lp.is_scaled_) {
    message = "Scaled";
  } else {
    message = "Unscaled";
  }
  highsLogDev(log_options, HighsLogType::kInfo, "\n%s model data: Analysis\n",
              message.c_str());
  if (lp.is_scaled_) {
    const HighsScale& scale = lp.scale_;
    analyseVectorValues(&log_options, "Column scaling factors", lp.num_col_,
                        scale.col, true, lp.model_name_);
    analyseVectorValues(&log_options, "Row    scaling factors", lp.num_row_,
                        scale.row, true, lp.model_name_);
  }
  analyseVectorValues(&log_options, "Column costs", lp.num_col_, lp.col_cost_,
                      true, lp.model_name_);
  analyseVectorValues(&log_options, "Column lower bounds", lp.num_col_,
                      lp.col_lower_, true, lp.model_name_);
  analyseVectorValues(&log_options, "Column upper bounds", lp.num_col_,
                      lp.col_upper_, true, lp.model_name_);
  //  analyseVectorValues(&log_options, "Column min abs bound", lp.num_col_,
  //                      min_colBound, true, lp.model_name_);
  //  analyseVectorValues(&log_options, "Column range", lp.num_col_, colRange,
  //  true,
  //                      lp.model_name_);
  analyseVectorValues(&log_options, "Row lower bounds", lp.num_row_,
                      lp.row_lower_, true, lp.model_name_);
  analyseVectorValues(&log_options, "Row upper bounds", lp.num_row_,
                      lp.row_upper_, true, lp.model_name_);
  //  analyseVectorValues(&log_options, "Row min abs bound", lp.num_row_,
  //                      min_rowBound, true, lp.model_name_);
  //  analyseVectorValues(&log_options, "Row range", lp.num_row_, rowRange,
  //  true,
  //                      lp.model_name_);
  analyseVectorValues(&log_options, "Matrix sparsity",
                      lp.a_matrix_.start_[lp.num_col_], lp.a_matrix_.value_,
                      true, lp.model_name_);
  analyseMatrixSparsity(log_options, "Constraint matrix", lp.num_col_,
                        lp.num_row_, lp.a_matrix_.start_, lp.a_matrix_.index_);
  analyseModelBounds(log_options, "Column", lp.num_col_, lp.col_lower_,
                     lp.col_upper_);
  analyseModelBounds(log_options, "Row", lp.num_row_, lp.row_lower_,
                     lp.row_upper_);
}

HighsStatus readSolutionFile(const std::string& filename,
                             const HighsOptions& options, HighsLp& lp,
                             HighsBasis& basis, HighsSolution& solution,
                             const HighsInt style) {
  HighsStatus return_status = HighsStatus::kOk;
  const HighsLogOptions& log_options = options.log_options;
  if (style != kSolutionStyleRaw && style != kSolutionStyleSparse &&
      style != kSolutionStylePartial) {
    highsLogUser(log_options, HighsLogType::kError,
                 "readSolutionFile: Cannot read file of style %d\n",
                 (int)style);
    return HighsStatus::kError;
  }
  std::ifstream in_file(filename);
  if (in_file.fail()) {
    highsLogUser(log_options, HighsLogType::kError,
                 "readSolutionFile: Cannot open readable file \"%s\"\n",
                 filename.c_str());
    return HighsStatus::kError;
  }
  std::string from_method = "readSolutionFile";
  std::string hash;
  std::string keyword;
  std::string value_string;
  std::string qualifier_string;
  std::string name;
  double value;
  HighsInt num_col = -1;
  HighsInt num_row = -1;
  const HighsInt lp_num_col = lp.num_col_;
  const HighsInt lp_num_row = lp.num_row_;
  // Define identifiers for reading in
  HighsSolution read_solution = solution;
  HighsBasis read_basis = basis;
  read_solution.clear();
  read_basis.clear();
  read_solution.col_value.resize(lp_num_col);
  read_solution.row_value.resize(lp_num_row);
  read_solution.col_dual.resize(lp_num_col);
  read_solution.row_dual.resize(lp_num_row);
  read_basis.col_status.resize(lp_num_col);
  read_basis.row_status.resize(lp_num_row);
  std::string section_name;
  if (!readSolutionFileIdIgnoreLineOk(section_name, in_file))
    return readSolutionFileErrorReturn(
        in_file);  // Model (status) or =obj= (value)
  const bool miplib_sol = section_name == "=obj=";
  const bool have_col_names =
      lp.col_names_.size() == static_cast<size_t>(lp_num_col);
  const bool have_row_names =
      lp.row_names_.size() == static_cast<size_t>(lp_num_row);
  if (miplib_sol) {
    // A MIPLIB solution file has nonzero solution values for a subset
    // of the variables identified by name, so there must be column
    // names,
    if (!have_col_names) {
      highsLogUser(log_options, HighsLogType::kError,
                   "readSolutionFile: Cannot read a MIPLIB solution file "
                   "without column names in the model\n");
      return HighsStatus::kError;
    }
  }
  if (have_col_names) {
    // Ensure that the col name hash table has been formed
    if (!lp.col_hash_.size()) lp.col_hash_.form(lp.col_names_);
  }
  if (have_row_names) {
    // Ensure that the row name hash table has been formed
    if (!lp.row_hash_.size()) lp.row_hash_.form(lp.row_names_);
  }
  bool sparse = false;
  bool partial = false;
  if (!miplib_sol) {
    if (!readSolutionFileIgnoreLineOk(in_file))
      return readSolutionFileErrorReturn(in_file);  // Optimal
    if (!readSolutionFileIgnoreLineOk(in_file))
      return readSolutionFileErrorReturn(in_file);  //
    if (!readSolutionFileIgnoreLineOk(in_file))
      return readSolutionFileErrorReturn(in_file);  // # Primal solution values
    if (!readSolutionFileKeywordLineOk(keyword, in_file))
      return readSolutionFileErrorReturn(in_file);
    // Read in the primal solution values: return warning if there is none
    if (keyword == "None")
      return readSolutionFileReturn(HighsStatus::kWarning, solution, basis,
                                    read_solution, read_basis, in_file);
    // If there are primal solution values then keyword is the status
    // and the next line is objective
    if (!readSolutionFileIgnoreLineOk(in_file))
      return readSolutionFileErrorReturn(in_file);  // EOL
    if (!readSolutionFileIgnoreLineOk(in_file))
      return readSolutionFileErrorReturn(in_file);  // Objective
    // Next line, possibly after ignoring lines giving the
    // "Objective_bound"and "MIP_rel_gap", should be "# Columns",
    // correct number and possibly a qualifier string, so can't be
    // read (easily) as a std::ifstream
    std::string column_section_line;
    for (;;) {
      std::getline(in_file, column_section_line);
      std::stringstream column_section_line_ss(column_section_line);
      if (!readSolutionFileColumnHeaderLineOk(hash, keyword, value_string,
                                              num_col, qualifier_string,
                                              column_section_line_ss)) {
        if (hash != "#") {
          // If # hasn't been read and is "Objective_bound" or
          // "MIP_rel_gap", look at the next line
          if (hash == "Objective_bound" || hash == "MIP_rel_gap") continue;
        }
        // Something illegal has been read, or the file has ended
        // after the objective line (#3499) - which is perversely
        // invalid
        highsLogUser(log_options, HighsLogType::kError,
                     "readSolutionFile: Error reading line \"%s\"\n",
                     column_section_line.c_str());
        return readSolutionFileErrorReturn(in_file);
      }
      // "# Columns" has been read correctly, so break
      break;
    }
    assert(keyword == "Columns");
    // The default style parameter is kSolutionStyleRaw, and this
    // still allows sparse or partial files to be read. Recognise
    // sparse or partial files from num_col <= 0, where sparse is the
    // default (for back-compatibility). Doesn't matter if num_col =
    // 0, since there's nothing to read either way
    sparse = num_col <= 0;
    partial = num_col <= 0 && qualifier_string == kHighsPartialString;
    if (partial) sparse = false;
    if (style == kSolutionStyleSparse) assert(sparse);
    if (style == kSolutionStylePartial) assert(partial);
    if (sparse || partial) {
      num_col = -num_col;
      assert(num_col <= lp_num_col);
    } else {
      if (num_col != lp_num_col) {
        highsLogUser(
            log_options, HighsLogType::kError,
            "readSolutionFile: Solution file is for %d columns, not %d\n",
            int(num_col), int(lp_num_col));
        return readSolutionFileErrorReturn(in_file);
      }
    }
  }
  bool is_col = true;
  HighsInt iCol;
  if (miplib_sol) {
    HighsInt num_value = 0;
    read_solution.col_value.assign(lp_num_col, 0);
    assert(is_col);
    for (;;) {
      // Only false return is for encountering EOF
      if (!readSolutionFileIdDoubleLineOk(name, value, in_file)) break;
      return_status = getIndexFromName(log_options, from_method, is_col, name,
                                       lp.col_hash_, iCol, lp.col_names_);
      if (return_status != HighsStatus::kOk) return return_status;
      read_solution.col_value[iCol] = value;
      num_value++;
      if (in_file.eof()) break;
    }
  } else if (sparse || partial) {
    read_solution.col_value.assign(lp_num_col, sparse ? 0.0 : kHighsUndefined);
    for (HighsInt iX = 0; iX < num_col; iX++) {
      // Each line should be a column name, a value and possibly an
      // index - which is used if the model has no column names
      std::string column_section_line;
      std::getline(in_file, column_section_line);
      std::stringstream column_section_line_ss(column_section_line);
      if (column_section_line.empty() ||
          !readSolutionFileColumnLineOk(name, value, iCol,
                                        column_section_line_ss)) {
        highsLogUser(
            log_options, HighsLogType::kError,
            "readSolutionFile: Error reading column section line \"%s\"\n",
            column_section_line.c_str());
        return readSolutionFileErrorReturn(in_file);
      }
      if (have_col_names) {
        // Use the column name if possible
        return_status = getIndexFromName(log_options, from_method, is_col, name,
                                         lp.col_hash_, iCol, lp.col_names_);
        if (return_status != HighsStatus::kOk) return return_status;
      } else if (iCol < 0) {
        // No column names or index, so cannot assign the value
        highsLogUser(log_options, HighsLogType::kError,
                     "readSolutionFile: Error reading sparse/partial column "
                     "value line \"%s\": column "
                     "index is undefined and model has no column names so name "
                     "\"%s\" cannot be used to assign the value\n",
                     column_section_line.c_str(), name.c_str());
        return readSolutionFileErrorReturn(in_file);
      }
      assert(0 <= iCol && iCol < lp.num_col_);
      read_solution.col_value[iCol] = value;
    }
  } else {
    assert(is_col);
    for (HighsInt iX = 0; iX < num_col; iX++) {
      if (!readSolutionFileIdDoubleLineOk(name, value, in_file)) {
        highsLogUser(
            log_options, HighsLogType::kError,
            "readSolutionFile: Error reading column value line \"%s %g\"\n",
            name.c_str(), value);
        return readSolutionFileErrorReturn(in_file);
      }
      if (have_col_names) {
        // Use the column name if possible
        return_status = getIndexFromName(log_options, from_method, is_col, name,
                                         lp.col_hash_, iCol, lp.col_names_);
        if (return_status != HighsStatus::kOk) return return_status;
      } else {
        // Have to assume column solution values are in the right order
        iCol = iX;
      }
      read_solution.col_value[iCol] = value;
    }
  }
  read_solution.value_valid = true;
  if (miplib_sol || sparse || partial) {
    // MIPLIB, sparse and partial format only have column values, so
    // compute the row values and return
    if (calculateRowValuesQuad(lp, read_solution.col_value,
                               read_solution.row_value) != HighsStatus::kOk)
      return readSolutionFileErrorReturn(in_file);
    return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                  read_solution, read_basis, in_file);
  }
  // Read in the row values: OK to have none, otherwise next line
  // should be "Rows" and correct number
  if (!readSolutionFileHashKeywordIntLineOk(hash, keyword, value_string,
                                            num_row, in_file)) {
    // Compute the row values since there are none to read
    if (calculateRowValuesQuad(lp, read_solution.col_value,
                               read_solution.row_value) != HighsStatus::kOk)
      return readSolutionFileErrorReturn(in_file);
    return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                  read_solution, read_basis, in_file);
  }
  assert(keyword == "Rows");
  // OK to read from a file with different number of rows, since the
  // primal solution is all that's important. For example, see #1284,
  // where the user is solving a sequence of MIPs with the same number
  // of variables, but increasing numbers of constraints, and wants to
  // used the solution from one MIP as the starting solution for the
  // next.
  HighsInt iRow;
  const bool num_row_ok = num_row == lp_num_row;
  assert(!have_row_names || lp.row_hash_.size());
  is_col = false;
  assert(!is_col);
  for (HighsInt iX = 0; iX < num_row; iX++) {
    if (!readSolutionFileIdDoubleLineOk(name, value, in_file)) {
      highsLogUser(log_options, HighsLogType::kError,
                   "readSolutionFile: Error reading row value line \"%s %g\"\n",
                   name.c_str(), value);
      return readSolutionFileErrorReturn(in_file);
    }
    if (num_row_ok) {
      if (have_row_names) {
        // Use the row name if possible
        return_status = getIndexFromName(log_options, from_method, is_col, name,
                                         lp.row_hash_, iRow, lp.row_names_);
        if (return_status != HighsStatus::kOk) return return_status;
      } else {
        // Have to assume row solution values are in the right order
        iRow = iX;
      }
      read_solution.row_value[iRow] = value;
    }
  }
  if (!num_row_ok) {
    highsLogUser(log_options, HighsLogType::kWarning,
                 "readSolutionFile: Solution file is for %d rows, not %d: row "
                 "values ignored\n",
                 int(num_row), int(lp_num_row));
    // Calculate the row values
    if (calculateRowValuesQuad(lp, read_solution.col_value,
                               read_solution.row_value) != HighsStatus::kOk)
      return readSolutionFileErrorReturn(in_file);
  }
  // OK to have no EOL
  if (!readSolutionFileIgnoreLineOk(in_file))
    return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                  read_solution, read_basis,
                                  in_file);  // EOL
  // OK to have no blank line
  if (!readSolutionFileIgnoreLineOk(in_file))
    return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                  read_solution, read_basis,
                                  in_file);  //
  // OK to have no reference to dual solution values
  if (!readSolutionFileIgnoreLineOk(in_file))
    return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                  read_solution, read_basis,
                                  in_file);  // # Dual solution values
  // If there's a reference to dual solution values, there's a keyword
  // to indicate their status
  if (!readSolutionFileKeywordLineOk(keyword, in_file))
    return readSolutionFileErrorReturn(in_file);
  if (keyword != "None") {
    if (!readSolutionFileIgnoreLineOk(in_file))
      return readSolutionFileErrorReturn(in_file);  // EOL
    // Next line should be "Columns" and correct number
    if (!readSolutionFileHashKeywordIntLineOk(hash, keyword, value_string,
                                              num_col, in_file)) {
      highsLogUser(log_options, HighsLogType::kError,
                   "readSolutionFile: Error reading line \"%s %s %s\"\n",
                   hash.c_str(), keyword.c_str(), value_string.c_str());
      return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                    read_solution, read_basis, in_file);
    }
    assert(keyword == "Columns");
    double value;
    is_col = true;
    assert(is_col);
    for (HighsInt iX = 0; iX < num_col; iX++) {
      if (!readSolutionFileIdDoubleLineOk(name, value, in_file)) {
        highsLogUser(
            log_options, HighsLogType::kError,
            "readSolutionFile: Error reading column dual line \"%s %g\"\n",
            name.c_str(), value);
        return readSolutionFileErrorReturn(in_file);
      }
      if (have_col_names) {
        // Use the column name if possible
        return_status = getIndexFromName(log_options, from_method, is_col, name,
                                         lp.col_hash_, iCol, lp.col_names_);
        if (return_status != HighsStatus::kOk) return return_status;
      } else {
        // Have to assume column dual values are in the right order
        iCol = iX;
      }
      read_solution.col_dual[iCol] = value;
    }
    // Read in the col dual values: next line should be "Rows" and correct
    // number
    if (!readSolutionFileHashKeywordIntLineOk(hash, keyword, value_string,
                                              num_col, in_file)) {
      highsLogUser(log_options, HighsLogType::kError,
                   "readSolutionFile: Error reading line \"%s %s %s\"\n",
                   hash.c_str(), keyword.c_str(), value_string.c_str());
      return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                    read_solution, read_basis, in_file);
    }
    assert(keyword == "Rows");
    is_col = false;
    assert(!is_col);
    for (HighsInt iX = 0; iX < num_row; iX++) {
      if (!readSolutionFileIdDoubleLineOk(name, value, in_file)) {
        highsLogUser(
            log_options, HighsLogType::kError,
            "readSolutionFile: Error reading row dual line \"%s %g\"\n",
            name.c_str(), value);
        return readSolutionFileErrorReturn(in_file);
      }
      if (have_row_names) {
        // Use the row name if possible
        return_status = getIndexFromName(log_options, from_method, is_col, name,
                                         lp.row_hash_, iRow, lp.row_names_);
        if (return_status != HighsStatus::kOk) return return_status;
      } else {
        // Have to assume row dual values are in the right order
        iRow = iX;
      }
      read_solution.row_dual[iRow] = value;
    }
  }
  // OK to have no EOL
  if (!readSolutionFileIgnoreLineOk(in_file))
    return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                  read_solution, read_basis,
                                  in_file);  // EOL
  // OK to have no blank line
  if (!readSolutionFileIgnoreLineOk(in_file))
    return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                  read_solution, read_basis,
                                  in_file);  //
  // OK to have no reference to basis
  if (!readSolutionFileIgnoreLineOk(in_file))
    return readSolutionFileReturn(HighsStatus::kOk, solution, basis,
                                  read_solution, read_basis,
                                  in_file);  // # Basis
  HighsStatus basis_read_status =
      readBasisStream(log_options, lp, read_basis, in_file);
  // Return with basis read status
  return readSolutionFileReturn(basis_read_status, solution, basis,
                                read_solution, read_basis, in_file);
}

HighsStatus readSolutionFileErrorReturn(std::ifstream& in_file) {
  in_file.close();
  return HighsStatus::kError;
}

HighsStatus readSolutionFileReturn(const HighsStatus status,
                                   HighsSolution& solution, HighsBasis& basis,
                                   const HighsSolution& read_solution,
                                   const HighsBasis& read_basis,
                                   std::ifstream& in_file) {
  in_file.close();
  if (status != HighsStatus::kOk) {
    return status;
  }
  solution = read_solution;
  basis = read_basis;
  return HighsStatus::kOk;
}

bool readSolutionFileIgnoreLineOk(std::ifstream& in_file) {
  if (in_file.eof()) return false;
  in_file.ignore(kMaxLineLength, '\n');
  return true;
}

bool readSolutionFileKeywordLineOk(std::string& keyword,
                                   std::ifstream& in_file) {
  if (in_file.eof()) return false;
  in_file >> keyword;
  return true;
}

bool readSolutionFileHashKeywordIntLineOk(std::string& hash,
                                          std::string& keyword,
                                          std::string& value_string,
                                          HighsInt& value,
                                          std::ifstream& in_file) {
  hash = "";
  keyword = "";
  value_string = "";
  // Read the hash symbol
  if (in_file.eof()) return false;
  in_file >> hash;  // #
  if (hash != "#") return false;

  // Read the keyword
  if (in_file.eof()) return false;
  in_file >> keyword;  // keyword

  // Read the value
  if (in_file.eof()) return false;
  // Read as a string, and then check it only contains digits
  in_file >> value_string;
  if (value_string[std::strspn(value_string.c_str(), kLegalIntegerCharacters)])
    return false;
  value = std::stoi(value_string);  // integer value
  return true;
}

bool readSolutionFileColumnHeaderLineOk(
    std::string& hash, std::string& keyword, std::string& value_string,
    HighsInt& value, std::string& qualifier_string,
    std::stringstream& column_section_line_ss) {
  hash = "";
  keyword = "";
  value_string = "";
  qualifier_string = "";
  // Read the hash symbol
  if (column_section_line_ss.eof()) return false;
  column_section_line_ss >> hash;  // #
  if (hash != "#") return false;

  // Read the keyword
  if (column_section_line_ss.eof()) return false;
  column_section_line_ss >> keyword;  // keyword

  // Read the value
  if (column_section_line_ss.eof()) return false;
  // Read as a string, and then check it only contains digits
  column_section_line_ss >> value_string;
  if (value_string[std::strspn(value_string.c_str(), kLegalIntegerCharacters)])
    return false;
  value = std::stoi(value_string);  // integer value
  // See whether there is a qualifier string
  if (!column_section_line_ss.eof()) {
    column_section_line_ss >> qualifier_string;
    tolower(qualifier_string);
    if (qualifier_string != kHighsSparseString &&
        qualifier_string != kHighsPartialString)
      qualifier_string = "";
  }

  return true;
}

bool readSolutionFileColumnLineOk(std::string& id, double& value,
                                  HighsInt& index,
                                  std::stringstream& column_section_line_ss) {
  if (column_section_line_ss.eof()) return false;
  column_section_line_ss >> id;  // Id
  if (column_section_line_ss.eof()) return false;
  column_section_line_ss >> value;  // double value
  // See whether there is an index
  index = -1;
  if (!column_section_line_ss.eof()) {
    std::string index_string;
    // Read as a string, and then check it only contains digits
    column_section_line_ss >> index_string;
    if (index_string[std::strspn(index_string.c_str(),
                                 kLegalIntegerCharacters)])
      return false;
    index = std::stoi(index_string);
  }
  return true;
}

bool readSolutionFileIdIgnoreLineOk(std::string& id, std::ifstream& in_file) {
  if (in_file.eof()) return false;
  in_file >> id;  // Id
  in_file.ignore(kMaxLineLength, '\n');
  return true;
}

bool readSolutionFileIdDoubleLineOk(std::string& id, double& value,
                                    std::ifstream& in_file) {
  if (in_file.eof()) return false;
  in_file >> id;  // Id
  if (in_file.eof()) return false;
  in_file >> value;  // double value
  return true;
}

void assessColPrimalSolution(const bool assess_as_mip,
                             const HighsOptions& options, const double primal,
                             const double lower, const double upper,
                             const HighsVarType type, double& col_infeasibility,
                             double& integer_infeasibility) {
  const double kPrimalFeasibilityTolerance =
      assess_as_mip ? options.mip_feasibility_tolerance
                    : options.primal_feasibility_tolerance;
  // @primal_infeasibility calculation
  col_infeasibility = 0;
  if (primal < lower - kPrimalFeasibilityTolerance) {
    col_infeasibility = lower - primal;
  } else if (primal > upper + kPrimalFeasibilityTolerance) {
    col_infeasibility = primal - upper;
  }
  integer_infeasibility = 0;
  if (assess_as_mip) {
    if (type == HighsVarType::kInteger || type == HighsVarType::kSemiInteger) {
      integer_infeasibility = fractionality(primal);
      assert(integer_infeasibility >= 0);
      if (integer_infeasibility <= options.mip_feasibility_tolerance)
        integer_infeasibility = 0;
    }
    if (col_infeasibility > 0 && (type == HighsVarType::kSemiContinuous ||
                                  type == HighsVarType::kSemiInteger)) {
      // Semi-variables at zero will have positive col infeasibility,
      // so possibly zero this
      if (std::fabs(primal) <= options.mip_feasibility_tolerance)
        col_infeasibility = 0;
      // If there is (still) column infeasibility due to value being
      // off zero or below lower bound, then this counts as an integer
      // infeasibility
      if (col_infeasibility && primal < upper)
        integer_infeasibility =
            std::max(col_infeasibility, integer_infeasibility);
    }
  }
}

// Determine validity, primal feasibility and (when relevant) integer
// feasibility of a solution
HighsStatus assessLpPrimalSolution(const std::string& message,
                                   const HighsOptions& options,
                                   const HighsLp& lp,
                                   const HighsSolution& solution, bool& valid,
                                   bool& integral, bool& feasible) {
  valid = false;
  integral = false;
  feasible = false;
  HighsInt num_col_infeasibilities = 0;
  double max_col_infeasibility = 0;
  double sum_col_infeasibilities = 0;
  HighsInt num_integer_infeasibilities = 0;
  double max_integer_infeasibility = 0;
  double sum_integer_infeasibilities = 0;
  HighsInt num_row_infeasibilities = 0;
  double max_row_infeasibility = 0;
  double sum_row_infeasibilities = 0;
  HighsInt num_row_residuals = 0;
  double max_row_residual = 0;
  double sum_row_residuals = 0;
  const double kRowResidualTolerance =
      options.primal_feasibility_tolerance;  // 1e-12;
  const bool assess_as_mip = lp.isMip() && !options.solve_relaxation;
  const double kPrimalFeasibilityTolerance =
      assess_as_mip ? options.mip_feasibility_tolerance
                    : options.primal_feasibility_tolerance;
  // Warning logging is done if an infeasibility or row residual
  // exceeding the tolerance is found.  Ultimately, if no warnings
  // have been issued, warning is set to false, and a one-line OK
  // logging message is issued.
  bool warning = true;
  bool logged_header = false;

  auto logHeader = [&]() {
    if (logged_header) return;
    highsLogUser(
        options.log_options,
        warning ? HighsLogType::kWarning : HighsLogType::kInfo,
        "%sAssessing feasibility of %s of %.4g%s\n", message.c_str(),
        assess_as_mip
            ? "MIP solution using primal feasibility and integrality tolerances"
            : "LP solution using primal feasibility tolerance",
        kPrimalFeasibilityTolerance, feasible ? ": feasible" : "");
    logged_header = true;
  };

  vector<double> row_value;
  row_value.assign(lp.num_row_, 0);
  const bool have_integrality = (lp.integrality_.size() != 0);
  if (!solution.value_valid) return HighsStatus::kError;
  const bool have_col_names =
      lp.col_names_.size() == static_cast<size_t>(lp.num_col_);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const std::string name_string =
        have_col_names ? (" (" + lp.col_names_[iCol] + ")") : "";
    const double primal = solution.col_value[iCol];
    const double lower = lp.col_lower_[iCol];
    const double upper = lp.col_upper_[iCol];
    const HighsVarType type =
        have_integrality ? lp.integrality_[iCol] : HighsVarType::kContinuous;

    double col_infeasibility = 0;
    double integer_infeasibility = 0;

    assessColPrimalSolution(assess_as_mip, options, primal, lower, upper, type,
                            col_infeasibility, integer_infeasibility);
    if (col_infeasibility > 0) {
      if (col_infeasibility > kPrimalFeasibilityTolerance) {
        if (col_infeasibility > 2 * max_col_infeasibility) {
          logHeader();
          highsLogUser(options.log_options, HighsLogType::kWarning,
                       "Col %6d%s has         infeasibility of %11.4g from "
                       "[lower, value, upper] = [%15.8g; %15.8g; %15.8g]\n",
                       int(iCol), name_string.c_str(), col_infeasibility, lower,
                       primal, upper);
        }
        num_col_infeasibilities++;
      }
      max_col_infeasibility =
          std::max(col_infeasibility, max_col_infeasibility);
      sum_col_infeasibilities += col_infeasibility;
    }
    if (integer_infeasibility > 0) {
      if (integer_infeasibility > options.mip_feasibility_tolerance) {
        if (integer_infeasibility > 2 * max_integer_infeasibility) {
          logHeader();
          highsLogUser(options.log_options, HighsLogType::kWarning,
                       "Col %6d%s has integer infeasibility of %11.4g\n",
                       (int)iCol, name_string.c_str(), integer_infeasibility);
        }
        num_integer_infeasibilities++;
      }
      max_integer_infeasibility =
          std::max(integer_infeasibility, max_integer_infeasibility);
      sum_integer_infeasibilities += integer_infeasibility;
    }
  }
  HighsStatus return_status =
      calculateRowValuesQuad(lp, solution.col_value, row_value);
  if (return_status != HighsStatus::kOk) return return_status;
  const bool have_row_names =
      lp.row_names_.size() >= static_cast<size_t>(lp.num_row_);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    const double primal = solution.row_value[iRow];
    const double lower = lp.row_lower_[iRow];
    const double upper = lp.row_upper_[iRow];
    const std::string name_string =
        have_row_names ? (" (" + lp.row_names_[iRow] + ")") : "";
    // @primal_infeasibility calculation
    double row_infeasibility = 0;
    if (primal < lower - kPrimalFeasibilityTolerance) {
      row_infeasibility = lower - primal;
    } else if (primal > upper + kPrimalFeasibilityTolerance) {
      row_infeasibility = primal - upper;
    }
    if (row_infeasibility > 0) {
      if (row_infeasibility > kPrimalFeasibilityTolerance) {
        if (row_infeasibility > 2 * max_row_infeasibility) {
          logHeader();
          highsLogUser(options.log_options, HighsLogType::kWarning,
                       "Row %6d%s has         infeasibility of %11.4g from "
                       "[lower, value, upper] = [%15.8g; %15.8g; %15.8g]\n",
                       int(iRow), name_string.c_str(), row_infeasibility, lower,
                       primal, upper);
        }
        num_row_infeasibilities++;
      }
      max_row_infeasibility =
          std::max(row_infeasibility, max_row_infeasibility);
      sum_row_infeasibilities += row_infeasibility;
    }
    double row_residual = fabs(primal - row_value[iRow]);
    if (row_residual > kRowResidualTolerance) {
      if (row_residual > 2 * max_row_residual) {
        logHeader();
        highsLogUser(options.log_options, HighsLogType::kWarning,
                     "Row %6d%s has         residual      of %11.4g\n",
                     int(iRow), name_string.c_str(), row_residual);
      }
      num_row_residuals++;
    } else {
      row_residual = 0;
    }
    max_row_residual = std::max(row_residual, max_row_residual);
    sum_row_residuals += row_residual;
  }
  valid = num_row_residuals == 0;
  integral = valid && num_integer_infeasibilities == 0;
  feasible = valid && num_col_infeasibilities == 0 &&
             num_integer_infeasibilities == 0 && num_row_infeasibilities == 0;
  warning = !(integral && feasible);
  if (!warning) {
    // Issue the "OK" logging message
    logHeader();
    return HighsStatus::kOk;
  }
  assert(logged_header);
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "Solution has               num          max          sum\n");
  if (num_col_infeasibilities > 0)
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Col     infeasibilities %6d  %11.4g  %11.4g\n",
                 int(num_col_infeasibilities), max_col_infeasibility,
                 sum_col_infeasibilities);
  if (lp.isMip() && num_integer_infeasibilities > 0)
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Integer infeasibilities %6d  %11.4g  %11.4g\n",
                 int(num_integer_infeasibilities), max_integer_infeasibility,
                 sum_integer_infeasibilities);
  if (num_row_infeasibilities > 0)
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Row     infeasibilities %6d  %11.4g  %11.4g\n",
                 int(num_row_infeasibilities), max_row_infeasibility,
                 sum_row_infeasibilities);
  if (num_row_residuals > 0)
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Row     residuals       %6d  %11.4g  %11.4g\n",
                 int(num_row_residuals), max_row_residual, sum_row_residuals);
  return HighsStatus::kWarning;
}

void writeBasisFile(FILE*& file, const HighsOptions& options, const HighsLp& lp,
                    const HighsBasis& basis) {
  const HighsLogOptions& log_options = options.log_options;
  std::stringstream ss;
  // Basis version line
  ss.str(std::string());
  ss << highsFormatToString("HiGHS_basis_file %s\n", kHighsBasisFileV2.c_str());
  highsFprintfString(file, log_options, ss.str());
  // Basis validity line
  ss.str(std::string());
  if (basis.valid == false) {
    ss << highsFormatToString("None\n");
    highsFprintfString(file, log_options, ss.str());
    return;
  }
  assert(basis.col_status.size() == static_cast<size_t>(lp.num_col_));
  assert(basis.row_status.size() == static_cast<size_t>(lp.num_row_));
  assert(lp.col_names_.size() == static_cast<size_t>(lp.num_col_));
  assert(lp.row_names_.size() == static_cast<size_t>(lp.num_row_));
  ss << highsFormatToString("Valid\n");
  highsFprintfString(file, log_options, ss.str());
  // Column count line
  ss.str(std::string());
  ss << highsFormatToString("# Columns %d\n", int(lp.num_col_));
  highsFprintfString(file, log_options, ss.str());
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const std::string& name = lp.col_names_[iCol];
    assert(name.length() != 0);
    ss.str(std::string());
    ss << highsFormatToString("%s %d\n", name.c_str(),
                              int(basis.col_status[iCol]));
    highsFprintfString(file, log_options, ss.str());
  }
  // Row count line
  ss.str(std::string());
  ss << highsFormatToString("# Rows %d\n", int(lp.num_row_));
  highsFprintfString(file, log_options, ss.str());
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    const std::string& name = lp.row_names_[iRow];
    assert(name.length() != 0);
    ss.str(std::string());
    ss << highsFormatToString("%s %d\n", name.c_str(),
                              int(basis.row_status[iRow]));
    highsFprintfString(file, log_options, ss.str());
  }
}

HighsStatus getIndexFromName(const HighsLogOptions& log_options,
                             std::string& from_method, const bool is_column,
                             const std::string& name, HighsNameHash& name_hash,
                             HighsInt& index,
                             const std::vector<std::string>& names) {
  std::unordered_map<std::string, int>& name2index = name_hash.name2index;
  assert(name2index.size());

  auto nameIsUnique = [&]() {
    auto search = name2index.find(name);
    if (search == name2index.end()) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s: %s name \"%s\" is not found\n", from_method.c_str(),
                   is_column ? "column" : "row", name.c_str());
      return false;
    }
    if (search->second < 0) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s: %s name \"%s\" is duplicated\n", from_method.c_str(),
                   is_column ? "column" : "row", name.c_str());
      return false;
    }
    index = search->second;
    return true;
  };

  if (!nameIsUnique()) return HighsStatus::kError;

  assert(static_cast<size_t>(index) < names.size());

  if (names[index] != name) {
    // If index has been used to store a duplication count that has
    // been reduced to zero, then index will be zero, so is unlikely
    // to be the index of name: name2index must be reformed
    name_hash.form(names);
    const bool name_is_unique = nameIsUnique();
    // If name must be found and not duplicate, since these cases
    // would have been handled after first call to nameIsUnique
    assert(name_is_unique);
    if (!name_is_unique) return HighsStatus::kError;
    const bool index_ok = names[index] == name;
    assert(index_ok);
    if (!index_ok) return HighsStatus::kError;
  }
  return HighsStatus::kOk;
}

HighsStatus readBasisFile(const HighsLogOptions& log_options, HighsLp& lp,
                          HighsBasis& basis, const std::string& filename) {
  // Opens a basis file as an ifstream
  HighsStatus return_status = HighsStatus::kOk;
  std::ifstream in_file;
  in_file.open(filename.c_str(), std::ios::in);
  if (in_file.is_open()) {
    return_status = readBasisStream(log_options, lp, basis, in_file);
    in_file.close();
  } else {
    highsLogUser(log_options, HighsLogType::kError,
                 "readBasisFile: Cannot open readable file \"%s\"\n",
                 filename.c_str());
    return_status = HighsStatus::kError;
  }
  return return_status;
}

HighsStatus readBasisStream(const HighsLogOptions& log_options, HighsLp& lp,
                            HighsBasis& basis, std::ifstream& in_file) {
  // Reads a basis as an ifstream, returning an error if what's read is
  // inconsistent with the sizes of the HighsBasis passed in
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status = HighsStatus::kOk;
  std::string from_method = "readBasisStream";
  std::string string_highs, string_version;
  in_file >> string_highs >> string_version;
  const bool v1 = string_version == kHighsBasisFileV1;
  const bool v2 = string_version == kHighsBasisFileV2;
  const bool have_col_names =
      lp.col_names_.size() == static_cast<size_t>(lp.num_col_);
  const bool have_row_names =
      lp.row_names_.size() == static_cast<size_t>(lp.num_row_);
  if (have_col_names) {
    // Ensure that the col name hash table has been formed
    if (!lp.col_hash_.size()) lp.col_hash_.form(lp.col_names_);
  }
  if (have_row_names) {
    // Ensure that the row name hash table has been formed
    if (!lp.row_hash_.size()) lp.row_hash_.form(lp.row_names_);
  }
  if (v1 || v2) {
    if (v1) {
      // Ability to read v1 basis files is deprecated
      highsLogUser(log_options, HighsLogType::kWarning,
                   "readBasisFile: Basis file format %s is deprecated\n",
                   kHighsBasisFileV1.c_str());
      return_status = HighsStatus::kWarning;
    }
    std::string keyword;
    in_file >> keyword;
    if (keyword == "None") {
      basis.valid = false;
      return return_status;
    }
    const HighsInt basis_num_col = (HighsInt)basis.col_status.size();
    const HighsInt basis_num_row = (HighsInt)basis.row_status.size();
    HighsInt int_status;
    std::string name;
    assert(keyword == "Valid");
    HighsInt num_col, num_row;
    // Read in the columns section
    in_file >> keyword >> keyword;
    assert(keyword == "Columns");
    in_file >> num_col;
    if (num_col != basis_num_col) {
      highsLogUser(log_options, HighsLogType::kError,
                   "readBasisFile: Basis file is for %" HIGHSINT_FORMAT
                   " columns, not %" HIGHSINT_FORMAT "\n",
                   num_col, basis_num_col);
      return HighsStatus::kError;
    }
    bool is_column = true;
    if (v1) {
      for (HighsInt iCol = 0; iCol < num_col; iCol++) {
        in_file >> int_status;
        basis.col_status[iCol] = (HighsBasisStatus)int_status;
      }
    } else {
      HighsInt iCol;
      for (HighsInt iX = 0; iX < num_col; iX++) {
        in_file >> name >> int_status;
        if (have_col_names) {
          // Use the column name if possible
          call_status =
              getIndexFromName(log_options, from_method, is_column, name,
                               lp.col_hash_, iCol, lp.col_names_);
          if (call_status != HighsStatus::kOk) return call_status;
        } else {
          // Have to assume column basis status are in the right order
          iCol = iX;
        }
        basis.col_status[iCol] = HighsBasisStatus(int_status);
      }
    }
    // Read in the rows section
    in_file >> keyword >> keyword;
    assert(keyword == "Rows");
    in_file >> num_row;
    if (num_row != basis_num_row) {
      highsLogUser(log_options, HighsLogType::kError,
                   "readBasisFile: Basis file is for %" HIGHSINT_FORMAT
                   " rows, not %" HIGHSINT_FORMAT "\n",
                   num_row, basis_num_row);
      return HighsStatus::kError;
    }
    is_column = false;
    if (v1) {
      for (HighsInt iRow = 0; iRow < num_row; iRow++) {
        in_file >> int_status;
        basis.row_status[iRow] = (HighsBasisStatus)int_status;
      }
    } else {
      HighsInt iRow;
      for (HighsInt iX = 0; iX < num_row; iX++) {
        in_file >> name >> int_status;
        if (have_row_names) {
          // Use the row name if possible
          call_status =
              getIndexFromName(log_options, from_method, is_column, name,
                               lp.row_hash_, iRow, lp.row_names_);
          if (call_status != HighsStatus::kOk) return call_status;
        } else {
          // Have to assume row solution values are in the right order
          iRow = iX;
        }
        basis.row_status[iRow] = HighsBasisStatus(int_status);
      }
    }
  } else {
    highsLogUser(log_options, HighsLogType::kError,
                 "readBasisFile: Cannot read basis file for HiGHS %s\n",
                 string_version.c_str());
    return_status = HighsStatus::kError;
  }
  return return_status;
}

HighsStatus calculateColDualsQuad(const HighsLp& lp, HighsSolution& solution) {
  const bool correct_size = int(solution.row_dual.size()) == lp.num_row_;
  const bool is_colwise = lp.a_matrix_.isColwise();
  const bool data_error = !correct_size || !is_colwise;
  assert(!data_error);
  if (data_error) return HighsStatus::kError;

  std::vector<HighsCDouble> col_dual_quad;
  col_dual_quad.assign(lp.num_col_, HighsCDouble{0.0});

  for (HighsInt col = 0; col < lp.num_col_; col++) {
    for (HighsInt i = lp.a_matrix_.start_[col];
         i < lp.a_matrix_.start_[col + 1]; i++) {
      const HighsInt row = lp.a_matrix_.index_[i];
      assert(row >= 0);
      assert(row < lp.num_row_);
      col_dual_quad[col] += static_cast<HighsCDouble>(solution.row_dual[row]) *
                            lp.a_matrix_.value_[i];
    }
    col_dual_quad[col] += lp.col_cost_[col];
  }

  // assign quad values to double vector
  solution.col_dual.resize(lp.num_col_);
  std::transform(col_dual_quad.begin(), col_dual_quad.end(),
                 solution.col_dual.begin(),
                 [](HighsCDouble x) { return double(x); });

  return HighsStatus::kOk;
}

HighsStatus calculateRowValuesQuad(const HighsLp& lp,
                                   const std::vector<double>& col_value,
                                   std::vector<double>& row_value,
                                   const HighsInt report_row) {
  const bool correct_size = int(col_value.size()) == lp.num_col_;
  const bool is_colwise = lp.a_matrix_.isColwise();
  const bool data_error = !correct_size || !is_colwise;
  assert(!data_error);
  if (data_error) return HighsStatus::kError;

  std::vector<HighsCDouble> row_value_quad;
  row_value_quad.assign(lp.num_row_, HighsCDouble{0.0});

  for (HighsInt col = 0; col < lp.num_col_; col++) {
    for (HighsInt i = lp.a_matrix_.start_[col];
         i < lp.a_matrix_.start_[col + 1]; i++) {
      const HighsInt row = lp.a_matrix_.index_[i];
      assert(row >= 0);
      assert(row < lp.num_row_);
      row_value_quad[row] +=
          static_cast<HighsCDouble>(col_value[col]) * lp.a_matrix_.value_[i];
      if (row == report_row) {
        printf(
            "calculateRowValuesQuad: Row %d becomes %g due to contribution of "
            ".col_value[%d] = %g\n",
            int(row), double(row_value_quad[row]), int(col), col_value[col]);
      }
    }
  }

  // assign quad values to double vector
  row_value.resize(lp.num_row_);
  std::transform(row_value_quad.begin(), row_value_quad.end(),
                 row_value.begin(), [](HighsCDouble x) { return double(x); });

  return HighsStatus::kOk;
}

HighsStatus calculateRowValuesQuad(const HighsLp& lp, HighsSolution& solution,
                                   const HighsInt report_row) {
  return calculateRowValuesQuad(lp, solution.col_value, solution.row_value,
                                report_row);
}

bool isColDataNull(const HighsLogOptions& log_options,
                   const double* usr_col_cost, const double* usr_col_lower,
                   const double* usr_col_upper) {
  bool null_data = false;
  null_data =
      doubleUserDataNotNull(log_options, usr_col_cost, "column costs") ||
      null_data;
  null_data = doubleUserDataNotNull(log_options, usr_col_lower,
                                    "column lower bounds") ||
              null_data;
  null_data = doubleUserDataNotNull(log_options, usr_col_upper,
                                    "column upper bounds") ||
              null_data;
  return null_data;
}

bool isRowDataNull(const HighsLogOptions& log_options,
                   const double* usr_row_lower, const double* usr_row_upper) {
  bool null_data = false;
  null_data =
      doubleUserDataNotNull(log_options, usr_row_lower, "row lower bounds") ||
      null_data;
  null_data =
      doubleUserDataNotNull(log_options, usr_row_upper, "row upper bounds") ||
      null_data;
  return null_data;
}

bool isMatrixDataNull(const HighsLogOptions& log_options,
                      const HighsInt* usr_matrix_start,
                      const HighsInt* usr_matrix_index,
                      const double* usr_matrix_value) {
  bool null_data = false;
  null_data =
      intUserDataNotNull(log_options, usr_matrix_start, "matrix starts") ||
      null_data;
  null_data =
      intUserDataNotNull(log_options, usr_matrix_index, "matrix indices") ||
      null_data;
  null_data =
      doubleUserDataNotNull(log_options, usr_matrix_value, "matrix values") ||
      null_data;
  return null_data;
}

void reportPresolveReductions(const HighsLogOptions& log_options,
                              HighsPresolveStatus presolve_status,
                              const HighsLp& lp, const HighsLp& presolved_lp) {
  const HighsInt num_col_from = lp.num_col_;
  const HighsInt num_row_from = lp.num_row_;
  const HighsInt num_nz_from = lp.numNz();
  HighsInt num_col_to = 0;
  HighsInt num_row_to = 0;
  HighsInt num_nz_to = 0;
  std::string message = "";

  switch (presolve_status) {
    case HighsPresolveStatus::kNotPresolved:
    case HighsPresolveStatus::kInfeasible:
    case HighsPresolveStatus::kUnboundedOrInfeasible:
      return;
    case HighsPresolveStatus::kNotReduced: {
      num_col_to = num_col_from;
      num_row_to = num_row_from;
      num_nz_to = num_nz_from;
      message = "- Not reduced";
      break;
    }
    case HighsPresolveStatus::kReduced:
    case HighsPresolveStatus::kTimeout: {
      num_col_to = presolved_lp.num_col_;
      num_row_to = presolved_lp.num_row_;
      num_nz_to = presolved_lp.numNz();
      message =
          presolve_status == HighsPresolveStatus::kTimeout ? "- Timeout" : "";
      break;
    }
    case HighsPresolveStatus::kReducedToEmpty: {
      num_col_to = 0;
      num_row_to = 0;
      num_nz_to = 0;
      message = "- Reduced to empty";
      break;
    }
    default: {
      assert(presolve_status == HighsPresolveStatus::kOutOfMemory ||
             presolve_status == HighsPresolveStatus::kException ||
             presolve_status == HighsPresolveStatus::kMatrixValueError);
      return;
    }
  }
  char nz_sign_char = '-';
  HighsInt delta_nz = num_nz_from - num_nz_to;
  if (num_nz_to > num_nz_from) {
    delta_nz = -delta_nz;
    nz_sign_char = '+';
  }
  highsLogUser(log_options, HighsLogType::kInfo,
               "Presolve reductions: rows %" HIGHSINT_FORMAT
               "(-%" HIGHSINT_FORMAT "); columns %" HIGHSINT_FORMAT
               "(-%" HIGHSINT_FORMAT
               "); "
               "nonzeros %" HIGHSINT_FORMAT "(%c%" HIGHSINT_FORMAT ") %s\n",
               num_row_to, (num_row_from - num_row_to), num_col_to,
               (num_col_from - num_col_to), num_nz_to, nz_sign_char, delta_nz,
               message.c_str());
}

bool isLessInfeasibleDSECandidate(const HighsLogOptions& log_options,
                                  const HighsLp& lp) {
  HighsInt max_col_num_en = -1;
  const HighsInt max_allowed_col_num_en = 24;
  const HighsInt max_assess_col_num_en =
      std::max(HighsInt{9}, max_allowed_col_num_en);
  const HighsInt max_average_col_num_en = 6;
  vector<HighsInt> col_length_k;
  col_length_k.resize(1 + max_assess_col_num_en, 0);
  bool LiDSE_candidate = true;
  for (HighsInt col = 0; col < lp.num_col_; col++) {
    // Check limit on number of entries in the column has not been breached
    HighsInt col_num_en =
        lp.a_matrix_.start_[col + 1] - lp.a_matrix_.start_[col];
    max_col_num_en = std::max(col_num_en, max_col_num_en);
    if (col_num_en > max_assess_col_num_en) return false;
    col_length_k[col_num_en]++;
    for (HighsInt en = lp.a_matrix_.start_[col];
         en < lp.a_matrix_.start_[col + 1]; en++) {
      double value = lp.a_matrix_.value_[en];
      // All nonzeros must be +1 or -1
      if (fabs(value) != 1) return false;
    }
  }
  double average_col_num_en = lp.a_matrix_.start_[lp.num_col_];
  average_col_num_en = average_col_num_en / lp.num_col_;
  LiDSE_candidate =
      LiDSE_candidate && average_col_num_en <= max_average_col_num_en;
  highsLogDev(log_options, HighsLogType::kInfo,
              "LP %s has all |entries|=1; max column count = %" HIGHSINT_FORMAT
              " (limit %" HIGHSINT_FORMAT
              "); average "
              "column count = %0.2g (limit %" HIGHSINT_FORMAT
              "): LP %s a candidate for LiDSE\n",
              lp.model_name_.c_str(), max_col_num_en, max_allowed_col_num_en,
              average_col_num_en, max_average_col_num_en,
              LiDSE_candidate ? "is" : "is not");
  return LiDSE_candidate;
}

HighsLp withoutSemiVariables(const HighsLp& lp_, HighsSolution& solution,
                             const double mip_feasibility_tolerance) {
  HighsLp lp = lp_;
  HighsInt num_col = lp.num_col_;
  HighsInt num_row = lp.num_row_;
  HighsInt num_semi_variables = 0;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
        lp.integrality_[iCol] == HighsVarType::kSemiInteger)
      num_semi_variables++;
  }
  assert(num_semi_variables);
  // Insert spaces for index/value of new coefficients for
  // semi-variables
  vector<HighsInt>& start = lp.a_matrix_.start_;
  vector<HighsInt>& index = lp.a_matrix_.index_;
  vector<double>& value = lp.a_matrix_.value_;
  HighsInt num_nz = start[num_col];
  HighsInt new_num_nz = num_nz + 2 * num_semi_variables;
  HighsInt new_el = new_num_nz;
  index.resize(new_num_nz);
  value.resize(new_num_nz);
  for (HighsInt iCol = num_col - 1; iCol >= 0; iCol--) {
    HighsInt from_el = start[iCol + 1] - 1;
    start[iCol + 1] = new_el;
    if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
        lp.integrality_[iCol] == HighsVarType::kSemiInteger)
      new_el -= 2;
    for (HighsInt iEl = from_el; iEl >= start[iCol]; iEl--) {
      new_el--;
      index[new_el] = index[iEl];
      value[new_el] = value[iEl];
    }
  }
  assert(new_el == 0);
  // Insert the new coefficients for semi-variables
  HighsInt row_num = num_row;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
        lp.integrality_[iCol] == HighsVarType::kSemiInteger) {
      HighsInt iEl = start[iCol + 1] - 2;
      index[iEl] = row_num++;
      value[iEl] = 1;
      iEl++;
      index[iEl] = row_num++;
      value[iEl] = 1;
    }
  }
  num_nz = start[num_col];
  new_num_nz = num_nz + 2 * num_semi_variables;
  row_num = num_row;
  HighsInt semi_col_num = 0;
  HighsInt semi_row_num = 0;
  // Insert the new variables and their coefficients
  std::stringstream ss;
  const bool have_col_names =
      lp.col_names_.size() == static_cast<size_t>(lp.num_col_);
  const bool have_row_names =
      lp.row_names_.size() == static_cast<size_t>(lp.num_row_);
  const bool have_solution = solution.value_valid;
  if (have_solution) {
    // Create zeroed row values for the new rows
    assert((HighsInt)solution.row_value.size() == lp_.num_row_);
    for (HighsInt iCol = 0; iCol < 2 * num_semi_variables; iCol++)
      solution.row_value.push_back(0);
    assert((HighsInt)solution.col_value.size() == lp_.num_col_);
    assert((HighsInt)solution.row_value.size() ==
           lp_.num_row_ + 2 * num_semi_variables);
  }
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
        lp.integrality_[iCol] == HighsVarType::kSemiInteger) {
      const double semi_lower_bound = lp.col_lower_[iCol];
      const double semi_upper_bound = lp.col_upper_[iCol];
      // Add a binary variable with zero cost
      lp.col_cost_.push_back(0);
      lp.col_lower_.push_back(0);
      lp.col_upper_.push_back(1);
      // Complete x - l*y >= 0
      lp.row_lower_.push_back(0);
      lp.row_upper_.push_back(kHighsInf);
      if (have_col_names) {
        // Create a column name
        ss.str(std::string());
        ss << "semi_binary_" << semi_col_num++;
        lp.col_names_.push_back(ss.str());
      }
      if (have_row_names) {
        // Create a row name
        ss.str(std::string());
        ss << "semi_lb_" << semi_row_num;
        lp.row_names_.push_back(ss.str());
      }
      index.push_back(row_num++);
      value.push_back(-semi_lower_bound);
      // Accommodate any primal solution
      if (have_solution) {
        // Record the previous solution value so any change can be
        // determined
        const double prev_primal = solution.col_value[iCol];
        if (solution.col_value[iCol] <= mip_feasibility_tolerance) {
          // Currently at or below zero, so binary is 0
          solution.col_value[iCol] = 0;
          solution.col_value.push_back(0);
        } else {
          // Otherwise, solution is at least lower bound, and binary
          // is 1
          solution.col_value[iCol] =
              std::max(semi_lower_bound, solution.col_value[iCol]);
          solution.col_value.push_back(1);
        }
        const double dl_primal = solution.col_value[iCol] - prev_primal;
        if (dl_primal) {
          // Change in primal value, so update row values. NB start
          // has been extended to incorporate the values in this
          // column for the new rows. Their solution values are zero,
          // but will be set later
          for (HighsInt iEl = start[iCol]; iEl < start[iCol + 1]; iEl++)
            solution.row_value[index[iEl]] += dl_primal * value[iEl];
        }
        const HighsInt new_col = lp.col_cost_.size() - 1;
        const double binary_value = solution.col_value[new_col];
        solution.row_value[row_num - 1] =
            solution.col_value[iCol] - semi_lower_bound * binary_value;
        solution.row_value[row_num] =
            solution.col_value[iCol] - semi_upper_bound * binary_value;
      }
      // Complete x - u*y <= 0
      lp.row_lower_.push_back(-kHighsInf);
      lp.row_upper_.push_back(0);
      if (have_row_names) {
        // Create a row name
        ss.str(std::string());
        ss << "semi_ub_" << semi_row_num++;
        lp.row_names_.push_back(ss.str());
      }
      index.push_back(row_num++);
      value.push_back(-semi_upper_bound);
      // Add the next start
      start.push_back(index.size());
      lp.integrality_.push_back(HighsVarType::kInteger);
      if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous) {
        lp.integrality_[iCol] = HighsVarType::kContinuous;
      } else if (lp.integrality_[iCol] == HighsVarType::kSemiInteger) {
        lp.integrality_[iCol] = HighsVarType::kInteger;
      }
      // Change the lower bound on the semi-variable to zero. Cannot
      // do this earlier, as its original value is used in constraint
      // 0 <= x-l*y
      lp.col_lower_[iCol] = 0;
    }
  }
  num_col += num_semi_variables;
  lp.num_col_ += num_semi_variables;
  lp.num_row_ += 2 * num_semi_variables;
  assert((HighsInt)index.size() == new_num_nz);
  // Ensure that the matrix dimensions are consistent with the LP
  // dimensions
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  // Clear any modifications inherited from lp_
  lp.mods_.clear();
  return lp;
}

void removeRowsOfCountOne(const HighsLogOptions& log_options, HighsLp& lp) {
  vector<HighsInt>& a_start = lp.a_matrix_.start_;
  vector<HighsInt>& a_index = lp.a_matrix_.index_;
  vector<double>& a_value = lp.a_matrix_.value_;
  vector<HighsInt> a_count;
  vector<HighsInt> ar_count;
  vector<HighsInt> ar_start;
  vector<HighsInt> ar_index;
  vector<double> ar_value;
  const bool have_name = lp.row_names_.size() > 0;
  HighsInt num_nz = a_start[lp.num_col_];
  const HighsInt original_num_nz = num_nz;
  const HighsInt original_num_row = lp.num_row_;
  HighsInt num_row_count_1 = 0;
  ar_count.assign(lp.num_row_, 0);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    for (HighsInt iEl = a_start[iCol]; iEl < a_start[iCol + 1]; iEl++)
      ar_count[a_index[iEl]]++;
  }
  ar_start.push_back(0);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    ar_start.push_back(ar_start[iRow] + ar_count[iRow]);
    ar_count[iRow] = ar_start[iRow];
  }
  ar_index.resize(num_nz);
  ar_value.resize(num_nz);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    for (HighsInt iEl = a_start[iCol]; iEl < a_start[iCol + 1]; iEl++) {
      HighsInt iRow = a_index[iEl];
      ar_index[ar_count[iRow]] = iCol;
      ar_value[ar_count[iRow]] = a_value[iEl];
      ar_count[iRow]++;
    }
  }
  HighsInt newRow = 0;
  HighsInt newEl = 0;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    const HighsInt row_count = ar_start[iRow + 1] - ar_start[iRow];
    if (row_count == 1) {
      HighsInt iCol = ar_index[ar_start[iRow]];
      double value = ar_value[ar_start[iRow]];
      assert(value);
      if (value > 0) {
        if (lp.row_lower_[iRow] > -kHighsInf)
          lp.col_lower_[iCol] =
              std::max(lp.row_lower_[iRow] / value, lp.col_lower_[iCol]);
        if (lp.row_upper_[iRow] < kHighsInf)
          lp.col_upper_[iCol] =
              std::min(lp.row_upper_[iRow] / value, lp.col_upper_[iCol]);
      } else {
        if (lp.row_lower_[iRow] > -kHighsInf)
          lp.col_upper_[iCol] =
              std::min(lp.row_lower_[iRow] / value, lp.col_upper_[iCol]);
        if (lp.row_upper_[iRow] < kHighsInf)
          lp.col_lower_[iCol] =
              std::max(lp.row_upper_[iRow] / value, lp.col_lower_[iCol]);
      }
      num_row_count_1++;
      continue;
    }
    lp.row_lower_[newRow] = lp.row_lower_[iRow];
    lp.row_upper_[newRow] = lp.row_upper_[iRow];
    if (have_name) lp.row_names_[newRow] = lp.row_names_[iRow];
    ar_start[newRow] = newEl;
    for (HighsInt iEl = ar_start[iRow]; iEl < ar_start[iRow + 1]; iEl++) {
      ar_index[newEl] = ar_index[iEl];
      ar_value[newEl] = ar_value[iEl];
      newEl++;
    }
    newRow++;
  }
  ar_start[newRow] = newEl;
  lp.num_row_ = newRow;
  lp.row_lower_.resize(newRow);
  lp.row_upper_.resize(newRow);
  if (have_name) lp.row_names_.resize(newRow);

  num_nz = ar_start[lp.num_row_];
  a_count.assign(lp.num_col_, 0);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    for (HighsInt iEl = ar_start[iRow]; iEl < ar_start[iRow + 1]; iEl++)
      a_count[ar_index[iEl]]++;
  }
  a_start[0] = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    a_start[iCol + 1] = a_start[iCol] + a_count[iCol];
    a_count[iCol] = a_start[iCol];
  }
  a_index.resize(num_nz);
  a_value.resize(num_nz);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    for (HighsInt iEl = ar_start[iRow]; iEl < ar_start[iRow + 1]; iEl++) {
      HighsInt iCol = ar_index[iEl];
      a_index[a_count[iCol]] = iRow;
      a_value[a_count[iCol]] = ar_value[iEl];
      a_count[iCol]++;
    }
  }
  assert(original_num_row - lp.num_row_ == num_row_count_1);
  assert(original_num_nz - num_nz == num_row_count_1);
  highsLogUser(log_options, HighsLogType::kWarning,
               "Removed %d rows of count 1\n", int(num_row_count_1));
}

void getSubVectors(const HighsIndexCollection& index_collection,
                   const HighsInt data_dim, const double* data0,
                   const double* data1, const double* data2,
                   const HighsSparseMatrix& matrix, HighsInt& num_sub_vector,
                   double* sub_vector_data0, double* sub_vector_data1,
                   double* sub_vector_data2, HighsInt& sub_matrix_num_nz,
                   HighsInt* sub_matrix_start, HighsInt* sub_matrix_index,
                   double* sub_matrix_value) {
  // Ensure that if there's no data0 then it's not required in the
  // sub-vector
  if (data0 == nullptr) assert(sub_vector_data0 == nullptr);
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  // Surely this is checked elsewhere
  assert(0 <= from_k && to_k < data_dim);
  assert(from_k <= to_k);
  HighsInt out_from_vector;
  HighsInt out_to_vector;
  HighsInt in_from_vector;
  HighsInt in_to_vector = -1;
  HighsInt current_set_entry = 0;
  num_sub_vector = 0;
  sub_matrix_num_nz = 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, out_from_vector, out_to_vector,
                     in_from_vector, in_to_vector, current_set_entry);
    assert(out_to_vector < data_dim);
    assert(in_to_vector < data_dim);
    for (HighsInt iVector = out_from_vector; iVector <= out_to_vector;
         iVector++) {
      if (sub_vector_data0 != nullptr)
        sub_vector_data0[num_sub_vector] = data0[iVector];
      if (sub_vector_data1 != nullptr)
        sub_vector_data1[num_sub_vector] = data1[iVector];
      if (sub_vector_data2 != nullptr)
        sub_vector_data2[num_sub_vector] = data2[iVector];
      if (sub_matrix_start != nullptr)
        sub_matrix_start[num_sub_vector] = sub_matrix_num_nz +
                                           matrix.start_[iVector] -
                                           matrix.start_[out_from_vector];
      num_sub_vector++;
    }
    for (HighsInt iEl = matrix.start_[out_from_vector];
         iEl < matrix.start_[out_to_vector + 1]; iEl++) {
      if (sub_matrix_index != nullptr)
        sub_matrix_index[sub_matrix_num_nz] = matrix.index_[iEl];
      if (sub_matrix_value != nullptr)
        sub_matrix_value[sub_matrix_num_nz] = matrix.value_[iEl];
      sub_matrix_num_nz++;
    }
    if (out_to_vector == data_dim - 1 || in_to_vector == data_dim - 1) break;
  }
}

void getSubVectorsTranspose(const HighsIndexCollection& index_collection,
                            const HighsInt data_dim, const double* data0,
                            const double* data1, const double* data2,
                            const HighsSparseMatrix& matrix,
                            HighsInt& num_sub_vector, double* sub_vector_data0,
                            double* sub_vector_data1, double* sub_vector_data2,
                            HighsInt& sub_matrix_num_nz,
                            HighsInt* sub_matrix_start,
                            HighsInt* sub_matrix_index,
                            double* sub_matrix_value) {
  // Ensure that if there's no data0 then it's not required in the
  // sub-vector
  if (data0 == nullptr) assert(sub_vector_data0 == nullptr);
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  // Surely this is checked elsewhere
  assert(0 <= from_k && to_k < data_dim);
  assert(from_k <= to_k);
  // "Out" means not in the set to be extracted
  // "In" means in the set to be extracted
  HighsInt out_from_vector;
  HighsInt out_to_vector;
  HighsInt in_from_vector;
  HighsInt in_to_vector = -1;
  HighsInt current_set_entry = 0;
  // Set up a mask so that entries to be got from the matrix can be
  // identified and have their correct index.
  vector<HighsInt> new_index;
  new_index.resize(data_dim);

  num_sub_vector = 0;
  sub_matrix_num_nz = 0;
  if (!index_collection.is_mask_) {
    out_to_vector = -1;
    current_set_entry = 0;
    for (HighsInt k = from_k; k <= to_k; k++) {
      updateOutInIndex(index_collection, in_from_vector, in_to_vector,
                       out_from_vector, out_to_vector, current_set_entry);
      if (k == from_k) {
        // Account for any initial vectors not being extracted
        for (HighsInt iVector = 0; iVector < in_from_vector; iVector++) {
          new_index[iVector] = -1;
        }
      }
      for (HighsInt iVector = in_from_vector; iVector <= in_to_vector;
           iVector++) {
        new_index[iVector] = num_sub_vector;
        num_sub_vector++;
      }
      for (HighsInt iVector = out_from_vector; iVector <= out_to_vector;
           iVector++) {
        new_index[iVector] = -1;
      }
      if (out_to_vector >= data_dim - 1) break;
    }
  } else {
    for (HighsInt iVector = 0; iVector < data_dim; iVector++) {
      if (index_collection.mask_[iVector]) {
        new_index[iVector] = num_sub_vector;
        num_sub_vector++;
      } else {
        new_index[iVector] = -1;
      }
    }
  }

  // Bail out if no vectors are to be extracted
  if (num_sub_vector == 0) return;

  for (HighsInt iVector = 0; iVector < data_dim; iVector++) {
    HighsInt new_iVector = new_index[iVector];
    if (new_iVector >= 0) {
      assert(new_iVector < num_sub_vector);
      if (sub_vector_data0 != NULL)
        sub_vector_data0[new_iVector] = data0[iVector];
      if (sub_vector_data1 != NULL)
        sub_vector_data1[new_iVector] = data1[iVector];
      if (sub_vector_data2 != NULL)
        sub_vector_data2[new_iVector] = data2[iVector];
    }
  }
  const bool extract_start = sub_matrix_start != NULL;
  const bool extract_index = sub_matrix_index != NULL;
  const bool extract_value = sub_matrix_value != NULL;
  const bool extract_matrix = extract_index || extract_value;
  // Allocate an array of lengths for the sub-matrix to be
  // extracted: necessary even if just the number of nonzeros is
  // required
  vector<HighsInt> sub_matrix_length;
  sub_matrix_length.assign(num_sub_vector, 0);
  // Identify the lengths of the vectors in the sub-matrix to be extracted
  HighsInt num_vector = matrix.start_.size() - 1;
  for (HighsInt vector = 0; vector < num_vector; vector++) {
    for (HighsInt iEl = matrix.start_[vector]; iEl < matrix.start_[vector + 1];
         iEl++) {
      HighsInt iVector = matrix.index_[iEl];
      HighsInt new_iVector = new_index[iVector];
      if (new_iVector >= 0) sub_matrix_length[new_iVector]++;
    }
  }
  if (!extract_start) {
    // bail out if no matrix starts are to be extracted, but only after
    // computing the number of nonzeros
    for (HighsInt iVector = 0; iVector < num_sub_vector; iVector++)
      sub_matrix_num_nz += sub_matrix_length[iVector];
    return;
  }
  // Allocate an array of lengths for the sub-matrix to be extracted
  sub_matrix_start[0] = 0;
  for (HighsInt iVector = 0; iVector < num_sub_vector - 1; iVector++) {
    sub_matrix_start[iVector + 1] =
        sub_matrix_start[iVector] + sub_matrix_length[iVector];
    sub_matrix_length[iVector] = sub_matrix_start[iVector];
  }
  HighsInt iVector = num_sub_vector - 1;
  sub_matrix_num_nz = sub_matrix_start[iVector] + sub_matrix_length[iVector];
  // Bail out if matrix indices and values are not required
  if (!extract_matrix) return;
  sub_matrix_length[iVector] = sub_matrix_start[iVector];
  // Fill the row-wise matrix with indices and values
  for (HighsInt vector = 0; vector < num_vector; vector++) {
    for (HighsInt iEl = matrix.start_[vector]; iEl < matrix.start_[vector + 1];
         iEl++) {
      HighsInt iVector = matrix.index_[iEl];
      HighsInt new_iVector = new_index[iVector];
      if (new_iVector >= 0) {
        HighsInt row_iEl = sub_matrix_length[new_iVector];
        if (extract_index) sub_matrix_index[row_iEl] = vector;
        if (extract_value) sub_matrix_value[row_iEl] = matrix.value_[iEl];
        sub_matrix_length[new_iVector]++;
      }
    }
  }
}

std::string highsVarTypeToString(const HighsVarType type) {
  switch (type) {
    case HighsVarType::kContinuous:
      return "continuous";
    case HighsVarType::kInteger:
      return "integer";
    case HighsVarType::kSemiContinuous:
      return "semi continuous";
    case HighsVarType::kSemiInteger:
      return "semi integer";
    case HighsVarType::kImplicitInteger:
      return "implicit integer";
    default:
      return "unknown";
  }
}

std::string highsVarTypeToString(const HighsInt type) {
  if (type < HighsInt(HighsVarType::kContinuous) ||
      type > HighsInt(HighsVarType::kImplicitInteger))
    return "unknown";
  HighsVarType type_ = HighsVarType(uint8_t(type));
  return highsVarTypeToString(type_);
}

void initialiseUserScaleData(const HighsOptions& options,
                             HighsUserScaleData& user_scale_data) {
  user_scale_data.initialise(options.user_objective_scale,
                             options.user_bound_scale, options.infinite_cost,
                             options.infinite_bound, options.small_matrix_value,
                             options.large_matrix_value);
}

void HighsUserScaleData::initialise(const HighsInt& user_objective_scale_,
                                    const HighsInt& user_bound_scale_,
                                    const double& infinite_cost_,
                                    const double& infinite_bound_,
                                    const double& small_matrix_value_,
                                    const double& large_matrix_value_) {
  this->user_objective_scale = user_objective_scale_;
  this->user_bound_scale = user_bound_scale_;
  this->infinite_cost = infinite_cost_;
  this->infinite_bound = infinite_bound_;
  this->small_matrix_value = small_matrix_value_;
  this->large_matrix_value = large_matrix_value_;
  this->num_infinite_costs = 0;
  this->num_infinite_hessian_values = 0;
  this->num_infinite_col_bounds = 0;
  this->num_infinite_row_bounds = 0;
  this->num_small_matrix_values = 0;
  this->num_large_matrix_values = 0;
  this->suggested_user_objective_scale = 0;
  this->suggested_user_bound_scale = 0;
  this->applied = false;
}

bool HighsUserScaleData::scaleError(std::string& message) const {
  if (this->num_infinite_costs + this->num_infinite_hessian_values +
          this->num_infinite_col_bounds + this->num_infinite_row_bounds +
          this->num_large_matrix_values ==
      0)
    return false;
  assert(this->user_objective_scale != 0 || this->user_bound_scale != 0);
  std::stringstream ss;
  ss.str(std::string());
  ss << "User scaling of";
  if (this->user_objective_scale != 0) {
    ss << " 2**(" << this->user_objective_scale << ") for costs";
  }
  if (this->user_bound_scale != 0) {
    if (this->user_objective_scale != 0) ss << " and";
    ss << " 2**(" << this->user_bound_scale << ") for bounds";
  }
  ss << " yields";
  if (this->num_infinite_costs) {
    ss << " " << this->num_infinite_costs << " infinite cost";
    if (this->num_infinite_costs > 1) ss << "s";
  }
  if (this->num_infinite_hessian_values) {
    if (this->num_infinite_costs) {
      if (this->num_infinite_col_bounds || this->num_infinite_row_bounds) {
        ss << ",";
      } else {
        ss << " and";
      }
    }
    ss << " " << this->num_infinite_hessian_values
       << " infinite Hessian values";
    if (this->num_infinite_hessian_values > 1) ss << "s";
  }
  if (this->num_infinite_col_bounds) {
    if (this->num_infinite_costs || this->num_infinite_hessian_values) {
      if (this->num_infinite_row_bounds) {
        ss << ",";
      } else {
        ss << " and";
      }
    }
    ss << " " << this->num_infinite_col_bounds << " infinite column bound";
    if (this->num_infinite_col_bounds > 1) ss << "s";
  }
  if (this->num_infinite_row_bounds) {
    if (this->num_infinite_costs || this->num_infinite_hessian_values ||
        this->num_infinite_col_bounds)
      ss << " and";
    ss << " " << this->num_infinite_row_bounds << " infinite row bound";
    if (this->num_infinite_row_bounds > 1) ss << "s";
  }
  if (this->num_large_matrix_values) {
    if (this->num_infinite_costs + this->num_infinite_hessian_values +
            this->num_infinite_col_bounds + this->num_infinite_row_bounds >
        0)
      ss << ", and";
    ss << " " << this->num_large_matrix_values << " large matrix value";
    if (this->num_large_matrix_values > 1) ss << "s";
  }
  ss << "\n";
  message = ss.str();
  return true;
}

bool HighsUserScaleData::scaleWarning(std::string& message) const {
  if (this->num_small_matrix_values == 0) return false;
  assert(this->user_bound_scale != 0);
  std::stringstream ss;
  ss.str(std::string());
  ss << "User scaling of 2**(" << this->user_bound_scale
     << ") for bounds yields " << this->num_small_matrix_values
     << " small matrix value";
  if (this->num_small_matrix_values > 1) ss << "s";
  ss << "\n";
  message = ss.str();
  return true;
}

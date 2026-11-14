/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsInterface.cpp
 * @brief
 */
#include <sstream>
#include <unordered_set>

#include "Highs.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "mip/MipTimer.h"
#include "model/HighsHessianUtils.h"
#include "parallel/HighsParallel.h"
#include "simplex/HSimplex.h"
#include "util/HighsMatrixUtils.h"
#include "util/HighsSort.h"

void Highs::reportModelStats() const {
  const HighsLp& lp = this->model_.lp_;
  const HighsHessian& hessian = this->model_.hessian_;
  const HighsLogOptions& log_options = this->options_.log_options;
  if (!*log_options.output_flag) return;
  HighsInt num_integer = 0;
  HighsInt num_binary = 0;
  HighsInt num_semi_continuous = 0;
  HighsInt num_semi_integer = 0;
  for (HighsInt iCol = 0; iCol < static_cast<HighsInt>(lp.integrality_.size());
       iCol++) {
    switch (lp.integrality_[iCol]) {
      case HighsVarType::kInteger:
        num_integer++;
        if (lp.col_lower_[iCol] == 0 && lp.col_upper_[iCol] == 1) num_binary++;
        break;
      case HighsVarType::kSemiContinuous:
        num_semi_continuous++;
        break;
      case HighsVarType::kSemiInteger:
        num_semi_integer++;
        break;
      default:
        break;
    }
  }
  std::string problem_type;
  const bool non_continuous =
      num_integer + num_semi_continuous + num_semi_integer;
  if (model_.isQp()) {
    if (non_continuous) {
      problem_type = "MIQP";
    } else {
      problem_type = "QP";
    }
  } else {
    if (non_continuous) {
      problem_type = "MIP";
    } else {
      problem_type = "LP";
    }
  }
  const HighsInt a_num_nz = lp.a_matrix_.numNz();
  const HighsInt q_num_nz = hessian.dim_ > 0 ? hessian.numNz() : 0;
  if (*log_options.log_dev_level) {
    highsLogDev(log_options, HighsLogType::kInfo, "%-4s      : %s\n",
                problem_type.c_str(), lp.model_name_.c_str());
    highsLogDev(log_options, HighsLogType::kInfo,
                "Rows      : %" HIGHSINT_FORMAT "\n", lp.num_row_);
    highsLogDev(log_options, HighsLogType::kInfo,
                "Cols      : %" HIGHSINT_FORMAT "\n", lp.num_col_);
    if (q_num_nz) {
      highsLogDev(log_options, HighsLogType::kInfo,
                  "Matrix Nz : %" HIGHSINT_FORMAT "\n", a_num_nz);
      highsLogDev(log_options, HighsLogType::kInfo,
                  "Hessian Nz: %" HIGHSINT_FORMAT "\n", q_num_nz);
    } else {
      highsLogDev(log_options, HighsLogType::kInfo,
                  "Nonzeros  : %" HIGHSINT_FORMAT "\n", a_num_nz);
    }
    if (num_integer)
      highsLogDev(log_options, HighsLogType::kInfo,
                  "Integer   : %" HIGHSINT_FORMAT " (%" HIGHSINT_FORMAT
                  " binary)\n",
                  num_integer, num_binary);
    if (num_semi_continuous)
      highsLogDev(log_options, HighsLogType::kInfo,
                  "SemiConts : %" HIGHSINT_FORMAT "\n", num_semi_continuous);
    if (num_semi_integer)
      highsLogDev(log_options, HighsLogType::kInfo,
                  "SemiInt   : %" HIGHSINT_FORMAT "\n", num_semi_integer);
  } else {
    std::stringstream stats_line;
    stats_line << problem_type;
    if (lp.model_name_.length()) stats_line << " " << lp.model_name_;
    stats_line << " has " << lp.num_row_ << " row"
               << highsIntToPlural(lp.num_row_) << "; " << lp.num_col_ << " col"
               << highsIntToPlural(lp.num_col_);
    if (q_num_nz) {
      stats_line << "; " << a_num_nz << " matrix nonzero"
                 << highsIntToPlural(a_num_nz);
      stats_line << "; " << q_num_nz << " Hessian nonzero"
                 << highsIntToPlural(q_num_nz);
    } else {
      stats_line << "; " << a_num_nz << " nonzero"
                 << highsIntToPlural(a_num_nz);
    }
    if (hessian.isOracle()) stats_line << "; Hessian as oracle";
    if (num_integer)
      stats_line << "; " << num_integer << " integer variable"
                 << highsIntToPlural(num_integer) << " (" << num_binary
                 << " binary)";
    if (num_semi_continuous)
      stats_line << "; " << num_semi_continuous << " semi-continuous variable"
                 << highsIntToPlural(num_semi_continuous);
    if (num_semi_integer)
      stats_line << "; " << num_semi_integer << " semi-integer variable"
                 << highsIntToPlural(num_semi_integer);
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n",
                 stats_line.str().c_str());
  }
}

HighsStatus Highs::formStandardFormLp() {
  this->clearStandardFormLp();
  HighsLp& lp = this->model_.lp_;
  HighsSparseMatrix& matrix = lp.a_matrix_;
  // Ensure that the incumbent LP and standard form LP matrices are rowwise
  matrix.ensureRowwise();
  // Original rows are processed before columns, so that any original
  // boxed rows can be transformed to pairs of one-sided rows,
  // requiring the standard form matrix to be row-wise. The original
  // columns are assumed to come before any new columns, so their
  // costs (as a minimization) must be defined befor costs of new
  // columns.
  // Determine the objective scaling, and apply it to any offset
  HighsInt sense = HighsInt(lp.sense_);
  this->standard_form_offset_ = sense * lp.offset_;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    this->standard_form_cost_.push_back(sense * lp.col_cost_[iCol]);
  this->standard_form_matrix_.format_ = MatrixFormat::kRowwise;
  this->standard_form_matrix_.num_col_ = lp.num_col_;
  // Create a HighsSparseMatrix instance to store rows extracted from
  // the original constraint matrix
  HighsInt local_row_min_nnz = std::max(lp.num_col_, HighsInt(2));
  HighsSparseMatrix local_row;
  local_row.ensureRowwise();
  local_row.num_row_ = 1;
  local_row.num_col_ = lp.num_col_;
  local_row.index_.resize(local_row_min_nnz);
  local_row.value_.resize(local_row_min_nnz);
  local_row.start_.resize(2);
  HighsInt& num_nz = local_row.start_[1];
  local_row.start_[0] = 0;
  HighsInt num_fixed_row = 0;
  HighsInt num_boxed_row = 0;
  HighsInt num_lower_row = 0;
  HighsInt num_upper_row = 0;
  HighsInt num_free_row = 0;
  HighsInt num_fixed_col = 0;
  HighsInt num_boxed_col = 0;
  HighsInt num_lower_col = 0;
  HighsInt num_upper_col = 0;
  HighsInt num_free_col = 0;
  std::vector<HighsInt> slack_ix;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    double lower = lp.row_lower_[iRow];
    double upper = lp.row_upper_[iRow];
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      assert(0 == 1);
      // Free row
      num_free_row++;
      continue;
    }
    if (lower == upper) {
      // Equality row
      num_fixed_row++;
      matrix.getRow(iRow, num_nz, local_row.index_.data(),
                    local_row.value_.data());
      this->standard_form_matrix_.addRows(local_row);
      this->standard_form_rhs_.push_back(upper);
      continue;
    } else if (lower <= -kHighsInf) {
      // Upper bounded row, so record the slack
      num_upper_row++;
      assert(upper < kHighsInf);
      HighsInt standard_form_row = this->standard_form_rhs_.size();
      slack_ix.push_back(standard_form_row + 1);
      matrix.getRow(iRow, num_nz, local_row.index_.data(),
                    local_row.value_.data());
      this->standard_form_matrix_.addRows(local_row);
      this->standard_form_rhs_.push_back(upper);
    } else if (upper >= kHighsInf) {
      // Lower bounded row, so record the slack
      num_lower_row++;
      assert(lower > -kHighsInf);
      HighsInt standard_form_row = this->standard_form_rhs_.size();
      slack_ix.push_back(-(standard_form_row + 1));
      matrix.getRow(iRow, num_nz, local_row.index_.data(),
                    local_row.value_.data());
      this->standard_form_matrix_.addRows(local_row);
      this->standard_form_rhs_.push_back(lower);
    } else {
      // Boxed row, so record the lower slack
      assert(lower > -kHighsInf);
      assert(upper < kHighsInf);
      num_boxed_row++;
      HighsInt standard_form_row = this->standard_form_rhs_.size();
      slack_ix.push_back(-(standard_form_row + 1));
      matrix.getRow(iRow, num_nz, local_row.index_.data(),
                    local_row.value_.data());
      this->standard_form_matrix_.addRows(local_row);
      this->standard_form_rhs_.push_back(lower);
      // .. and upper slack, adding a copy of the row
      standard_form_row = this->standard_form_rhs_.size();
      slack_ix.push_back(standard_form_row + 1);
      this->standard_form_matrix_.addRows(local_row);
      this->standard_form_rhs_.push_back(upper);
    }
  }
  // Add rows corresponding to boxed columns
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    double lower = lp.col_lower_[iCol];
    double upper = lp.col_upper_[iCol];
    if (lower > -kHighsInf && upper < kHighsInf) {
      // Boxed column
      //
      // x will be replaced by x = l + X (below) with X >= 0
      //
      // Introduce variable s >= 0 so that (with x >= l still)
      //
      // x = u - s => x + s = u
      this->standard_form_cost_.push_back(0);
      this->standard_form_matrix_.num_col_++;
      local_row.num_col_++;
      local_row.index_[0] = iCol;
      local_row.index_[1] = this->standard_form_matrix_.num_col_ - 1;
      local_row.value_[0] = 1;
      local_row.value_[1] = 1;
      local_row.start_[1] = 2;
      this->standard_form_matrix_.addRows(local_row);
      this->standard_form_rhs_.push_back(upper);
    }
  }
  // Finished with both matrices, row-wise, so ensure that the
  // incumbent matrix leaves col-wise, and that the standard form
  // matrix is col-wise so RHS shifts can be applied and more columns
  // can be added
  matrix.ensureColwise();
  this->standard_form_matrix_.ensureColwise();
  // Work through the columns, ensuring that all have non-negativity bounds
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    double cost = sense * lp.col_cost_[iCol];
    double lower = lp.col_lower_[iCol];
    double upper = lp.col_upper_[iCol];
    if (lower > -kHighsInf) {
      // Finite lower bound
      if (upper < kHighsInf) {
        if (lower == upper) {
          // Fixed column
          num_fixed_col++;
        } else {
          // Boxed column
          num_boxed_col++;
        }
      } else {
        // Lower column
        num_lower_col++;
      }
      if (lower != 0) {
        // x >= l, so shift x-l = X >= 0, giving x = X + l
        //
        // Cost contribution c(X+l) = cX + cl
        this->standard_form_offset_ += cost * lower;
        // Constraint contribution a(X+l) = aX + al
        for (HighsInt iEl = this->standard_form_matrix_.start_[iCol];
             iEl < this->standard_form_matrix_.start_[iCol + 1]; iEl++)
          this->standard_form_rhs_[this->standard_form_matrix_.index_[iEl]] -=
              this->standard_form_matrix_.value_[iEl] * lower;
      }
    } else if (upper < kHighsInf) {
      // Upper column
      num_upper_col++;
      // Have to operate even if u=0, since cost and column values are negated
      //
      // x <= u, so shift u-x = X >= 0, giving x = u - X
      //
      // Cost contribution c(u-X) = cu - cX
      this->standard_form_offset_ += cost * upper;
      this->standard_form_cost_[iCol] = -cost;
      // Constraint contribution a(u-X) = -aX + au
      for (HighsInt iEl = this->standard_form_matrix_.start_[iCol];
           iEl < this->standard_form_matrix_.start_[iCol + 1]; iEl++) {
        this->standard_form_rhs_[this->standard_form_matrix_.index_[iEl]] -=
            this->standard_form_matrix_.value_[iEl] * upper;
        this->standard_form_matrix_.value_[iEl] =
            -this->standard_form_matrix_.value_[iEl];
      }
    } else {
      // Free column
      num_free_col++;
      // Represent as x = x+ - x-
      //
      // where original column is now x+ >= 0
      //
      // and x- >= 0 has negation of its cost and matrix column
      this->standard_form_cost_.push_back(-cost);
      for (HighsInt iEl = this->standard_form_matrix_.start_[iCol];
           iEl < this->standard_form_matrix_.start_[iCol + 1]; iEl++) {
        this->standard_form_matrix_.index_.push_back(
            this->standard_form_matrix_.index_[iEl]);
        this->standard_form_matrix_.value_.push_back(
            -this->standard_form_matrix_.value_[iEl]);
      }
      this->standard_form_matrix_.start_.push_back(
          HighsInt(this->standard_form_matrix_.index_.size()));
    }
  }
  // Now add the slack variables
  for (HighsInt iRow : slack_ix) {
    this->standard_form_cost_.push_back(0);
    if (iRow > 0) {
      this->standard_form_matrix_.index_.push_back(iRow - 1);
      this->standard_form_matrix_.value_.push_back(1);
    } else {
      this->standard_form_matrix_.index_.push_back(-iRow - 1);
      this->standard_form_matrix_.value_.push_back(-1);
    }
    this->standard_form_matrix_.start_.push_back(
        HighsInt(this->standard_form_matrix_.index_.size()));
  }
  // Now set correct values for the dimensions of this->standard_form_matrix_
  this->standard_form_matrix_.num_col_ = int(standard_form_cost_.size());
  this->standard_form_matrix_.num_row_ = int(standard_form_rhs_.size());
  this->standard_form_valid_ = true;
  highsLogUser(options_.log_options, HighsLogType::kInfo,
               "Standard form LP obtained for LP with (free / lower / upper / "
               "boxed / fixed) variables"
               " (%d / %d / %d / %d / %d) and constraints"
               " (%d / %d / %d / %d / %d) \n",
               int(num_free_col), int(num_lower_col), int(num_upper_col),
               int(num_boxed_col), int(num_fixed_col), int(num_free_row),
               int(num_lower_row), int(num_upper_row), int(num_boxed_row),
               int(num_fixed_row));
  return HighsStatus::kOk;
}

HighsStatus Highs::basisForSolution() {
  HighsLp& lp = model_.lp_;
  assert(!lp.isMip() || options_.solve_relaxation);
  assert(solution_.value_valid);
  invalidateBasis();
  HighsInt num_basic = 0;
  HighsBasis basis;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (std::fabs(lp.col_lower_[iCol] - solution_.col_value[iCol]) <=
        options_.primal_feasibility_tolerance) {
      basis.col_status.push_back(HighsBasisStatus::kLower);
    } else if (std::fabs(lp.col_upper_[iCol] - solution_.col_value[iCol]) <=
               options_.primal_feasibility_tolerance) {
      basis.col_status.push_back(HighsBasisStatus::kUpper);
    } else {
      num_basic++;
      basis.col_status.push_back(HighsBasisStatus::kBasic);
    }
  }
  const HighsInt num_basic_col = num_basic;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (std::fabs(lp.row_lower_[iRow] - solution_.row_value[iRow]) <=
        options_.primal_feasibility_tolerance) {
      basis.row_status.push_back(HighsBasisStatus::kLower);
    } else if (std::fabs(lp.row_upper_[iRow] - solution_.row_value[iRow]) <=
               options_.primal_feasibility_tolerance) {
      basis.row_status.push_back(HighsBasisStatus::kUpper);
    } else {
      num_basic++;
      basis.row_status.push_back(HighsBasisStatus::kBasic);
    }
  }
  const HighsInt num_basic_row = num_basic - num_basic_col;
  assert((int)basis.col_status.size() == lp.num_col_);
  assert((int)basis.row_status.size() == lp.num_row_);
  highsLogDev(options_.log_options, HighsLogType::kInfo,
              "LP has %d rows and solution yields %d possible basic variables "
              "(%d / %d; %d / %d)\n",
              (int)lp.num_row_, (int)num_basic, (int)num_basic_col,
              (int)lp.num_col_, (int)num_basic_row, (int)lp.num_row_);
  return this->setBasis(basis);
}

HighsStatus Highs::addColsInterface(
    HighsInt ext_num_new_col, const double* ext_col_cost,
    const double* ext_col_lower, const double* ext_col_upper,
    HighsInt ext_num_new_nz, const HighsInt* ext_a_start,
    const HighsInt* ext_a_index, const double* ext_a_value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsOptions& options = options_;
  if (ext_num_new_col < 0) return HighsStatus::kError;
  if (ext_num_new_nz < 0) return HighsStatus::kError;
  if (ext_num_new_col == 0) return HighsStatus::kOk;
  if (ext_num_new_col > 0)
    if (isColDataNull(options.log_options, ext_col_cost, ext_col_lower,
                      ext_col_upper))
      return HighsStatus::kError;
  if (ext_num_new_nz > 0)
    if (isMatrixDataNull(options.log_options, ext_a_start, ext_a_index,
                         ext_a_value))
      return HighsStatus::kError;

  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsScale& scale = lp.scale_;
  bool& useful_basis = basis.useful;
  bool& lp_has_scaling = lp.scale_.has_scaling;

  // Check that if nonzeros are to be added then the model has a positive number
  // of rows
  if (lp.num_row_ <= 0 && ext_num_new_nz > 0) return HighsStatus::kError;

  // Record the new number of columns
  HighsInt newNumCol = lp.num_col_ + ext_num_new_col;

  HighsIndexCollection index_collection;
  index_collection.dimension_ = ext_num_new_col;
  index_collection.is_interval_ = true;
  index_collection.from_ = 0;
  index_collection.to_ = ext_num_new_col - 1;

  // Take a copy of the cost and bounds that can be normalised
  std::vector<double> local_colCost{ext_col_cost,
                                    ext_col_cost + ext_num_new_col};
  std::vector<double> local_colLower{ext_col_lower,
                                     ext_col_lower + ext_num_new_col};
  std::vector<double> local_colUpper{ext_col_upper,
                                     ext_col_upper + ext_num_new_col};

  bool local_has_infinite_cost = false;
  return_status = interpretCallStatus(
      options_.log_options,
      assessCosts(options, lp.num_col_, index_collection, local_colCost,
                  local_has_infinite_cost, options.infinite_cost),
      return_status, "assessCosts");
  if (return_status == HighsStatus::kError) return return_status;
  // Assess the column bounds
  return_status = interpretCallStatus(
      options_.log_options,
      assessBounds(options, "Col", lp.num_col_, index_collection,
                   local_colLower, local_colUpper, options.infinite_bound,
                   nullptr),
      return_status, "assessBounds");
  if (return_status == HighsStatus::kError) return return_status;
  // Form a column-wise HighsSparseMatrix of the new matrix columns so
  // that is easy to handle and, if there are nonzeros, it can be
  // normalised
  HighsSparseMatrix local_a_matrix;
  local_a_matrix.name_ = "local_a_matrix";
  local_a_matrix.num_col_ = ext_num_new_col;
  local_a_matrix.num_row_ = lp.num_row_;
  local_a_matrix.format_ = MatrixFormat::kColwise;
  if (ext_num_new_nz) {
    local_a_matrix.start_ = {ext_a_start, ext_a_start + ext_num_new_col};
    local_a_matrix.start_.resize(ext_num_new_col + 1);
    local_a_matrix.start_[ext_num_new_col] = ext_num_new_nz;
    local_a_matrix.index_ = {ext_a_index, ext_a_index + ext_num_new_nz};
    local_a_matrix.value_ = {ext_a_value, ext_a_value + ext_num_new_nz};
    const bool sum_duplicates = false;
    // Assess the matrix rows
    return_status = interpretCallStatus(
        options_.log_options,
        local_a_matrix.assess(options.log_options, options.small_matrix_value,
                              options.large_matrix_value, sum_duplicates),
        return_status, "local_a_matrix.assess");
    if (return_status == HighsStatus::kError) return return_status;
  } else {
    // No nonzeros so, whether the constraint matrix is column-wise or
    // row-wise, adding the empty matrix is trivial. Complete the
    // setup of an empty column-wise HighsSparseMatrix of the new
    // matrix columns
    local_a_matrix.start_.assign(ext_num_new_col + 1, 0);
  }
  assert(return_status != HighsStatus::kError);
  // Now that data for new columns has been assessed and no errors can
  // occur, append the columns to the LP vectors and matrix
  appendColsToLpVectors(lp, ext_num_new_col, local_colCost, local_colLower,
                        local_colUpper);
  lp.a_matrix_.addCols(local_a_matrix);
  if (lp_has_scaling) {
    // Extend the column scaling factors
    scale.col.resize(newNumCol);
    for (HighsInt iCol = 0; iCol < ext_num_new_col; iCol++)
      scale.col[lp.num_col_ + iCol] = 1.0;
    scale.num_col = newNumCol;
    // Apply the existing row scaling to the new columns
    local_a_matrix.applyRowScale(scale);
    // Consider applying column scaling to the new columns.
    local_a_matrix.considerColScaling(options.allowed_matrix_scale_factor,
                                      &scale.col[lp.num_col_]);
  }
  // Update the basis corresponding to new nonbasic columns
  if (useful_basis) appendNonbasicColsToBasisInterface(ext_num_new_col);

  // Possibly add blank column names
  lp.addColNames("", ext_num_new_col);

  // Increase the number of columns in the LP
  lp.num_col_ += ext_num_new_col;
  assert(lp.legalDimensions("addCols", options.log_options));

  // Interpret possible introduction of infinite costs
  lp.has_infinite_cost_ = lp.has_infinite_cost_ || local_has_infinite_cost;
  assert(lp.has_infinite_cost_ == lp.hasInfiniteCost(options.infinite_cost));

  // Deduce the consequences of adding new columns
  invalidateModelStatusSolutionAndInfo();

  // Determine any implications for simplex data
  ekk_instance_.addCols(lp, local_a_matrix);

  // Extend any Hessian with zeros on the diagonal
  if (this->model_.hessian_.dim_)
    completeHessian(lp.num_col_, this->model_.hessian_);
  return return_status;
}

HighsStatus Highs::addRowsInterface(HighsInt ext_num_new_row,
                                    const double* ext_row_lower,
                                    const double* ext_row_upper,
                                    HighsInt ext_num_new_nz,
                                    const HighsInt* ext_ar_start,
                                    const HighsInt* ext_ar_index,
                                    const double* ext_ar_value) {
  // addRows is fundamentally different from addCols, since the new
  // matrix data are held row-wise, so we have to insert data into the
  // column-wise matrix of the LP.
  if (kExtendInvertWhenAddingRows) {
    if (ekk_instance_.status_.has_nla)
      ekk_instance_.debugNlaCheckInvert("Start of Highs::addRowsInterface",
                                        kHighsDebugLevelExpensive + 1);
  }
  HighsStatus return_status = HighsStatus::kOk;
  HighsOptions& options = options_;
  if (ext_num_new_row < 0) return HighsStatus::kError;
  if (ext_num_new_nz < 0) return HighsStatus::kError;
  if (ext_num_new_row == 0) return HighsStatus::kOk;
  if (ext_num_new_row > 0)
    if (isRowDataNull(options.log_options, ext_row_lower, ext_row_upper))
      return HighsStatus::kError;
  if (ext_num_new_nz > 0)
    if (isMatrixDataNull(options.log_options, ext_ar_start, ext_ar_index,
                         ext_ar_value))
      return HighsStatus::kError;

  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsScale& scale = lp.scale_;
  bool& useful_basis = basis.useful;

  bool& lp_has_scaling = lp.scale_.has_scaling;

  // Check that if nonzeros are to be added then the model has a positive number
  // of columns
  if (lp.num_col_ <= 0 && ext_num_new_nz > 0) return HighsStatus::kError;

  // Record the new number of rows
  HighsInt newNumRow = lp.num_row_ + ext_num_new_row;

  HighsIndexCollection index_collection;
  index_collection.dimension_ = ext_num_new_row;
  index_collection.is_interval_ = true;
  index_collection.from_ = 0;
  index_collection.to_ = ext_num_new_row - 1;
  // Take a copy of the bounds that can be normalised
  std::vector<double> local_rowLower{ext_row_lower,
                                     ext_row_lower + ext_num_new_row};
  std::vector<double> local_rowUpper{ext_row_upper,
                                     ext_row_upper + ext_num_new_row};

  return_status = interpretCallStatus(
      options_.log_options,
      assessBounds(options, "Row", lp.num_row_, index_collection,
                   local_rowLower, local_rowUpper, options.infinite_bound,
                   nullptr),
      return_status, "assessBounds");
  if (return_status == HighsStatus::kError) return return_status;
  // Form a row-wise HighsSparseMatrix of the new matrix rows so that
  // is easy to handle and, if there are nonzeros, it can be
  // normalised
  HighsSparseMatrix local_ar_matrix;
  local_ar_matrix.name_ = "local_ar_matrix";
  local_ar_matrix.num_col_ = lp.num_col_;
  local_ar_matrix.num_row_ = ext_num_new_row;
  local_ar_matrix.format_ = MatrixFormat::kRowwise;
  if (ext_num_new_nz) {
    local_ar_matrix.start_ = {ext_ar_start, ext_ar_start + ext_num_new_row};
    local_ar_matrix.start_.resize(ext_num_new_row + 1);
    local_ar_matrix.start_[ext_num_new_row] = ext_num_new_nz;
    local_ar_matrix.index_ = {ext_ar_index, ext_ar_index + ext_num_new_nz};
    local_ar_matrix.value_ = {ext_ar_value, ext_ar_value + ext_num_new_nz};
    // Assess the matrix columns
    const bool sum_duplicates = false;
    return_status = interpretCallStatus(
        options_.log_options,
        local_ar_matrix.assess(options.log_options, options.small_matrix_value,
                               options.large_matrix_value, sum_duplicates),
        return_status, "local_ar_matrix.assess");
    if (return_status == HighsStatus::kError) return return_status;
  } else {
    // No nonzeros so, whether the constraint matrix is row-wise or
    // column-wise, adding the empty matrix is trivial. Complete the
    // setup of an empty row-wise HighsSparseMatrix of the new matrix
    // rows
    local_ar_matrix.start_.assign(ext_num_new_row + 1, 0);
  }
  assert(return_status != HighsStatus::kError);
  // Now that data for new rows has been assessed and no errors can
  // occur, append the rows to the LP vectors and matrix
  appendRowsToLpVectors(lp, ext_num_new_row, local_rowLower, local_rowUpper);
  lp.a_matrix_.addRows(local_ar_matrix);
  if (lp_has_scaling) {
    // Extend the row scaling factors
    scale.row.resize(newNumRow);
    for (HighsInt iRow = 0; iRow < ext_num_new_row; iRow++)
      scale.row[lp.num_row_ + iRow] = 1.0;
    scale.num_row = newNumRow;
    // Apply the existing column scaling to the new rows
    local_ar_matrix.applyColScale(scale);
    // Consider applying row scaling to the new rows.
    local_ar_matrix.considerRowScaling(options.allowed_matrix_scale_factor,
                                       &scale.row[lp.num_row_]);
  }
  // Update the basis corresponding to new basic rows
  if (useful_basis) appendBasicRowsToBasisInterface(ext_num_new_row);

  // Possibly add blank row names
  lp.addRowNames("", ext_num_new_row);

  // Increase the number of rows in the LP
  lp.num_row_ += ext_num_new_row;
  assert(lp.legalDimensions("addRows", options.log_options));

  // Deduce the consequences of adding new rows
  invalidateModelStatusSolutionAndInfo();
  // Determine any implications for simplex data
  ekk_instance_.addRows(lp, local_ar_matrix);

  return return_status;
}

static void deleteBasisEntries(std::vector<HighsBasisStatus>& status,
                               bool& deleted_basic, bool& deleted_nonbasic,
                               const HighsIndexCollection& index_collection,
                               const HighsInt entry_dim) {
  assert(ok(index_collection));
  assert(static_cast<size_t>(entry_dim) == status.size());
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return;

  HighsInt delete_from_entry;
  HighsInt delete_to_entry;
  HighsInt keep_from_entry;
  HighsInt keep_to_entry = -1;
  HighsInt current_set_entry = 0;
  HighsInt new_num_entry = 0;
  deleted_basic = false;
  deleted_nonbasic = false;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_entry, delete_to_entry,
                     keep_from_entry, keep_to_entry, current_set_entry);
    // Account for the initial entries being kept
    if (k == from_k) new_num_entry = delete_from_entry;
    // Identify whether a basic or a nonbasic entry has been deleted
    for (HighsInt entry = delete_from_entry; entry <= delete_to_entry;
         entry++) {
      if (status[entry] == HighsBasisStatus::kBasic) {
        deleted_basic = true;
      } else {
        deleted_nonbasic = true;
      }
    }
    if (delete_to_entry >= entry_dim - 1) break;
    for (HighsInt entry = keep_from_entry; entry <= keep_to_entry; entry++) {
      status[new_num_entry] = status[entry];
      new_num_entry++;
    }
    if (keep_to_entry >= entry_dim - 1) break;
  }
  status.resize(new_num_entry);
}

static void deleteBasisCols(HighsBasis& basis,
                            const HighsIndexCollection& index_collection,
                            const HighsInt original_num_col) {
  bool deleted_basic;
  bool deleted_nonbasic;
  deleteBasisEntries(basis.col_status, deleted_basic, deleted_nonbasic,
                     index_collection, original_num_col);
  if (deleted_basic) basis.valid = false;
}

static void deleteBasisRows(HighsBasis& basis,
                            const HighsIndexCollection& index_collection,
                            const HighsInt original_num_row) {
  bool deleted_basic;
  bool deleted_nonbasic;
  deleteBasisEntries(basis.row_status, deleted_basic, deleted_nonbasic,
                     index_collection, original_num_row);
  if (deleted_nonbasic) basis.valid = false;
}

void Highs::deleteColsInterface(HighsIndexCollection& index_collection) {
  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  lp.ensureColwise();

  // Keep a copy of the original number of columns to check whether
  // any columns have been removed, and if there is mask to be updated
  HighsInt original_num_col = lp.num_col_;

  lp.deleteCols(index_collection);
  model_.hessian_.deleteCols(index_collection);
  // Bail out if no columns were actually deleted
  if (lp.num_col_ == original_num_col) return;

  assert(lp.num_col_ < original_num_col);

  // Nontrivial deletion so reset the model_status and update any
  // Highs basis
  model_status_ = HighsModelStatus::kNotset;
  if (basis_.useful) {
    assert(basis_.col_status.size() == static_cast<size_t>(original_num_col));
    // Have a full set of column basis status values, so maintain
    // them, and only invalidate the basis if a basic column has been
    // deleted
    deleteBasisCols(basis_, index_collection, original_num_col);
  } else {
    assert(!basis.valid);
  }

  if (lp.scale_.has_scaling) {
    deleteScale(lp.scale_.col, index_collection);
    lp.scale_.col.resize(lp.num_col_);
    lp.scale_.num_col = lp.num_col_;
  }
  // Deduce the consequences of deleting columns
  invalidateModelStatusSolutionAndInfo();

  // Determine any implications for simplex data
  ekk_instance_.deleteCols(index_collection);

  if (index_collection.is_mask_) {
    // Set the mask values to indicate the new index value of the
    // remaining columns
    HighsInt new_col = 0;
    for (HighsInt col = 0; col < original_num_col; col++) {
      if (!index_collection.mask_[col]) {
        index_collection.mask_[col] = new_col;
        new_col++;
      } else {
        index_collection.mask_[col] = -1;
      }
    }
    assert(new_col == lp.num_col_);
  }
  assert(lp.legalDimensions("deleteCols", options_.log_options));
  lp.col_hash_.name2index.clear();
}

void Highs::deleteRowsInterface(HighsIndexCollection& index_collection) {
  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  lp.ensureColwise();
  // Keep a copy of the original number of rows to check whether
  // any rows have been removed, and if there is mask to be updated
  HighsInt original_num_row = lp.num_row_;

  lp.deleteRows(index_collection);
  // Bail out if no rows were actually deleted
  if (lp.num_row_ == original_num_row) return;

  assert(lp.num_row_ < original_num_row);

  // Nontrivial deletion so reset the model_status and update any
  // Highs basis
  model_status_ = HighsModelStatus::kNotset;
  if (basis_.useful) {
    assert(basis_.row_status.size() == static_cast<size_t>(original_num_row));
    // Have a full set of row basis status values, so maintain them,
    // and only invalidate the basis if a nonbasic row has been
    // deleted
    deleteBasisRows(basis_, index_collection, original_num_row);
  } else {
    assert(!basis.valid);
  }

  if (lp.scale_.has_scaling) {
    deleteScale(lp.scale_.row, index_collection);
    lp.scale_.row.resize(lp.num_row_);
    lp.scale_.num_row = lp.num_row_;
  }
  // Deduce the consequences of deleting rows
  invalidateModelStatusSolutionAndInfo();

  // Determine any implications for simplex data
  ekk_instance_.deleteRows(index_collection);
  if (index_collection.is_mask_) {
    HighsInt new_row = 0;
    for (HighsInt row = 0; row < original_num_row; row++) {
      if (!index_collection.mask_[row]) {
        index_collection.mask_[row] = new_row;
        new_row++;
      } else {
        index_collection.mask_[row] = -1;
      }
    }
    assert(new_row == lp.num_row_);
  }
  assert(lp.legalDimensions("deleteRows", options_.log_options));
  lp.row_hash_.name2index.clear();
}

void Highs::getColsInterface(const HighsIndexCollection& index_collection,
                             HighsInt& num_col, double* cost, double* lower,
                             double* upper, HighsInt& num_nz, HighsInt* start,
                             HighsInt* index, double* value) const {
  const HighsLp& lp = model_.lp_;
  if (lp.a_matrix_.isColwise()) {
    getSubVectors(index_collection, lp.num_col_, lp.col_cost_.data(),
                  lp.col_lower_.data(), lp.col_upper_.data(), lp.a_matrix_,
                  num_col, cost, lower, upper, num_nz, start, index, value);
  } else {
    getSubVectorsTranspose(index_collection, lp.num_col_, lp.col_cost_.data(),
                           lp.col_lower_.data(), lp.col_upper_.data(),
                           lp.a_matrix_, num_col, cost, lower, upper, num_nz,
                           start, index, value);
  }
}

void Highs::getRowsInterface(const HighsIndexCollection& index_collection,
                             HighsInt& num_row, double* lower, double* upper,
                             HighsInt& num_nz, HighsInt* start, HighsInt* index,
                             double* value) const {
  const HighsLp& lp = model_.lp_;
  if (lp.a_matrix_.isColwise()) {
    getSubVectorsTranspose(index_collection, lp.num_row_, nullptr,
                           lp.row_lower_.data(), lp.row_upper_.data(),
                           lp.a_matrix_, num_row, nullptr, lower, upper, num_nz,
                           start, index, value);
  } else {
    getSubVectors(index_collection, lp.num_row_, nullptr, lp.row_lower_.data(),
                  lp.row_upper_.data(), lp.a_matrix_, num_row, nullptr, lower,
                  upper, num_nz, start, index, value);
  }
}

void Highs::getCoeffInterface(const HighsInt ext_row, const HighsInt ext_col,
                              double& value) const {
  const HighsLp& lp = model_.lp_;
  assert(0 <= ext_row && ext_row < lp.num_row_);
  assert(0 <= ext_col && ext_col < lp.num_col_);
  value = 0;

  if (lp.a_matrix_.isColwise()) {
    for (HighsInt el = lp.a_matrix_.start_[ext_col];
         el < lp.a_matrix_.start_[ext_col + 1]; el++) {
      if (lp.a_matrix_.index_[el] == ext_row) {
        value = lp.a_matrix_.value_[el];
        break;
      }
    }
  } else {
    for (HighsInt el = lp.a_matrix_.start_[ext_row];
         el < lp.a_matrix_.start_[ext_row + 1]; el++) {
      if (lp.a_matrix_.index_[el] == ext_col) {
        value = lp.a_matrix_.value_[el];
        break;
      }
    }
  }
}

HighsStatus Highs::changeIntegralityInterface(
    HighsIndexCollection& index_collection, const HighsVarType* integrality) {
  HighsInt num_integrality = dataSize(index_collection);
  // If a non-positive number of integrality (may) need changing nothing needs
  // to be done
  if (num_integrality <= 0) return HighsStatus::kOk;
  if (highsVarTypeUserDataNotNull(options_.log_options, integrality,
                                  "column integrality"))
    return HighsStatus::kError;
  // Take a copy of the integrality that can be normalised
  std::vector<HighsVarType> local_integrality{integrality,
                                              integrality + num_integrality};
  // If changing the integrality for a set of columns, verify that the
  // set entries are in ascending order
  if (index_collection.is_set_)
    assert(increasingSetOk(index_collection.set_, 0,
                           index_collection.dimension_, true));
  HighsStatus return_status = changeLpIntegrality(model_.lp_, index_collection,
                                                  local_integrality, options_);
  assert(return_status != HighsStatus::kError);
  // Deduce the consequences of new integrality
  invalidateModelStatus();
  return return_status;
}

HighsStatus Highs::changeCostsInterface(HighsIndexCollection& index_collection,
                                        const double* cost) {
  HighsInt num_cost = dataSize(index_collection);
  // If a non-positive number of costs (may) need changing nothing needs to be
  // done
  if (num_cost <= 0) return HighsStatus::kOk;
  if (doubleUserDataNotNull(options_.log_options, cost, "column costs"))
    return HighsStatus::kError;
  // Take a copy of the cost that can be normalised
  std::vector<double> local_colCost{cost, cost + num_cost};
  HighsStatus return_status = HighsStatus::kOk;
  bool local_has_infinite_cost = false;
  return_status = interpretCallStatus(
      options_.log_options,
      assessCosts(options_, 0, index_collection, local_colCost,
                  local_has_infinite_cost, options_.infinite_cost),
      return_status, "assessCosts");
  if (return_status == HighsStatus::kError) return return_status;
  HighsLp& lp = model_.lp_;
  changeLpCosts(lp, index_collection, local_colCost, options_.infinite_cost);

  // Interpret possible introduction of infinite costs
  lp.has_infinite_cost_ = lp.has_infinite_cost_ || local_has_infinite_cost;
  assert(lp.has_infinite_cost_ == lp.hasInfiniteCost(options_.infinite_cost));

  // Deduce the consequences of new costs
  invalidateModelStatusSolutionAndInfo();
  // Determine any implications for simplex data
  ekk_instance_.updateStatus(LpAction::kNewCosts);
  return HighsStatus::kOk;
}

bool Highs::feasibleWrtBounds(const bool columns) const {
  if (this->info_.primal_solution_status != kSolutionStatusFeasible)
    return false;
  const HighsLp& lp = model_.lp_;
  const double primal_feasibility_tolerance =
      this->options_.primal_feasibility_tolerance;
  const std::vector<double>& value =
      columns ? this->solution_.col_value : this->solution_.row_value;
  const std::vector<double>& lower = columns ? lp.col_lower_ : lp.row_lower_;
  const std::vector<double>& upper = columns ? lp.col_upper_ : lp.row_upper_;
  HighsInt dim = columns ? lp.num_col_ : lp.num_row_;
  for (HighsInt iX = 0; iX < dim; iX++) {
    if (value[iX] < lower[iX] - primal_feasibility_tolerance) return false;
    if (value[iX] > upper[iX] + primal_feasibility_tolerance) return false;
  }
  return true;
}

HighsStatus Highs::changeColBoundsInterface(
    HighsIndexCollection& index_collection, const double* col_lower,
    const double* col_upper) {
  HighsInt num_col_bounds = dataSize(index_collection);
  // If a non-positive number of costs (may) need changing nothing needs to be
  // done
  if (num_col_bounds <= 0) return HighsStatus::kOk;
  bool null_data = false;
  null_data = doubleUserDataNotNull(options_.log_options, col_lower,
                                    "column lower bounds") ||
              null_data;
  null_data = doubleUserDataNotNull(options_.log_options, col_upper,
                                    "column upper bounds") ||
              null_data;
  if (null_data) return HighsStatus::kError;
  // Take a copy of the cost that can be normalised
  std::vector<double> local_colLower{col_lower, col_lower + num_col_bounds};
  std::vector<double> local_colUpper{col_upper, col_upper + num_col_bounds};
  // If changing the bounds for a set of columns, ensure that the
  // set and data are in ascending order
  if (index_collection.is_set_)
    sortSetData(index_collection.set_num_entries_, index_collection.set_,
                col_lower, col_upper, NULL, local_colLower.data(),
                local_colUpper.data(), NULL);
  HighsStatus return_status = HighsStatus::kOk;
  return_status = interpretCallStatus(
      options_.log_options,
      assessBounds(options_, "col", 0, index_collection, local_colLower,
                   local_colUpper, options_.infinite_bound, nullptr),
      return_status, "assessBounds");
  if (return_status == HighsStatus::kError) return return_status;
  HighsLp& lp = model_.lp_;

  changeLpColBounds(lp, index_collection, local_colLower, local_colUpper);
  // Update HiGHS basis status and (any) simplex move status of
  // nonbasic variables whose bounds have changed
  setNonbasicStatusInterface(index_collection, true);
  // Deduce the consequences of new col bounds
  if (!this->basis_.useful && feasibleWrtBounds()) {
    // Retain the solution if there's no basis, and the solution is
    // feasible
    invalidateModelStatusAndInfo();
  } else {
    // Invalidate the solution
    invalidateModelStatusSolutionAndInfo();
  }
  // Determine any implications for simplex data
  ekk_instance_.updateStatus(LpAction::kNewBounds);
  return return_status;
}

HighsStatus Highs::changeRowBoundsInterface(
    HighsIndexCollection& index_collection, const double* lower,
    const double* upper) {
  HighsInt num_row_bounds = dataSize(index_collection);
  // If a non-positive number of costs (may) need changing nothing needs to be
  // done
  if (num_row_bounds <= 0) return HighsStatus::kOk;
  bool null_data = false;
  null_data =
      doubleUserDataNotNull(options_.log_options, lower, "row lower bounds") ||
      null_data;
  null_data =
      doubleUserDataNotNull(options_.log_options, upper, "row upper bounds") ||
      null_data;
  if (null_data) return HighsStatus::kError;
  // Take a copy of the cost that can be normalised
  std::vector<double> local_rowLower{lower, lower + num_row_bounds};
  std::vector<double> local_rowUpper{upper, upper + num_row_bounds};
  // If changing the bounds for a set of rows, ensure that the
  // set and data are in ascending order
  if (index_collection.is_set_)
    sortSetData(index_collection.set_num_entries_, index_collection.set_, lower,
                upper, NULL, local_rowLower.data(), local_rowUpper.data(),
                NULL);
  HighsStatus return_status = HighsStatus::kOk;
  return_status = interpretCallStatus(
      options_.log_options,
      assessBounds(options_, "row", 0, index_collection, local_rowLower,
                   local_rowUpper, options_.infinite_bound, nullptr),
      return_status, "assessBounds");
  if (return_status == HighsStatus::kError) return return_status;
  HighsLp& lp = model_.lp_;

  changeLpRowBounds(lp, index_collection, local_rowLower, local_rowUpper);
  // Update HiGHS basis status and (any) simplex move status of
  // nonbasic variables whose bounds have changed
  setNonbasicStatusInterface(index_collection, false);
  // Deduce the consequences of new row bounds
  if (!this->basis_.useful && feasibleWrtBounds(false)) {
    // Retain the solution if there's no basis, and the solution is
    // feasible
    invalidateModelStatusAndInfo();
  } else {
    // Invalidate the solution
    invalidateModelStatusSolutionAndInfo();
  }
  // Determine any implications for simplex data
  ekk_instance_.updateStatus(LpAction::kNewBounds);
  return return_status;
}

// Change a single coefficient in the matrix
HighsStatus Highs::changeCoeffInterface(const HighsInt ext_row,
                                        const HighsInt ext_col,
                                        const double ext_new_value) {
  HighsLp& lp = model_.lp_;
  // Ensure that the LP is column-wise
  lp.ensureColwise();
  // Should alread yhave checked that the row and column indices are
  // legal, and that the value is not large
  assert(0 <= ext_row && ext_row < lp.num_row_);
  assert(0 <= ext_col && ext_col < lp.num_col_);
  const double abs_ext_new_value = std::fabs(ext_new_value);
  assert(abs_ext_new_value < options_.large_matrix_value);

  HighsStatus return_status =
      changeLpMatrixCoeff(options_.log_options, lp, ext_row, ext_col,
                          ext_new_value, options_.small_matrix_value);
  assert(return_status != HighsStatus::kError);
  // Deduce the consequences of a changed element
  //
  // ToDo: Can do something more intelligent if element is in nonbasic
  // column
  //
  const bool basic_column =
      this->basis_.col_status[ext_col] == HighsBasisStatus::kBasic;
  //
  // For now, treat it as if it's a new row
  invalidateModelStatusSolutionAndInfo();

  if (basic_column) {
    // Basis is retained, but is has to be viewed as alien, since the
    // basis matrix has changed
    this->basis_.was_alien = true;
    this->basis_.alien = true;
  }

  // Determine any implications for simplex data
  ekk_instance_.updateStatus(LpAction::kNewRows);
  return return_status;
}

HighsStatus Highs::scaleColInterface(const HighsInt col,
                                     const double scale_value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsSimplexStatus& simplex_status = ekk_instance_.status_;

  // Ensure that the LP is column-wise
  lp.ensureColwise();
  if (col < 0) return HighsStatus::kError;
  if (col >= lp.num_col_) return HighsStatus::kError;
  if (!scale_value) return HighsStatus::kError;

  return_status = interpretCallStatus(options_.log_options,
                                      applyScalingToLpCol(lp, col, scale_value),
                                      return_status, "applyScalingToLpCol");
  if (return_status == HighsStatus::kError) return return_status;

  if (scale_value < 0 && basis.valid) {
    // Negative, so flip any nonbasic status
    if (basis.col_status[col] == HighsBasisStatus::kLower) {
      basis.col_status[col] = HighsBasisStatus::kUpper;
    } else if (basis.col_status[col] == HighsBasisStatus::kUpper) {
      basis.col_status[col] = HighsBasisStatus::kLower;
    }
  }
  if (simplex_status.initialised_for_solve) {
    SimplexBasis& simplex_basis = ekk_instance_.basis_;
    if (scale_value < 0 && simplex_status.has_basis) {
      // Negative, so flip any nonbasic status
      if (simplex_basis.nonbasicMove_[col] == kNonbasicMoveUp) {
        simplex_basis.nonbasicMove_[col] = kNonbasicMoveDn;
      } else if (simplex_basis.nonbasicMove_[col] == kNonbasicMoveDn) {
        simplex_basis.nonbasicMove_[col] = kNonbasicMoveUp;
      }
    }
  }
  // Deduce the consequences of a scaled column
  invalidateModelStatusSolutionAndInfo();

  // Determine any implications for simplex data
  ekk_instance_.updateStatus(LpAction::kScaledCol);
  return HighsStatus::kOk;
}

HighsStatus Highs::scaleRowInterface(const HighsInt row,
                                     const double scale_value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsSimplexStatus& simplex_status = ekk_instance_.status_;

  // Ensure that the LP is column-wise
  lp.ensureColwise();

  if (row < 0) return HighsStatus::kError;
  if (row >= lp.num_row_) return HighsStatus::kError;
  if (!scale_value) return HighsStatus::kError;

  return_status = interpretCallStatus(options_.log_options,
                                      applyScalingToLpRow(lp, row, scale_value),
                                      return_status, "applyScalingToLpRow");
  if (return_status == HighsStatus::kError) return return_status;

  if (scale_value < 0 && basis.valid) {
    // Negative, so flip any nonbasic status
    if (basis.row_status[row] == HighsBasisStatus::kLower) {
      basis.row_status[row] = HighsBasisStatus::kUpper;
    } else if (basis.row_status[row] == HighsBasisStatus::kUpper) {
      basis.row_status[row] = HighsBasisStatus::kLower;
    }
  }
  if (simplex_status.initialised_for_solve) {
    SimplexBasis& simplex_basis = ekk_instance_.basis_;
    if (scale_value < 0 && simplex_status.has_basis) {
      // Negative, so flip any nonbasic status
      const HighsInt var = lp.num_col_ + row;
      if (simplex_basis.nonbasicMove_[var] == kNonbasicMoveUp) {
        simplex_basis.nonbasicMove_[var] = kNonbasicMoveDn;
      } else if (simplex_basis.nonbasicMove_[var] == kNonbasicMoveDn) {
        simplex_basis.nonbasicMove_[var] = kNonbasicMoveUp;
      }
    }
  }
  // Deduce the consequences of a scaled row
  invalidateModelStatusSolutionAndInfo();

  // Determine any implications for simplex data
  ekk_instance_.updateStatus(LpAction::kScaledRow);
  return HighsStatus::kOk;
}

void Highs::setNonbasicStatusInterface(
    const HighsIndexCollection& index_collection, const bool columns) {
  HighsBasis& highs_basis = basis_;
  if (!highs_basis.valid) return;
  const bool has_simplex_basis = ekk_instance_.status_.has_basis;
  SimplexBasis& simplex_basis = ekk_instance_.basis_;
  HighsLp& lp = model_.lp_;

  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  HighsInt ix_dim;
  if (columns) {
    ix_dim = lp.num_col_;
  } else {
    ix_dim = lp.num_row_;
  }
  // Surely this is checked elsewhere
  assert(0 <= from_k && to_k < ix_dim);
  assert(from_k <= to_k);
  HighsInt set_from_ix;
  HighsInt set_to_ix;
  HighsInt ignore_from_ix;
  HighsInt ignore_to_ix = -1;
  HighsInt current_set_entry = 0;
  // Given a basic-nonbasic partition, all status settings are defined
  // by the bounds unless boxed, in which case any definitive (ie not
  // just kNonbasic) existing status is retained. Otherwise, set to
  // bound nearer to zero
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, set_from_ix, set_to_ix, ignore_from_ix,
                     ignore_to_ix, current_set_entry);
    assert(set_to_ix < ix_dim);
    assert(ignore_to_ix < ix_dim);
    if (columns) {
      for (HighsInt iCol = set_from_ix; iCol <= set_to_ix; iCol++) {
        if (highs_basis.col_status[iCol] == HighsBasisStatus::kBasic) continue;
        // Nonbasic column
        double lower = lp.col_lower_[iCol];
        double upper = lp.col_upper_[iCol];
        HighsBasisStatus status = highs_basis.col_status[iCol];
        int8_t move = kIllegalMoveValue;
        if (lower == upper) {
          if (status == HighsBasisStatus::kNonbasic)
            status = HighsBasisStatus::kLower;
          move = kNonbasicMoveZe;
        } else if (!highs_isInfinity(-lower)) {
          // Finite lower bound so boxed or lower
          if (!highs_isInfinity(upper)) {
            // Finite upper bound so boxed
            if (status == HighsBasisStatus::kNonbasic) {
              // No definitive status, so set to bound nearer to zero
              if (fabs(lower) < fabs(upper)) {
                status = HighsBasisStatus::kLower;
                move = kNonbasicMoveUp;
              } else {
                status = HighsBasisStatus::kUpper;
                move = kNonbasicMoveDn;
              }
            } else if (status == HighsBasisStatus::kLower) {
              move = kNonbasicMoveUp;
            } else {
              move = kNonbasicMoveDn;
            }
          } else {
            // Lower (since upper bound is infinite)
            status = HighsBasisStatus::kLower;
            move = kNonbasicMoveUp;
          }
        } else if (!highs_isInfinity(upper)) {
          // Upper
          status = HighsBasisStatus::kUpper;
          move = kNonbasicMoveDn;
        } else {
          // FREE
          status = HighsBasisStatus::kZero;
          move = kNonbasicMoveZe;
        }
        highs_basis.col_status[iCol] = status;
        if (has_simplex_basis) {
          assert(move != kIllegalMoveValue);
          simplex_basis.nonbasicFlag_[iCol] = kNonbasicFlagTrue;
          simplex_basis.nonbasicMove_[iCol] = move;
        }
      }
    } else {
      for (HighsInt iRow = set_from_ix; iRow <= set_to_ix; iRow++) {
        if (highs_basis.row_status[iRow] == HighsBasisStatus::kBasic) continue;
        // Nonbasic column
        double lower = lp.row_lower_[iRow];
        double upper = lp.row_upper_[iRow];
        HighsBasisStatus status = highs_basis.row_status[iRow];
        int8_t move = kIllegalMoveValue;
        if (lower == upper) {
          if (status == HighsBasisStatus::kNonbasic)
            status = HighsBasisStatus::kLower;
          move = kNonbasicMoveZe;
        } else if (!highs_isInfinity(-lower)) {
          // Finite lower bound so boxed or lower
          if (!highs_isInfinity(upper)) {
            // Finite upper bound so boxed
            if (status == HighsBasisStatus::kNonbasic) {
              // No definitive status, so set to bound nearer to zero
              if (fabs(lower) < fabs(upper)) {
                status = HighsBasisStatus::kLower;
                move = kNonbasicMoveDn;
              } else {
                status = HighsBasisStatus::kUpper;
                move = kNonbasicMoveUp;
              }
            } else if (status == HighsBasisStatus::kLower) {
              move = kNonbasicMoveDn;
            } else {
              move = kNonbasicMoveUp;
            }
          } else {
            // Lower (since upper bound is infinite)
            status = HighsBasisStatus::kLower;
            move = kNonbasicMoveDn;
          }
        } else if (!highs_isInfinity(upper)) {
          // Upper
          status = HighsBasisStatus::kUpper;
          move = kNonbasicMoveUp;
        } else {
          // FREE
          status = HighsBasisStatus::kZero;
          move = kNonbasicMoveZe;
        }
        highs_basis.row_status[iRow] = status;
        if (has_simplex_basis) {
          assert(move != kIllegalMoveValue);
          simplex_basis.nonbasicFlag_[lp.num_col_ + iRow] = kNonbasicFlagTrue;
          simplex_basis.nonbasicMove_[lp.num_col_ + iRow] = move;
        }
      }
    }
    if (ignore_to_ix >= ix_dim - 1) break;
  }
}

void Highs::appendNonbasicColsToBasisInterface(const HighsInt ext_num_new_col) {
  if (ext_num_new_col == 0) return;
  HighsBasis& highs_basis = basis_;
  if (!highs_basis.useful) return;
  const bool has_simplex_basis = ekk_instance_.status_.has_basis;
  SimplexBasis& simplex_basis = ekk_instance_.basis_;
  HighsLp& lp = model_.lp_;

  assert(highs_basis.col_status.size() == static_cast<size_t>(lp.num_col_));
  assert(highs_basis.row_status.size() == static_cast<size_t>(lp.num_row_));

  // Add nonbasic structurals
  HighsInt newNumCol = lp.num_col_ + ext_num_new_col;
  HighsInt newNumTot = newNumCol + lp.num_row_;
  highs_basis.col_status.resize(newNumCol);
  if (has_simplex_basis) {
    simplex_basis.nonbasicFlag_.resize(newNumTot);
    simplex_basis.nonbasicMove_.resize(newNumTot);
    // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
    // necessary
    for (HighsInt iRow = lp.num_row_ - 1; iRow >= 0; iRow--) {
      HighsInt iCol = simplex_basis.basicIndex_[iRow];
      if (iCol >= lp.num_col_) {
        // This basic variable is a row, so shift its index
        simplex_basis.basicIndex_[iRow] += ext_num_new_col;
      }
      simplex_basis.nonbasicFlag_[newNumCol + iRow] =
          simplex_basis.nonbasicFlag_[lp.num_col_ + iRow];
      simplex_basis.nonbasicMove_[newNumCol + iRow] =
          simplex_basis.nonbasicMove_[lp.num_col_ + iRow];
    }
  }
  // Make any new columns nonbasic
  for (HighsInt iCol = lp.num_col_; iCol < newNumCol; iCol++) {
    double lower = lp.col_lower_[iCol];
    double upper = lp.col_upper_[iCol];
    HighsBasisStatus status = HighsBasisStatus::kNonbasic;
    int8_t move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      status = HighsBasisStatus::kLower;
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (fabs(lower) < fabs(upper)) {
          status = HighsBasisStatus::kLower;
          move = kNonbasicMoveUp;
        } else {
          status = HighsBasisStatus::kUpper;
          move = kNonbasicMoveDn;
        }
      } else {
        // Lower (since upper bound is infinite)
        status = HighsBasisStatus::kLower;
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      status = HighsBasisStatus::kUpper;
      move = kNonbasicMoveDn;
    } else {
      // FREE
      status = HighsBasisStatus::kZero;
      move = kNonbasicMoveZe;
    }
    assert(status != HighsBasisStatus::kNonbasic);
    highs_basis.col_status[iCol] = status;
    if (has_simplex_basis) {
      assert(move != kIllegalMoveValue);
      simplex_basis.nonbasicFlag_[iCol] = kNonbasicFlagTrue;
      simplex_basis.nonbasicMove_[iCol] = move;
    }
  }
}

void Highs::appendBasicRowsToBasisInterface(const HighsInt ext_num_new_row) {
  if (ext_num_new_row == 0) return;
  HighsBasis& highs_basis = basis_;
  if (!highs_basis.useful) return;
  const bool has_simplex_basis = ekk_instance_.status_.has_basis;
  SimplexBasis& simplex_basis = ekk_instance_.basis_;
  HighsLp& lp = model_.lp_;

  assert(highs_basis.col_status.size() == static_cast<size_t>(lp.num_col_));
  assert(highs_basis.row_status.size() == static_cast<size_t>(lp.num_row_));

  // Add basic logicals
  // Add the new rows to the Highs basis
  HighsInt newNumRow = lp.num_row_ + ext_num_new_row;
  highs_basis.row_status.resize(newNumRow);
  for (HighsInt iRow = lp.num_row_; iRow < newNumRow; iRow++)
    highs_basis.row_status[iRow] = HighsBasisStatus::kBasic;
  if (has_simplex_basis) {
    // Add the new rows to the simplex basis
    HighsInt newNumTot = lp.num_col_ + newNumRow;
    simplex_basis.nonbasicFlag_.resize(newNumTot);
    simplex_basis.nonbasicMove_.resize(newNumTot);
    simplex_basis.basicIndex_.resize(newNumRow);
    for (HighsInt iRow = lp.num_row_; iRow < newNumRow; iRow++) {
      simplex_basis.nonbasicFlag_[lp.num_col_ + iRow] = kNonbasicFlagFalse;
      simplex_basis.nonbasicMove_[lp.num_col_ + iRow] = 0;
      simplex_basis.basicIndex_[iRow] = lp.num_col_ + iRow;
    }
  }
}

// Get the basic variables, performing INVERT if necessary
HighsStatus Highs::getBasicVariablesInterface(HighsInt* basic_variables) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsInt num_row = lp.num_row_;
  HighsInt num_col = lp.num_col_;
  HighsSimplexStatus& ekk_status = ekk_instance_.status_;
  // For an LP with no rows the solution is vacuous
  if (num_row == 0) return return_status;
  if (!basis_.valid) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasicVariables called without a HiGHS basis\n");
    return HighsStatus::kError;
  }
  if (!ekk_status.has_invert) {
    // The LP has no invert to use, so have to set one up, but only
    // for the current basis, so return_value is the rank deficiency.
    HighsLpSolverObject solver_object(lp, basis_, solution_, info_,
                                      ekk_instance_, callback_, options_,
                                      timer_);
    solver_object.setProfiling(this->profiling_);
    const bool only_from_known_basis = true;
    return_status = interpretCallStatus(
        options_.log_options,
        formSimplexLpBasisAndFactor(solver_object, only_from_known_basis),
        return_status, "formSimplexLpBasisAndFactor");
    if (return_status != HighsStatus::kOk) return return_status;
  }
  assert(ekk_status.has_invert);

  for (HighsInt row = 0; row < num_row; row++) {
    HighsInt var = ekk_instance_.basis_.basicIndex_[row];
    if (var < num_col) {
      basic_variables[row] = var;
    } else {
      basic_variables[row] = -(1 + var - num_col);
    }
  }
  return return_status;
}

// Solve (transposed) system involving the basis matrix

HighsStatus Highs::basisSolveInterface(const vector<double>& rhs,
                                       double* solution_vector,
                                       HighsInt* solution_num_nz,
                                       HighsInt* solution_indices,
                                       bool transpose) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsInt num_row = lp.num_row_;
  // For an LP with no rows the solution is vacuous
  if (num_row == 0) return return_status;
  // EKK must have an INVERT, but simplex NLA may need the pointer to
  // its LP to be refreshed so that it can use its scale factors
  assert(ekk_instance_.status_.has_invert);
  // Reset the simplex NLA LP and scale pointers for the unscaled LP
  ekk_instance_.setNlaPointersForLpAndScale(lp);
  assert(!lp.is_moved_);
  // Set up solve vector with suitably scaled RHS
  HVector solve_vector;
  solve_vector.setup(num_row);
  solve_vector.clear();
  HighsInt rhs_num_nz = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (rhs[iRow]) {
      solve_vector.index[rhs_num_nz++] = iRow;
      solve_vector.array[iRow] = rhs[iRow];
    }
  }
  solve_vector.count = rhs_num_nz;
  //
  // Note that solve_vector.count is just used to determine whether
  // hyper-sparse solves should be used. The indices of the nonzeros
  // in the solution are always accumulated. There's no switch (such
  // as setting solve_vector.count = num_row+1) to not do this.
  //
  // Get expected_density from analysis during simplex solve.
  const double expected_density = 1;
  if (transpose) {
    ekk_instance_.btran(solve_vector, expected_density);
  } else {
    ekk_instance_.ftran(solve_vector, expected_density);
  }
  // Extract the solution
  if (solution_indices == NULL) {
    // Nonzeros in the solution not required
    if (solve_vector.count > num_row) {
      // Solution nonzeros not known
      for (HighsInt iRow = 0; iRow < num_row; iRow++) {
        solution_vector[iRow] = solve_vector.array[iRow];
      }
    } else {
      // Solution nonzeros are known
      for (HighsInt iRow = 0; iRow < num_row; iRow++) solution_vector[iRow] = 0;
      for (HighsInt iX = 0; iX < solve_vector.count; iX++) {
        HighsInt iRow = solve_vector.index[iX];
        solution_vector[iRow] = solve_vector.array[iRow];
      }
    }
  } else {
    // Nonzeros in the solution are required
    if (solve_vector.count > num_row) {
      // Solution nonzeros not known
      solution_num_nz = 0;
      for (HighsInt iRow = 0; iRow < num_row; iRow++) {
        solution_vector[iRow] = 0;
        if (solve_vector.array[iRow]) {
          solution_vector[iRow] = solve_vector.array[iRow];
          solution_indices[*solution_num_nz++] = iRow;
        }
      }
    } else {
      // Solution nonzeros are known
      for (HighsInt iRow = 0; iRow < num_row; iRow++) solution_vector[iRow] = 0;
      for (HighsInt iX = 0; iX < solve_vector.count; iX++) {
        HighsInt iRow = solve_vector.index[iX];
        solution_vector[iRow] = solve_vector.array[iRow];
        solution_indices[iX] = iRow;
      }
      *solution_num_nz = solve_vector.count;
    }
  }
  return HighsStatus::kOk;
}

void Highs::zeroIterationCounts() {
  info_.simplex_iteration_count = 0;
  info_.ipm_iteration_count = 0;
  info_.crossover_iteration_count = 0;
  info_.pdlp_iteration_count = 0;
  info_.qp_iteration_count = 0;
}

HighsStatus Highs::getDualRayInterface(bool& has_dual_ray,
                                       double* dual_ray_value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsInt num_row = lp.num_row_;
  // For an LP with no rows the dual ray is vacuous
  if (num_row == 0) return return_status;
  bool has_invert = ekk_instance_.status_.has_invert;
  assert(!lp.is_moved_);
  has_dual_ray = ekk_instance_.dual_ray_record_.index != kNoRayIndex;

  // Declare identifiers to save column costs, integrality, any
  // Hessian and the presolve setting, and a flag to know when they
  // should be recovered
  std::vector<double> col_cost;
  HighsHessian hessian;
  bool solve_relaxation;
  std::string presolve;
  std::string solver;
  bool allow_unbounded_or_infeasible;
  bool solve_feasibility_problem = false;
  const bool is_qp = model_.isQp();

  if (dual_ray_value != NULL) {
    // User wants a dual ray whatever
    if (!has_dual_ray || !has_invert) {
      // No dual ray is known, or no INVERT to compute it
      //
      // No point in trying to get a dual ray if the model status is
      // optimal
      if (this->model_status_ == HighsModelStatus::kOptimal) {
        highsLogUser(options_.log_options, HighsLogType::kInfo,
                     "Model status is optimal, so no dual ray is available\n");
        return return_status;
      }
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Solving LP to try to compute dual ray\n");
      // Save the column costs, any Hessian, the relaxation, presolve
      // and solver settings, so that the feasibility problem is
      // solved to yield a basic solution
      col_cost = lp.col_cost_;
      if (is_qp) hessian = model_.hessian_;
      solve_relaxation = this->options_.solve_relaxation;
      presolve = this->options_.presolve;
      solver = this->options_.solver;
      allow_unbounded_or_infeasible =
          this->options_.allow_unbounded_or_infeasible;
      solve_feasibility_problem = true;
      // Zero the costs and Hessian
      std::vector<double> zero_costs;
      zero_costs.assign(lp.num_col_, 0);
      // Take a copy of the primal ray record, since this will be
      // cleared by calling changeColsCost
      HighsRayRecord primal_ray_record =
          this->ekk_instance_.primal_ray_record_.getRayRecord();
      HighsStatus status =
          this->changeColsCost(0, lp.num_col_ - 1, zero_costs.data());
      assert(status == HighsStatus::kOk);
      // Reinstate the primal ray record
      this->ekk_instance_.primal_ray_record_.setRayRecord(primal_ray_record);
      if (is_qp) {
        HighsHessian zero_hessian;
        this->passHessian(zero_hessian);
      }
      //
      this->options_.solve_relaxation = true;
      this->options_.presolve = kHighsOffString;
      this->options_.solver = kSimplexString;
      this->options_.allow_unbounded_or_infeasible = false;
      HighsStatus call_status = this->run();
      // Restore the column costs, any Hessian, the relaxation,
      // presolve and solver settings
      this->options_.solve_relaxation = solve_relaxation;
      this->options_.presolve = presolve;
      this->options_.solver = solver;
      this->options_.allow_unbounded_or_infeasible =
          allow_unbounded_or_infeasible;

      lp.col_cost_ = col_cost;
      if (is_qp) model_.hessian_ = hessian;
      if (call_status != HighsStatus::kOk) return_status = call_status;
      has_dual_ray = ekk_instance_.dual_ray_record_.index != kNoRayIndex;
      has_invert = ekk_instance_.status_.has_invert;
      assert(has_invert);
    }
    if (has_dual_ray) {
      if (ekk_instance_.dual_ray_record_.value.size()) {
        // Dual ray is already computed
        highsLogUser(options_.log_options, HighsLogType::kInfo,
                     "Copying known dual ray\n");
        for (HighsInt iRow = 0; iRow < num_row; iRow++)
          dual_ray_value[iRow] = ekk_instance_.dual_ray_record_.value[iRow];
      } else if (has_invert) {
        // Dual ray is known and can be calculated
        highsLogUser(options_.log_options, HighsLogType::kInfo,
                     "Solving linear system to compute dual ray\n");
        vector<double> rhs;
        HighsInt iRow = ekk_instance_.dual_ray_record_.index;
        rhs.assign(num_row, 0);
        rhs[iRow] = ekk_instance_.dual_ray_record_.sign;
        HighsInt* dual_ray_num_nz = 0;
        basisSolveInterface(rhs, dual_ray_value, dual_ray_num_nz, NULL, true);
        // Now save the dual ray itself
        ekk_instance_.dual_ray_record_.value.resize(num_row);
        for (HighsInt iRow = 0; iRow < num_row; iRow++)
          ekk_instance_.dual_ray_record_.value[iRow] = dual_ray_value[iRow];
      } else {
        assert(!has_invert);
        // Dual ray is known but cannot be calculated
        highsLogUser(options_.log_options, HighsLogType::kError,
                     "No LP invertible representation to compute dual ray\n");
        return_status = HighsStatus::kError;
      }
    } else {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "No dual ray found\n");
      return_status = HighsStatus::kOk;
    }
  }
  if (solve_feasibility_problem) {
    // Feasibility problem has been solved, so any objective-related
    // information has been lost. Reverting the objective function via
    // Highs calls clears info_, so better to just copy the data
    // directly and set the info_ entries that are no longer valid
    //
    // The relaxation for an infeasible MIP may be feasible - so no
    // ray is generated - so make sure (#2415) that the primal
    // solution status is reset
    this->info_.primal_solution_status = SolutionStatus::kSolutionStatusNone;
    // Modify the objective-related information
    this->info_.dual_solution_status = SolutionStatus::kSolutionStatusNone;
    this->info_.objective_function_value = 0;
    this->info_.invalidateDualKkt();
    if (has_dual_ray) {
      assert(this->info_.num_primal_infeasibilities > 0 ||
             this->info_.num_semi_infeasibilities > 0);
      assert(this->model_status_ == HighsModelStatus::kInfeasible);
    } else {
      // If someone has tried to get a dual ray for a feasible problem
      // - or if the relaxation is feasible - then any model and
      // primal KKT status of the original model has been lost
      this->info_.invalidatePrimalKkt();
      this->model_status_ = HighsModelStatus::kNotset;
    }
  }
  return return_status;
}

HighsStatus Highs::getPrimalRayInterface(bool& has_primal_ray,
                                         double* primal_ray_value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsInt num_row = lp.num_row_;
  HighsInt num_col = lp.num_col_;
  // For an LP with no rows the primal ray is vacuous
  if (num_row == 0) return return_status;
  if (model_.isQp()) {
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Cannot find primal ray for unbounded QP\n");
    return HighsStatus::kError;
  }
  bool has_invert = ekk_instance_.status_.has_invert;
  assert(!lp.is_moved_);
  has_primal_ray = ekk_instance_.primal_ray_record_.index != kNoRayIndex;

  bool solve_relaxation;
  std::string presolve;
  std::string solver;
  bool allow_unbounded_or_infeasible;
  bool solve_unboundedness_problem = false;

  if (primal_ray_value != NULL) {
    // User wants a primal ray whatever
    if (!has_primal_ray || !has_invert) {
      // No primal ray is known, or no INVERT to compute it
      //
      // No point in trying to get a primal ray if the model status is
      // optimal
      if (this->model_status_ == HighsModelStatus::kOptimal) {
        highsLogUser(
            options_.log_options, HighsLogType::kInfo,
            "Model status is optimal, so no primal ray is available\n");
        return return_status;
      }
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Solving LP to try to compute primal ray\n");
      // Save the relaxation, presolve and solver settings
      solve_relaxation = this->options_.solve_relaxation;
      presolve = this->options_.presolve;
      solver = this->options_.solver;
      allow_unbounded_or_infeasible =
          this->options_.allow_unbounded_or_infeasible;
      solve_unboundedness_problem = true;
      this->options_.solve_relaxation = true;
      this->options_.presolve = kHighsOffString;
      this->options_.solver = kSimplexString;
      this->options_.allow_unbounded_or_infeasible = false;
      HighsStatus call_status = this->run();
      // Restore the relaxation, presolve and solver settings
      this->options_.solve_relaxation = solve_relaxation;
      this->options_.presolve = presolve;
      this->options_.solver = solver;
      this->options_.allow_unbounded_or_infeasible =
          allow_unbounded_or_infeasible;
      if (call_status != HighsStatus::kOk) return_status = call_status;
      has_primal_ray = ekk_instance_.primal_ray_record_.index != kNoRayIndex;
      has_invert = ekk_instance_.status_.has_invert;
      assert(has_invert);
    }
    if (has_primal_ray) {
      if (ekk_instance_.primal_ray_record_.value.size()) {
        // Primal ray is already computed
        highsLogUser(options_.log_options, HighsLogType::kInfo,
                     "Copying known primal ray\n");
        for (HighsInt iCol = 0; iCol < num_col; iCol++)
          primal_ray_value[iCol] = ekk_instance_.primal_ray_record_.value[iCol];
        return return_status;
      } else if (has_invert) {
        // Primal ray is known and can be calculated
        highsLogUser(options_.log_options, HighsLogType::kInfo,
                     "Solving linear system to compute primal ray\n");
        HighsInt col = ekk_instance_.primal_ray_record_.index;
        assert(ekk_instance_.basis_.nonbasicFlag_[col] == kNonbasicFlagTrue);
        // Get this pivotal column
        vector<double> rhs;
        vector<double> column;
        column.assign(num_row, 0);
        rhs.assign(num_row, 0);
        lp.ensureColwise();
        HighsInt primal_ray_sign = ekk_instance_.primal_ray_record_.sign;
        if (col < num_col) {
          for (HighsInt iEl = lp.a_matrix_.start_[col];
               iEl < lp.a_matrix_.start_[col + 1]; iEl++)
            rhs[lp.a_matrix_.index_[iEl]] =
                primal_ray_sign * lp.a_matrix_.value_[iEl];
        } else {
          rhs[col - num_col] = primal_ray_sign;
        }
        HighsInt* column_num_nz = 0;
        basisSolveInterface(rhs, column.data(), column_num_nz, NULL, false);
        // Now zero primal_ray_value and scatter the column according to
        // the basic variables.
        for (HighsInt iCol = 0; iCol < num_col; iCol++)
          primal_ray_value[iCol] = 0;
        for (HighsInt iRow = 0; iRow < num_row; iRow++) {
          HighsInt iCol = ekk_instance_.basis_.basicIndex_[iRow];
          if (iCol < num_col) primal_ray_value[iCol] = column[iRow];
        }
        if (col < num_col) primal_ray_value[col] = -primal_ray_sign;
        // Now save the primal ray itself
        ekk_instance_.primal_ray_record_.value.resize(num_col);
        for (HighsInt iCol = 0; iCol < num_col; iCol++)
          ekk_instance_.primal_ray_record_.value[iCol] = primal_ray_value[iCol];
      }
    } else {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "No primal ray found\n");
      return_status = HighsStatus::kOk;
    }
  }
  const bool is_mip = this->model_.isMip();
  if (solve_unboundedness_problem) {
    if (is_mip) {
      // Unboundedness LP has been solved, but that will give dual
      // solution status kInfeasible which, for a MIP is not correct
      this->info_.dual_solution_status = SolutionStatus::kSolutionStatusNone;
      this->info_.invalidateDualKkt();
    }
    if (has_primal_ray) {
      assert(is_mip || this->info_.num_dual_infeasibilities > 0);
      assert(this->model_status_ == HighsModelStatus::kUnbounded);
    }
  }
  return return_status;
}

HighsStatus Highs::getRangingInterface() {
  HighsLpSolverObject solver_object(model_.lp_, basis_, solution_, info_,
                                    ekk_instance_, callback_, options_, timer_);
  solver_object.setProfiling(this->profiling_);
  solver_object.model_status_ = model_status_;
  return getRangingData(this->ranging_, solver_object);
}

HighsStatus Highs::getIisInterfaceReturn(
    const HighsStatus return_status, const HighsOptions& original_options,
    const std::vector<HighsBool>& original_callback_active) {
  // Restore options and callbacks
  this->options_ = original_options;
  for (int i = kCallbackMin; i <= kCallbackMax; i++) {
    if (original_callback_active[i]) this->startCallback(i);
  }

  // Exit early if there was an error
  if (return_status == HighsStatus::kError) {
    // Report final results
    this->iis_.reportFinal(original_options);
    return return_status;
  }
  assert(this->iis_.valid_);

  // A valid HighsIis instance is one for which the information is
  // known to be correct
  //
  // In the case of a system that's feasible this will give empty
  // HighsIis::col_index_ and HighsIis::row_index_
  //
  // If the system is infeasible, then it's possible that only an IS
  // (not an IIS) has been formed. This will be identified by the
  // HighsIis::col_bound_ and HighsIis::row_bound_ of the columns and
  // rows in HighsIis::col_index_ and HighsIis::row_index_ being
  // kIisStatusMaybeInConflict, rather than kIisStatusInConflict
  //
  // The columns and rows not in HighsIis::col_index_ and
  // HighsIis::row_index_ have HighsIis::col_bound_ and
  // HighsIis::row_bound_ values set to kIisStatusNotInConflict

  // If the IIS process has identified infeasibility, then set
  if (this->iis_.status_ >= kIisModelStatusTimeLimit)
    this->model_status_ = HighsModelStatus::kInfeasible;

  HighsLp& lp = this->model_.lp_;
  HighsLp& iis_lp = this->iis_.model_.lp_;
  const bool has_is =
      this->iis_.col_index_.size() || this->iis_.row_index_.size();
  const bool has_iis = this->iis_.status_ == kIisModelStatusIrreducible;

  HighsOptions opts = this->options_;
  opts.time_limit = kHighsInf;
  opts.simplex_iteration_limit = kHighsIInf;
  opts.objective_bound = kHighsInf;

  if (has_iis) assert(has_is);
  if (has_is) {
    // Construct the HighsIis LP
    this->iis_.setLp(lp);
    // Check that the LP data are OK (correspond to original model
    // reduced to HighsIis col/row and bound data).
    bool lp_data_ok = this->iis_.lpDataOk(lp, opts);
    assert(lp_data_ok);
    if (!lp_data_ok) {
      // Report final results
      this->iis_.reportFinal(original_options);
      return HighsStatus::kError;
    }
    // Check that the HighsIis LP has the property of being infeasible
    // and, if a true IIS is claimed, optimal/unbounded if any bound
    // is relaxed
    bool lp_ok = this->iis_.lpOk(opts);
    if (!lp_ok) {
      // If we fail to prove infeasibility mark candidate set as invalid and
      // return with a warning
      this->iis_.valid_ = false;
      if (this->iis_.status_ != kIisModelStatusTimeLimit)
        this->iis_.status_ = kIisModelStatusUnknown;
      // Report final results
      this->iis_.reportFinal(original_options);
      return HighsStatus::kWarning;
    }
  } else {
    assert(this->iis_.status_ <= kIisModelStatusReducible);
    assert(!this->iis_.col_bound_.size());
    assert(!this->iis_.row_bound_.size());
  }
  // Construct the ISS status vectors for cols and rows of original
  // model
  this->iis_.setStatus(lp);
  // Check consistency of the col/row_index_ and col/row_status_
  bool index_status_ok = this->iis_.indexStatusOk(lp);
  assert(index_status_ok);
  if (!index_status_ok) {
    this->iis_.reportFinal(original_options);
    return HighsStatus::kError;
  }

  // Report final results
  this->iis_.reportFinal(original_options);
  return return_status;
}

HighsStatus Highs::getIisInterface() {
  if (this->model_status_ == HighsModelStatus::kOptimal ||
      this->model_status_ == HighsModelStatus::kUnbounded) {
    // Strange to call getIis for a model that's known to be feasible
    highsLogUser(
        options_.log_options, HighsLogType::kInfo,
        "Calling Highs::getIis for a model that is known to be feasible\n");
    this->iis_.clear();
    // No IIS exists, so validate the empty HighsIis instance
    this->iis_.valid_ = true;
    this->iis_.status_ = kIisModelStatusFeasible;
    return this->getIisInterfaceReturn(HighsStatus::kOk, options_,
                                       callback_.active);
  }
  // Early exit for existing valid IIS
  if (this->iis_.valid_)
    return this->getIisInterfaceReturn(HighsStatus::kOk, options_,
                                       callback_.active);
  // Clear IIS
  this->iis_.clear();
  // Check for trivial IIS: empty infeasible row or inconsistent bounds
  const HighsLp& lp = model_.lp_;
  if (this->iis_.trivial(lp, options_))
    return this->getIisInterfaceReturn(HighsStatus::kOk, options_,
                                       callback_.active);
  HighsInt num_row = lp.num_row_;
  if (num_row == 0) {
    // For an LP with no rows, the only scope for infeasibility is
    // inconsistent columns bounds - which has already been assessed,
    // so validate the empty HighsIis instance
    this->iis_.valid_ = true;
    return this->getIisInterfaceReturn(HighsStatus::kOk, options_,
                                       callback_.active);
  }
  // Look for infeasible rows based on row value bounds
  if (this->iis_.rowValueBounds(lp, options_))
    return this->getIisInterfaceReturn(HighsStatus::kOk, options_,
                                       callback_.active);
  // Don't continue with more expensive techniques if using the IIS
  // light strategy
  if (options_.iis_strategy == kIisStrategyLight)
    return this->getIisInterfaceReturn(HighsStatus::kOk, options_,
                                       callback_.active);
  // Clear IIS
  this->iis_.clear();
  // Save original options
  HighsOptions original_options = this->options_;
  // Save original active callbacks and disable all except for
  // kCallbackLogging and kCallbackSimplexInterrupt
  std::vector<HighsBool> original_callback_active = callback_.active;
  for (int i = kCallbackMin; i <= kCallbackMax; i++) {
    if (i != kCallbackLogging && i != kCallbackSimplexInterrupt &&
        callback_.active[i])
      this->stopCallback(i);
  }
  // Zero out all clocks and set time_limit to iis_time_limit
  this->zeroAllClocks();
  this->setOptionValue("time_limit", options_.iis_time_limit);
  // Disallow kIterationLimit and kObjectiveBound model statuses
  this->setOptionValue("simplex_iteration_limit", kHighsIInf);
  this->setOptionValue("objective_bound", kHighsInf);
  // Disallow kUnboundedOrInfeasible model status
  this->setOptionValue("allow_unbounded_or_infeasible", false);
  // Run model with new options if needed
  if (this->model_status_ == HighsModelStatus::kNotset ||
      this->model_status_ == HighsModelStatus::kIterationLimit ||
      this->model_status_ == HighsModelStatus::kObjectiveBound ||
      this->model_status_ == HighsModelStatus::kInterrupt ||
      this->model_status_ == HighsModelStatus::kTimeLimit ||
      this->model_status_ == HighsModelStatus::kUnboundedOrInfeasible) {
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Model status is %s. Resolving to establish infeasibility "
                 "before computing IIS\n\n",
                 this->modelStatusToString(this->model_status_).c_str());
    double simplex_time = -this->getRunTime();
    HighsInt simplex_iterations = -info_.simplex_iteration_count;
    HighsStatus run_status = this->optimizeModel();
    simplex_time += this->getRunTime();
    simplex_iterations += info_.simplex_iteration_count;
    this->iis_.info_.update(simplex_time, simplex_iterations);
    highsLogUser(options_.log_options, HighsLogType::kInfo, "\n");
  }
  if (this->model_status_ == HighsModelStatus::kOptimal ||
      this->model_status_ == HighsModelStatus::kUnbounded) {
    // Model became feasible
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Model became feasible\n");
    // No IIS exists, so validate the empty HighsIis instance
    this->iis_.valid_ = true;
    this->iis_.status_ = kIisModelStatusFeasible;
    this->iis_.strategy_ = options_.iis_strategy;
    return this->getIisInterfaceReturn(HighsStatus::kOk, original_options,
                                       original_callback_active);
  } else if (this->model_status_ == HighsModelStatus::kTimeLimit) {
    // Time limit reached
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Time limit reached prior to establishing infeasibility\n");
    this->iis_.status_ = kIisModelStatusTimeLimit;
    this->iis_.strategy_ = options_.iis_strategy;
    return this->getIisInterfaceReturn(HighsStatus::kError, original_options,
                                       original_callback_active);
  } else if (this->model_status_ != HighsModelStatus::kInfeasible) {
    // kModelEmpty, kModelError, kSolveError, kMemoryLimit, kUnknown model
    // status
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Can not compute IIS for a model with status %s\n",
                 this->modelStatusToString(this->model_status_).c_str());
    this->iis_.strategy_ = options_.iis_strategy;
    return this->getIisInterfaceReturn(HighsStatus::kError, original_options,
                                       original_callback_active);
  }
  assert(this->model_status_ == HighsModelStatus::kInfeasible);

  // Construct an IS from ray or lp
  bool ray_option =
      // kIisStrategyFromRay & options_.iis_strategy;
      false;
  const bool lp_option = kIisStrategyFromLp & options_.iis_strategy;
  if (ray_option && !ekk_instance_.status_.has_invert) {
    // Model is known to be infeasible, and a dual ray option is
    // chosen, but it has no INVERT, presumably because infeasibility
    // detected in presolve, so solve without presolve
    std::string presolve = options_.presolve;
    options_.presolve = kHighsOffString;

    double simplex_time = -this->getRunTime();
    HighsInt simplex_iterations = -info_.simplex_iteration_count;
    HighsStatus run_status = this->optimizeModel();
    options_.presolve = presolve;
    simplex_time += this->getRunTime();
    simplex_iterations += info_.simplex_iteration_count;
    this->iis_.info_.update(simplex_time, simplex_iterations);

    // Model should remain infeasible!
    if (this->model_status_ == HighsModelStatus::kTimeLimit) {
      // Time limit reached
      this->iis_.status_ = kIisModelStatusTimeLimit;
      return this->getIisInterfaceReturn(HighsStatus::kError, original_options,
                                         original_callback_active);
    } else if (this->model_status_ != HighsModelStatus::kInfeasible) {
      highsLogUser(
          options_.log_options, HighsLogType::kError,
          "Model status has switched from %s to %s when solving without "
          "presolve\n",
          this->modelStatusToString(HighsModelStatus::kInfeasible).c_str(),
          this->modelStatusToString(this->model_status_).c_str());
      return this->getIisInterfaceReturn(HighsStatus::kError, original_options,
                                         original_callback_active);
    }
  }
  const bool has_dual_ray = ekk_instance_.dual_ray_record_.index != kNoRayIndex;
  if (ray_option && !has_dual_ray)
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "No known dual ray from which to compute IIS\n");
  ray_option = ray_option && has_dual_ray;
  HighsStatus return_status = HighsStatus::kOk;
  if (ray_option) {
    assert(has_dual_ray);
    // Compute the dual ray to identify an infeasible subset of rows
    assert(ekk_instance_.status_.has_invert);
    assert(!lp.is_moved_);
    std::vector<double> rhs;
    HighsInt iRow = ekk_instance_.dual_ray_record_.index;
    rhs.assign(num_row, 0);
    rhs[iRow] = 1;
    std::vector<double> dual_ray_value(num_row);
    HighsInt* dual_ray_num_nz = 0;
    basisSolveInterface(rhs, dual_ray_value.data(), dual_ray_num_nz, NULL,
                        true);
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
      if (dual_ray_value[iRow]) this->iis_.row_index_.push_back(iRow);
  } else if (lp_option) {
    // Full LP option chosen or no dual ray to use
    // Apply the elasticity filter to the whole model in order to
    // determine an infeasible subset of rows
    return_status = this->elasticityFilter(-1.0, -1.0, 1.0, nullptr, nullptr,
                                           nullptr, true);
  }
  // Do not continue if not using kIisStrategyIrreducible
  if (!(kIisStrategyIrreducible & this->options_.iis_strategy))
    return this->getIisInterfaceReturn(return_status, original_options,
                                       original_callback_active);
  // If neither ray and lp options were requested or if they fail to produce a
  // valid IS, make one consisting of all constraints
  if (!this->iis_.valid_ || iis_.status_ != kIisModelStatusReducible) {
    this->iis_.valid_ = true;
    this->iis_.status_ = kIisModelStatusReducible;
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
      this->iis_.row_index_.push_back(iRow);
  }

  assert(this->iis_.row_index_.size());
  // A subset of infeasible rows, so at least have a reducible
  // infeasibility set

  //  Attempt to compute a true IIS
  //
  //  To do this the matrix must be column-wise
  model_.lp_.a_matrix_.ensureColwise();
  return_status = this->iis_.deduce(lp, options_, callback_, basis_);

  // Analyse the LP solution data
  const HighsInt num_lp_solved = this->iis_.info_.num_lp_solved;
  const double min_time = this->iis_.info_.min_simplex_time;
  const double sum_time = this->iis_.info_.sum_simplex_times;
  const double max_time = this->iis_.info_.max_simplex_time;
  const HighsInt min_iterations = this->iis_.info_.min_simplex_iteration_count;
  const HighsInt sum_iterations = this->iis_.info_.sum_simplex_iteration_counts;
  const HighsInt max_iterations = this->iis_.info_.max_simplex_iteration_count;

  if (kIisDevReport) {
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 " %d cols, %d rows, %d LPs solved"
                 " (min / average / max) iteration count (%6d / %6.2g / % 6d)"
                 " and time (%6.2f / %6.2f / % 6.2f) \n",
                 int(this->iis_.col_index_.size()),
                 int(this->iis_.row_index_.size()), int(num_lp_solved),
                 int(min_iterations),
                 num_lp_solved > 0 ? (1.0 * sum_iterations) / num_lp_solved : 0,
                 int(max_iterations), min_time,
                 num_lp_solved > 0 ? sum_time / num_lp_solved : 0, max_time);
  }
  return this->getIisInterfaceReturn(return_status, original_options,
                                     original_callback_active);
}

HighsStatus Highs::elasticityFilterReturn(
    const HighsStatus return_status, const std::string& original_model_name,
    const HighsModelStatus original_model_status,
    const HighsInt original_num_col, const HighsInt original_num_row,
    const ObjSense original_sense, const std::vector<double>& original_col_cost,
    const std::vector<double>& original_col_lower,
    const std::vector<double>& original_col_upper,
    const std::vector<HighsVarType>& original_integrality) {
  const HighsLp& lp = this->model_.lp_;
  // The IIS is cleared by restoring the original LP, so save it
  HighsIis iis = this->iis_;
  double objective_function_value = info_.objective_function_value;
  // Delete any additional rows and columns, and restore the original
  // column costs and bounds
  HighsStatus run_status;
  run_status = this->deleteRows(original_num_row, lp.num_row_ - 1);
  assert(run_status == HighsStatus::kOk);

  run_status = this->deleteCols(original_num_col, lp.num_col_ - 1);
  assert(run_status == HighsStatus::kOk);
  //
  // Now that deleteRows and deleteCols may yield a valid basis, the
  // lack of dual values triggers an assert in
  // getKktFailures. Ultimately (#2081) the dual values will be
  // available but, for now, make the basis invalid.
  basis_.valid = false;

  run_status =
      this->changeColsCost(0, original_num_col - 1, original_col_cost.data());
  assert(run_status == HighsStatus::kOk);

  run_status = this->changeObjectiveSense(original_sense);
  assert(run_status == HighsStatus::kOk);

  run_status =
      this->changeColsBounds(0, original_num_col - 1, original_col_lower.data(),
                             original_col_upper.data());
  assert(run_status == HighsStatus::kOk);

  if (kIisStrategyRelaxation & this->options_.iis_strategy &&
      original_integrality.size()) {
    this->changeColsIntegrality(0, original_num_col - 1,
                                original_integrality.data());
    assert(run_status == HighsStatus::kOk);
  }
  // Revert the model name
  this->passModelName(original_model_name);
  assert(lp.num_col_ == original_num_col);
  assert(lp.num_row_ == original_num_row);

  // Clear solver data
  this->invalidateSolverData();

  if (return_status == HighsStatus::kOk) {
    // Solution is invalidated by deleting rows and columns, but
    // primal values are correct. Have to recompute row activities,
    // though
    this->model_.lp_.a_matrix_.productQuad(this->solution_.row_value,
                                           this->solution_.col_value);
    this->solution_.value_valid = true;
    // Set the feasibility objective and any KKT failures
    info_.objective_function_value = objective_function_value;
    getKktFailures(options_, model_, solution_, basis_, info_);
    info_.valid = true;
  }

  // Revert model status
  this->model_status_ = original_model_status;

  // Restore IIS
  this->iis_ = iis;

  return return_status;
}

HighsStatus Highs::elasticityFilter(const double global_lower_penalty,
                                    const double global_upper_penalty,
                                    const double global_rhs_penalty,
                                    const double* local_lower_penalty,
                                    const double* local_upper_penalty,
                                    const double* local_rhs_penalty,
                                    const bool get_iis) {
  //  this->writeModel("infeasible.mps");
  //
  // Solve the feasibility relaxation problem for the given penalties,
  // continuing to act as the elasticity filter if get_iis
  // is true, resulting in an infeasibility subset for further
  // refinement as an IIS
  //
  // Construct the e-LP:
  //
  // Constraints L <= Ax <= U; l <= x <= u
  //
  // Transformed to
  //
  // L <= Ax + e_L - e_U <= U,
  //
  // l <=  x + e_l - e_u <= u,
  //
  // where the elastic variables are not used if the corresponding
  // bound is infinite or the local/global penalty is negative.
  //
  // x is free, and the objective is the linear function of the
  // elastic variables given by the local/global penalties
  //
  // col_of_ecol lists the column indices corresponding to the entries in
  // bound_of_col_of_ecol_is_lower so that the results can be interpreted
  //
  // row_of_ecol lists the row indices corresponding to the entries in
  // bound_of_row_of_ecol_is_lower so that the results can be interpreted
  std::vector<HighsInt> col_of_ecol;
  std::vector<HighsInt> row_of_ecol;
  std::vector<HighsBool> bound_of_row_of_ecol_is_lower;
  std::vector<HighsBool> bound_of_col_of_ecol_is_lower;
  std::vector<double> erow_lower;
  std::vector<double> erow_upper;
  std::vector<HighsInt> erow_start;
  std::vector<HighsInt> erow_index;
  std::vector<double> erow_value;
  // Accumulate names for ecols and erows, re-using ecol_name for the
  // names of row ecols after defining the names of col ecols
  std::vector<std::string> ecol_name;
  std::vector<std::string> erow_name;
  std::vector<double> ecol_cost;
  std::vector<double> ecol_lower;
  std::vector<double> ecol_upper;
  const HighsLp& lp = this->model_.lp_;
  HighsInt evar_ix = lp.num_col_;
  HighsStatus run_status;
  const bool write_model = false;
  // Take copies of the original model name, dimensions and column
  // data vectors, as they will be modified in forming the e-LP
  const std::string original_model_name = lp.model_name_;
  const HighsModelStatus original_model_status = this->getModelStatus();
  const ObjSense original_sense = lp.sense_;
  const HighsInt original_num_col = lp.num_col_;
  const HighsInt original_num_row = lp.num_row_;
  const std::vector<double> original_col_cost = lp.col_cost_;
  const std::vector<double> original_col_lower = lp.col_lower_;
  const std::vector<double> original_col_upper = lp.col_upper_;
  const std::vector<HighsVarType> original_integrality = lp.integrality_;
  HighsIis iis = this->iis_;
  // Clear solver data
  this->invalidateSolverData();
  // Give the model a new name to avoid confusing logging when the
  // elastic LP is solved
  this->passModelName(original_model_name + "_elastic");
  // Zero the column costs
  std::vector<double> zero_costs;
  zero_costs.assign(original_num_col, 0);
  run_status = this->changeColsCost(0, lp.num_col_ - 1, zero_costs.data());
  assert(run_status == HighsStatus::kOk);
  // Change the sense to minimize
  run_status = this->changeObjectiveSense(ObjSense::kMinimize);
  assert(run_status == HighsStatus::kOk);

  const bool mip_relaxation_iis =
      get_iis && lp.isMip() &&
      kIisStrategyRelaxation & this->options_.iis_strategy;
  if (mip_relaxation_iis) {
    // Set any integrality to continuous
    std::vector<HighsVarType> all_continuous;
    all_continuous.assign(original_num_col, HighsVarType::kContinuous);
    run_status =
        this->changeColsIntegrality(0, lp.num_col_ - 1, all_continuous.data());
    assert(run_status == HighsStatus::kOk);
  }

  // For the columns
  const bool has_local_lower_penalty = local_lower_penalty;
  const bool has_global_elastic_lower = global_lower_penalty >= 0;
  const bool has_elastic_lower =
      has_local_lower_penalty || has_global_elastic_lower;
  const bool has_local_upper_penalty = local_upper_penalty;
  const bool has_global_elastic_upper = global_upper_penalty >= 0;
  const bool has_elastic_upper =
      has_local_upper_penalty || has_global_elastic_upper;
  const bool has_elastic_columns = has_elastic_lower || has_elastic_upper;
  // For the rows
  const bool has_local_rhs_penalty = local_rhs_penalty;
  const bool has_global_elastic_rhs = global_rhs_penalty >= 0;
  const bool has_elastic_rows = has_local_rhs_penalty || has_global_elastic_rhs;
  assert(has_elastic_columns || has_elastic_rows);

  const HighsInt col_ecol_offset = lp.num_col_;
  if (has_elastic_columns) {
    // Accumulate bounds to be used for columns
    std::vector<double> col_lower;
    std::vector<double> col_upper;
    // When defining names, need to know the column number
    const bool has_col_names = lp.col_names_.size() > 0;
    erow_start.push_back(0);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      const double lower = lp.col_lower_[iCol];
      const double upper = lp.col_upper_[iCol];
      // Original bounds used unless e-variable introduced
      col_lower.push_back(lower);
      col_upper.push_back(upper);
      // Free columns have no erow
      if (lower <= -kHighsInf && upper >= kHighsInf) continue;

      // Get the penalty for violating the lower bounds on this column
      const double lower_penalty = has_local_lower_penalty
                                       ? local_lower_penalty[iCol]
                                       : global_lower_penalty;
      // Negative lower penalty and infinite upper bound implies that the
      // bounds cannot be violated
      if (lower_penalty < 0 && upper >= kHighsInf) continue;

      // Get the penalty for violating the upper bounds on this column
      const double upper_penalty = has_local_upper_penalty
                                       ? local_upper_penalty[iCol]
                                       : global_upper_penalty;
      // Infinite upper bound and negative lower penalty implies that the
      // bounds cannot be violated
      if (lower <= -kHighsInf && upper_penalty < 0) continue;
      erow_lower.push_back(lower);
      erow_upper.push_back(upper);
      if (has_col_names)
        erow_name.push_back("row_" + std::to_string(iCol) + "_" +
                            lp.col_names_[iCol] + "_erow");
      // Define the entry for x[iCol]
      erow_index.push_back(iCol);
      erow_value.push_back(1);
      if (lower > -kHighsInf && lower_penalty >= 0) {
        // New e_l variable
        col_of_ecol.push_back(iCol);
        if (has_col_names)
          ecol_name.push_back("col_" + std::to_string(iCol) + "_" +
                              lp.col_names_[iCol] + "_lower");
        // Save the original lower bound on this column and free its
        // lower bound
        bound_of_col_of_ecol_is_lower.push_back(true);
        col_lower[iCol] = -kHighsInf;
        erow_index.push_back(evar_ix);
        erow_value.push_back(1);
        ecol_cost.push_back(lower_penalty);
        evar_ix++;
      }
      if (upper < kHighsInf && upper_penalty >= 0) {
        // New e_u variable
        col_of_ecol.push_back(iCol);
        if (has_col_names)
          ecol_name.push_back("col_" + std::to_string(iCol) + "_" +
                              lp.col_names_[iCol] + "_upper");
        // Save the original upper bound on this column and free its
        // upper bound
        bound_of_col_of_ecol_is_lower.push_back(false);
        col_upper[iCol] = kHighsInf;
        erow_index.push_back(evar_ix);
        erow_value.push_back(-1);
        ecol_cost.push_back(upper_penalty);
        evar_ix++;
      }
      erow_start.push_back(erow_index.size());
      HighsInt row_nz =
          erow_start[erow_start.size() - 1] - erow_start[erow_start.size() - 2];
      //      printf("eRow for column %d has %d nonzeros\n", int(iCol),
      //      int(row_nz));
      assert(row_nz == 2 || row_nz == 3);
    }
    HighsInt num_new_col = col_of_ecol.size();
    HighsInt num_new_row = erow_start.size() - 1;
    HighsInt num_new_nz = erow_start[num_new_row];
    if (kIisDevReport)
      printf(
          "Elasticity filter: For columns there are %d variables and %d "
          "constraints\n",
          int(num_new_col), int(num_new_row));
    // Apply the original column bound changes
    assert(col_lower.size() == static_cast<size_t>(lp.num_col_));
    assert(col_upper.size() == static_cast<size_t>(lp.num_col_));
    run_status = this->changeColsBounds(0, lp.num_col_ - 1, col_lower.data(),
                                        col_upper.data());
    assert(run_status == HighsStatus::kOk);
    // Add the new columns
    ecol_lower.assign(num_new_col, 0);
    ecol_upper.assign(num_new_col, kHighsInf);
    run_status = this->addCols(num_new_col, ecol_cost.data(), ecol_lower.data(),
                               ecol_upper.data(), 0, nullptr, nullptr, nullptr);
    assert(run_status == HighsStatus::kOk);
    // Add the new rows
    assert(erow_start.size() == static_cast<size_t>(num_new_row + 1));
    assert(erow_index.size() == static_cast<size_t>(num_new_nz));
    assert(erow_value.size() == static_cast<size_t>(num_new_nz));
    run_status = this->addRows(num_new_row, erow_lower.data(),
                               erow_upper.data(), num_new_nz, erow_start.data(),
                               erow_index.data(), erow_value.data());
    assert(run_status == HighsStatus::kOk);
    if (has_col_names) {
      for (HighsInt iCol = 0; iCol < num_new_col; iCol++)
        this->passColName(col_ecol_offset + iCol, ecol_name[iCol]);
      for (HighsInt iRow = 0; iRow < num_new_row; iRow++)
        this->passRowName(original_num_row + iRow, erow_name[iRow]);
    }
    assert(ecol_cost.size() == static_cast<size_t>(num_new_col));
    assert(ecol_lower.size() == static_cast<size_t>(num_new_col));
    assert(ecol_upper.size() == static_cast<size_t>(num_new_col));
    if (write_model) {
      printf("\nAfter adding %d e-rows\n=============\n", int(num_new_col));
      bool output_flag;
      run_status = this->getOptionValue("output_flag", output_flag);
      this->setOptionValue("output_flag", true);
      this->writeModel("");
      this->setOptionValue("output_flag", output_flag);
    }
  }
  const HighsInt row_ecol_offset = lp.num_col_;
  if (has_elastic_rows) {
    // Add the columns corresponding to the e_L and e_U variables for
    // the constraints
    ecol_name.clear();
    ecol_cost.clear();
    std::vector<HighsInt> ecol_start;
    std::vector<HighsInt> ecol_index;
    std::vector<double> ecol_value;
    ecol_start.push_back(0);
    const bool has_row_names = lp.row_names_.size() > 0;
    for (HighsInt iRow = 0; iRow < original_num_row; iRow++) {
      // Get the penalty for violating the bounds on this row
      const double penalty =
          has_local_rhs_penalty ? local_rhs_penalty[iRow] : global_rhs_penalty;
      // Negative penalty implies that the bounds cannot be violated
      if (penalty < 0) continue;
      const double lower = lp.row_lower_[iRow];
      const double upper = lp.row_upper_[iRow];
      if (lower > -kHighsInf) {
        // Create an e-var for the row lower bound
        row_of_ecol.push_back(iRow);
        if (has_row_names)
          ecol_name.push_back("row_" + std::to_string(iRow) + "_" +
                              lp.row_names_[iRow] + "_lower");
        bound_of_row_of_ecol_is_lower.push_back(true);
        // Define the sub-matrix column
        ecol_index.push_back(iRow);
        ecol_value.push_back(1);
        ecol_start.push_back(ecol_index.size());
        ecol_cost.push_back(penalty);
        evar_ix++;
      }
      if (upper < kHighsInf) {
        // Create an e-var for the row upper bound
        row_of_ecol.push_back(iRow);
        if (has_row_names)
          ecol_name.push_back("row_" + std::to_string(iRow) + "_" +
                              lp.row_names_[iRow] + "_upper");
        bound_of_row_of_ecol_is_lower.push_back(false);
        // Define the sub-matrix column
        ecol_index.push_back(iRow);
        ecol_value.push_back(-1);
        ecol_start.push_back(ecol_index.size());
        ecol_cost.push_back(penalty);
        evar_ix++;
      }
    }
    HighsInt num_new_col = ecol_start.size() - 1;
    HighsInt num_new_nz = ecol_start[num_new_col];
    ecol_lower.assign(num_new_col, 0);
    ecol_upper.assign(num_new_col, kHighsInf);
    assert(ecol_cost.size() == static_cast<size_t>(num_new_col));
    assert(ecol_lower.size() == static_cast<size_t>(num_new_col));
    assert(ecol_upper.size() == static_cast<size_t>(num_new_col));
    assert(ecol_start.size() == static_cast<size_t>(num_new_col + 1));
    assert(ecol_index.size() == static_cast<size_t>(num_new_nz));
    assert(ecol_value.size() == static_cast<size_t>(num_new_nz));
    run_status = this->addCols(num_new_col, ecol_cost.data(), ecol_lower.data(),
                               ecol_upper.data(), num_new_nz, ecol_start.data(),
                               ecol_index.data(), ecol_value.data());
    assert(run_status == HighsStatus::kOk);
    if (has_row_names) {
      for (HighsInt iCol = 0; iCol < num_new_col; iCol++)
        this->passColName(row_ecol_offset + iCol, ecol_name[iCol]);
    }

    if (write_model) {
      bool output_flag;
      printf("\nAfter adding %d e-cols\n=============\n", int(num_new_col));
      run_status = this->getOptionValue("output_flag", output_flag);
      this->setOptionValue("output_flag", true);
      this->writeModel("");
      this->setOptionValue("output_flag", output_flag);
    }
  }

  if (write_model) this->writeModel("elastic.mps");
  // Initial logging
  if (get_iis) {
    highsLogUser(
        options_.log_options, HighsLogType::kInfo,
        "Running elasticity filter to identify an infeasible subset of rows\n");
    iis.reportIteration(options_, HighsInt(0), HighsInt(0), true);
  }
  // Working with a copy of the IIS so clear this->iis_
  this->iis_.clear();
  // Lambda for gathering data when solving an LP
  auto solveLp = [&]() -> HighsStatus {
    bool output_flag = options_.output_flag;
    options_.output_flag = false;
    double simplex_time = -this->getRunTime();
    HighsInt simplex_iterations = -info_.simplex_iteration_count;
    run_status = this->optimizeModel();
    simplex_time += this->getRunTime();
    simplex_iterations += info_.simplex_iteration_count;
    iis.info_.update(simplex_time, simplex_iterations);
    options_.output_flag = output_flag;
    return run_status;
  };

  // Solve the elastic LP
  run_status = solveLp();
  if (run_status != HighsStatus::kOk) {
    // Failed to establish feasibility or infeasibility due to unknown error
    // or hitting a time limit
    if (this->model_status_ == HighsModelStatus::kTimeLimit) {
      iis.status_ = kIisModelStatusTimeLimit;
      if (get_iis)
        highsLogUser(
            options_.log_options, HighsLogType::kError,
            "Elasticity filter failed because time limit was reached\n");
    } else {
      iis.status_ = kIisModelStatusUnknown;
      if (get_iis)
        highsLogUser(options_.log_options, HighsLogType::kError,
                     "Elasticity filter failed because it encountered an "
                     "unknown model status\n");
    }
    iis.valid_ = false;
    this->iis_ = iis;
    return elasticityFilterReturn(
        HighsStatus::kError, original_model_name, original_model_status,
        original_num_col, original_num_row, original_sense, original_col_cost,
        original_col_lower, original_col_upper, original_integrality);
  }
  if (kIisDevReport) this->writeSolution("", kSolutionStylePretty);
  // Model status should be optimal, unless model is unbounded
  assert(this->model_status_ == HighsModelStatus::kOptimal ||
         this->model_status_ == HighsModelStatus::kUnbounded);
  iis.valid_ = true;
  iis.status_ = this->info_.objective_function_value > 0
                    ? kIisModelStatusReducible
                    : kIisModelStatusFeasible;
  if (!get_iis) {
    this->iis_ = iis;
    return elasticityFilterReturn(
        HighsStatus::kOk, original_model_name, original_model_status,
        original_num_col, original_num_row, original_sense, original_col_cost,
        original_col_lower, original_col_upper, original_integrality);
  }
  // If getting an (I)IS, assume no elastic columns, so no additional rows
  assert(!has_elastic_columns);
  assert(has_elastic_rows);
  assert(original_num_row == lp.num_row_);

  // Get solution
  const HighsSolution& solution = this->getSolution();
  // Now fix e-variables that are positive and re-solve until e-LP is
  // infeasible
  HighsInt loop_k = 0;
  std::unordered_set<HighsInt> row_set;

  for (;;) {
    loop_k++;
    if (kIisDevReport)
      printf("\nElasticity filter pass %d\n==============\n", int(loop_k));
    HighsInt num_fixed = 0;
    if (has_elastic_columns) {
      for (size_t eCol = 0; eCol < col_of_ecol.size(); eCol++) {
        HighsInt iCol = col_of_ecol[eCol];
        if (solution.col_value[col_ecol_offset + eCol] >
            this->options_.primal_feasibility_tolerance) {
          if (kIisDevReport)
            printf(
                "E-col %2d (column %2d) corresponds to column %2d with %s "
                "bound "
                "%11.4g "
                "and has solution value %11.4g\n",
                int(eCol), int(col_ecol_offset + eCol), int(iCol),
                bound_of_col_of_ecol_is_lower[eCol] ? "lower" : "upper",
                bound_of_col_of_ecol_is_lower[eCol] ? lp.col_lower_[iCol]
                                                    : lp.col_upper_[iCol],
                solution.col_value[col_ecol_offset + eCol]);
          this->changeColBounds(col_ecol_offset + eCol, 0, 0);
          num_fixed++;
        }
      }
    }
    if (has_elastic_rows) {
      for (size_t eCol = 0; eCol < row_of_ecol.size(); eCol++) {
        HighsInt iRow = row_of_ecol[eCol];
        if (solution.col_value[row_ecol_offset + eCol] >
            this->options_.primal_feasibility_tolerance) {
          if (kIisDevReport)
            printf(
                "E-row %2d (column %2d) corresponds to    row %2d with %s "
                "bound "
                "%11.4g "
                "and has solution value %11.4g\n",
                int(eCol), int(row_ecol_offset + eCol), int(iRow),
                bound_of_row_of_ecol_is_lower[eCol] ? "lower" : "upper",
                bound_of_row_of_ecol_is_lower[eCol] ? lp.row_lower_[iRow]
                                                    : lp.row_upper_[iRow],
                solution.col_value[row_ecol_offset + eCol]);
          this->changeColBounds(row_ecol_offset + eCol, 0, 0);
          num_fixed++;
          row_set.insert(iRow);
        }
      }
    }

    if (num_fixed == 0) {
      // No elastic variables were positive, so problem is feasible
      iis.status_ = kIisModelStatusFeasible;
      iis.reportIteration(options_, loop_k, HighsInt(row_set.size()), true);
      break;
    }
    HighsStatus run_status = solveLp();
    HighsModelStatus model_status = this->getModelStatus();
    if (run_status != HighsStatus::kOk) {
      // Solve failed
      if (model_status == HighsModelStatus::kTimeLimit) {
        iis.status_ = kIisModelStatusTimeLimit;

        highsLogUser(
            options_.log_options, HighsLogType::kError,
            "Elasticity filter failed because time limit was reached\n");
      } else {
        iis.status_ = kIisModelStatusUnknown;

        highsLogUser(options_.log_options, HighsLogType::kError,
                     "Elasticity filter failed because it encountered an "
                     "unknown model status\n");
      }
      iis.valid_ = false;
      this->iis_ = iis;
      return elasticityFilterReturn(
          HighsStatus::kError, original_model_name, original_model_status,
          original_num_col, original_num_row, original_sense, original_col_cost,
          original_col_lower, original_col_upper, original_integrality);
    }
    if (kIisDevReport) this->writeSolution("", kSolutionStylePretty);
    bool terminate = (model_status == HighsModelStatus::kInfeasible);
    iis.reportIteration(options_, loop_k, HighsInt(row_set.size()), terminate);
    if (terminate) break;
  }
  HighsInt num_enforced_col_ecol = 0;
  HighsInt num_enforced_row_ecol = 0;
  if (has_elastic_columns) {
    for (size_t eCol = 0; eCol < col_of_ecol.size(); eCol++) {
      HighsInt iCol = col_of_ecol[eCol];
      if (lp.col_upper_[col_ecol_offset + eCol] == 0) {
        num_enforced_col_ecol++;
        if (kIisDevReport)
          printf(
              "Col e-col %2d (column %2d) corresponds to column %2d with %s "
              "bound "
              "%11.4g "
              "and is enforced\n",
              int(eCol), int(col_ecol_offset + eCol), int(iCol),
              bound_of_col_of_ecol_is_lower[eCol] ? "lower" : "upper",
              bound_of_col_of_ecol_is_lower[eCol] ? lp.col_lower_[iCol]
                                                  : lp.col_upper_[iCol]);
      }
    }
  }
  if (has_elastic_rows) {
    for (size_t eCol = 0; eCol < row_of_ecol.size(); eCol++) {
      HighsInt iRow = row_of_ecol[eCol];
      if (lp.col_upper_[row_ecol_offset + eCol] == 0) {
        num_enforced_row_ecol++;
        if (iis.row_index_.empty() || iis.row_index_.back() != iRow)
          iis.row_index_.push_back(iRow);

        /*
        iis.row_bound_.push_back(erow_value[eRow] > 0 ?
                                       kIisBoundStatusUpper :
                                       kIisBoundStatusLower);
        */
        if (kIisDevReport)
          printf(
              "Row e-col %2d (column %2d) corresponds to    row %2d with %s "
              "bound "
              "%11.4g and is enforced\n",
              int(eCol), int(row_ecol_offset + eCol), int(iRow),
              bound_of_row_of_ecol_is_lower[eCol] ? "lower" : "upper",
              bound_of_row_of_ecol_is_lower[eCol] ? lp.row_lower_[iRow]
                                                  : lp.row_upper_[iRow]);
      }
    }
  }
  HighsInt num_iis_row = iis.row_index_.size();
  if (iis.status_ == kIisModelStatusFeasible) {
    assert(num_enforced_col_ecol == 0 && num_enforced_row_ecol == 0);
    assert(num_iis_row == 0);
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Elasticity filter failed to reproduce infeasibility\n");
  } else {
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Elasticity filter after %d passes found an infeasible subset "
                 "of %d rows\n",
                 int(loop_k), row_set.size());
  }

  iis.valid_ = true;
  iis.strategy_ = this->options_.iis_strategy;
  if (iis.status_ == kIisModelStatusFeasible) {
    this->iis_ = iis;
    return elasticityFilterReturn(
        HighsStatus::kOk, original_model_name, original_model_status,
        original_num_col, original_num_row, original_sense, original_col_cost,
        original_col_lower, original_col_upper, original_integrality);
  }

  // Model is infeasible because there are (at least) a positive
  // number of rows in the infeasibility set. Hence the IIS status is
  // reducible
  assert(iis.status_ == kIisModelStatusReducible);
  assert(num_iis_row > 0);

  // Have to be able to map from original row indices to position in
  // iis.row_index_
  std::vector<HighsInt> in_row_index(original_num_row, -1);
  for (HighsInt iX = 0; iX < num_iis_row; iX++)
    in_row_index[iis.row_index_[iX]] = iX;

  // Determine the columns with nonzeros in the row subset
  std::vector<HighsBool> nonzero_in_row_index(original_num_col, false);
  if (lp.a_matrix_.isColwise()) {
    for (HighsInt iCol = 0; iCol < original_num_col; iCol++) {
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
           iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
        if (in_row_index[lp.a_matrix_.index_[iEl]] >= 0)
          nonzero_in_row_index[iCol] = true;
    }
  } else {
    for (HighsInt iX = 0; iX < num_iis_row; iX++) {
      HighsInt iRow = iis.row_index_[iX];
      for (HighsInt iEl = lp.a_matrix_.start_[iRow];
           iEl < lp.a_matrix_.start_[iRow + 1]; iEl++)
        nonzero_in_row_index[lp.a_matrix_.index_[iEl]] = true;
    }
  }
  // Now form iis.col_index_ and iis.col_bound_
  assert(iis.col_index_.size() == 0);
  for (HighsInt iCol = 0; iCol < original_num_col; iCol++) {
    if (nonzero_in_row_index[iCol]) {
      iis.col_index_.push_back(iCol);
      if (lp.col_lower_[iCol] > -kHighsInf) {
        if (lp.col_upper_[iCol] < kHighsInf) {
          // Boxed
          iis.col_bound_.push_back(kIisBoundStatusBoxed);
        } else {
          // LB
          iis.col_bound_.push_back(kIisBoundStatusLower);
        }
      } else {
        if (lp.col_upper_[iCol] < kHighsInf) {
          // UB
          iis.col_bound_.push_back(kIisBoundStatusUpper);
        } else {
          // Free
          iis.col_bound_.push_back(kIisBoundStatusFree);
        }
      }
    }
  }
  // Now, go through the elastic columns for the rows, identifying the bound
  // that is required
  iis.row_bound_.assign(num_iis_row, -1);
  for (size_t eCol = 0; eCol < row_of_ecol.size(); eCol++) {
    HighsInt iRow = row_of_ecol[eCol];
    if (lp.col_upper_[row_ecol_offset + eCol] == 0) {
      HighsInt iX = in_row_index[iRow];
      if (iis.row_bound_[iX] == -1) {
        // Unassigned
        iis.row_bound_[iX] = bound_of_row_of_ecol_is_lower[eCol]
                                 ? kIisBoundStatusLower
                                 : kIisBoundStatusUpper;
      } else {
        // Mark row as boxed when encountered twice, indicating both elastic
        // variables were enforced
        iis.row_bound_[iX] = kIisBoundStatusBoxed;
      }
    }
  }
  assert(iis.row_bound_.size() == iis.row_index_.size());
  for (HighsInt iX = 0; iX < num_iis_row; iX++) assert(iis.row_bound_[iX] >= 0);
  this->iis_ = iis;
  return elasticityFilterReturn(
      HighsStatus::kOk, original_model_name, original_model_status,
      original_num_col, original_num_row, original_sense, original_col_cost,
      original_col_lower, original_col_upper, original_integrality);
}

HighsStatus Highs::extractIis(HighsInt& num_iis_col, HighsInt& num_iis_row,
                              HighsInt* iis_col_index, HighsInt* iis_row_index,
                              HighsInt* iis_col_bound,
                              HighsInt* iis_row_bound) {
  assert(this->iis_.valid_);
  num_iis_col = this->iis_.col_index_.size();
  num_iis_row = this->iis_.row_index_.size();
  if (iis_col_index || iis_col_bound) {
    for (HighsInt iCol = 0; iCol < num_iis_col; iCol++) {
      if (iis_col_index) iis_col_index[iCol] = this->iis_.col_index_[iCol];
      if (iis_col_bound) iis_col_bound[iCol] = this->iis_.col_bound_[iCol];
    }
  }
  if (iis_row_index || iis_row_bound) {
    for (HighsInt iRow = 0; iRow < num_iis_row; iRow++) {
      if (iis_row_index) iis_row_index[iRow] = this->iis_.row_index_[iRow];
      if (iis_row_bound) iis_row_bound[iRow] = this->iis_.row_bound_[iRow];
    }
  }
  return HighsStatus::kOk;
}

bool Highs::aFormatOk(const HighsInt num_nz, const HighsInt format) {
  if (!num_nz) return true;
  const bool ok_format = format == (HighsInt)MatrixFormat::kColwise ||
                         format == (HighsInt)MatrixFormat::kRowwise;
  if (!ok_format)
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Non-empty Constraint matrix has illegal format = %d\n",
                 int(format));
  assert(ok_format);
  return ok_format;
}

bool Highs::qFormatOk(const HighsInt num_nz, const HighsInt format) {
  if (!num_nz) return true;
  const bool ok_format =
      format == static_cast<HighsInt>(HessianFormat::kTriangular) ||
      format == static_cast<HighsInt>(HessianFormat::kSquare);
  if (!ok_format)
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Non-empty Hessian matrix has illegal format = %d\n",
                 int(format));
  assert(ok_format);
  return ok_format;
}

void Highs::clearZeroHessian() {
  HighsHessian& hessian = model_.hessian_;
  if (hessian.dim_) {
    // Clear any zero Hessian
    if (hessian.numNz() == 0) {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Hessian has dimension %" HIGHSINT_FORMAT
                   " but no nonzeros, so is ignored\n",
                   hessian.dim_);
      hessian.clear();
    }
  }
}

void Highs::callLpKktCheck(const HighsLp& lp, const std::string& message) {
  lpKktCheck(this->model_status_, this->info_, lp, this->solution_,
             this->basis_, this->options_, message);
}

HighsStatus Highs::invertRequirementError(
    const std::string& method_name) const {
  assert(!ekk_instance_.status_.has_invert);
  highsLogUser(options_.log_options, HighsLogType::kError,
               "No invertible representation for %s\n", method_name.c_str());
  return HighsStatus::kError;
}

HighsStatus Highs::handleInfCost() {
  HighsLp& lp = this->model_.lp_;
  if (!lp.has_infinite_cost_) return HighsStatus::kOk;
  HighsLpMods& mods = lp.mods_;
  double inf_cost = this->options_.infinite_cost;
  for (HighsInt k = 0; k < 2; k++) {
    // Pass twice: first checking that infinite costs can be handled,
    // then handling them, so that model is unmodified if infinite
    // costs cannot be handled
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      double cost = lp.col_cost_[iCol];
      if (cost > -inf_cost && cost < inf_cost) continue;
      double lower = lp.col_lower_[iCol];
      double upper = lp.col_upper_[iCol];
      if (lp.isMip()) {
        if (lp.integrality_[iCol] == HighsVarType::kInteger) {
          lower = std::ceil(lower);
          upper = std::floor(upper);
        }
      }
      if (cost <= -inf_cost) {
        if (lp.sense_ == ObjSense::kMinimize) {
          // Minimizing with -inf cost so try to fix at upper bound
          if (upper < kHighsInf) {
            if (k) lp.col_lower_[iCol] = upper;
          } else {
            highsLogUser(options_.log_options, HighsLogType::kError,
                         "Cannot minimize with a cost on variable %d of %g and "
                         "upper bound of %g\n",
                         int(iCol), cost, upper);
            return HighsStatus::kError;
          }
        } else {
          // Maximizing with -inf cost so try to fix at lower bound
          if (lower > -kHighsInf) {
            if (k) lp.col_upper_[iCol] = lower;
          } else {
            highsLogUser(options_.log_options, HighsLogType::kError,
                         "Cannot maximize with a cost on variable %d of %g and "
                         "lower bound of %g\n",
                         int(iCol), cost, lower);
            return HighsStatus::kError;
          }
        }
      } else {
        if (lp.sense_ == ObjSense::kMinimize) {
          // Minimizing with inf cost so try to fix at lower bound
          if (lower > -kHighsInf) {
            if (k) lp.col_upper_[iCol] = lower;
          } else {
            highsLogUser(options_.log_options, HighsLogType::kError,
                         "Cannot minimize with a cost on variable %d of %g and "
                         "lower bound of %g\n",
                         int(iCol), cost, lower);
            return HighsStatus::kError;
          }
        } else {
          // Maximizing with inf cost so try to fix at upper bound
          if (upper < kHighsInf) {
            if (k) lp.col_lower_[iCol] = upper;
          } else {
            highsLogUser(options_.log_options, HighsLogType::kError,
                         "Cannot maximize with a cost on variable %d of %g and "
                         "upper bound of %g\n",
                         int(iCol), cost, upper);
            return HighsStatus::kError;
          }
        }
      }
      if (k) {
        mods.save_inf_cost_variable_index.push_back(iCol);
        mods.save_inf_cost_variable_cost.push_back(cost);
        mods.save_inf_cost_variable_lower.push_back(lower);
        mods.save_inf_cost_variable_upper.push_back(upper);
        lp.col_cost_[iCol] = 0;
      }
    }
  }
  // Infinite costs have been removed, but their presence in the
  // original model is known from mods.save_inf_cost_variable_*, so
  // set lp.has_infinite_cost_ to be false to avoid assert when run()
  // is called using copy of model in MIP solver (See #1446)
  lp.has_infinite_cost_ = false;

  return HighsStatus::kOk;
}

HighsStatus Highs::optionChangeAction() {
  if (this->iis_.valid_ && options_.iis_strategy != this->iis_.strategy_)
    this->iis_.clear();
  return HighsStatus::kOk;
}

void Highs::restoreInfCost(HighsStatus& return_status) {
  HighsLp& lp = this->model_.lp_;
  HighsBasis& basis = this->basis_;
  HighsLpMods& mods = lp.mods_;
  HighsInt num_inf_cost = mods.save_inf_cost_variable_index.size();
  if (num_inf_cost <= 0) return;
  assert(num_inf_cost);
  for (HighsInt ix = 0; ix < num_inf_cost; ix++) {
    HighsInt iCol = mods.save_inf_cost_variable_index[ix];
    double cost = mods.save_inf_cost_variable_cost[ix];
    double lower = mods.save_inf_cost_variable_lower[ix];
    double upper = mods.save_inf_cost_variable_upper[ix];
    double value = solution_.value_valid ? solution_.col_value[iCol] : 0;
    if (basis.valid) {
      assert(basis.col_status[iCol] != HighsBasisStatus::kBasic);
      if (lp.col_lower_[iCol] == lower) {
        basis.col_status[iCol] = HighsBasisStatus::kLower;
      } else {
        basis.col_status[iCol] = HighsBasisStatus::kUpper;
      }
    }
    assert(lp.col_cost_[iCol] == 0);
    if (value) this->info_.objective_function_value += value * cost;
    lp.col_cost_[iCol] = cost;
    lp.col_lower_[iCol] = lower;
    lp.col_upper_[iCol] = upper;
  }
  // Infinite costs have been reintroduced, so reset to true the flag
  // that was set false in Highs::handleInfCost() (See #1446)
  lp.has_infinite_cost_ = true;

  if (this->model_status_ == HighsModelStatus::kInfeasible) {
    // Model is infeasible with the infinite cost variables fixed at
    // appropriate values, so model status cannot be determined
    this->model_status_ = HighsModelStatus::kUnknown;
    setHighsModelStatusAndClearSolutionAndBasis(this->model_status_);
    return_status = highsStatusFromHighsModelStatus(model_status_);
  }
}

HighsStatus Highs::userScale(HighsUserScaleData& data) {
  if (!options_.user_objective_scale && !options_.user_bound_scale)
    return HighsStatus::kOk;
  if (this->model_.hessian_.isOracle()) {
    highsLogUser(this->options_.log_options, HighsLogType::kWarning,
                 "Hessian is represented via an oracle, so user scaling cannot "
                 "be applied\n");
    return HighsStatus::kWarning;
  }
  // User objective and bound scaling data are accumulated in the
  // HighsUserScaleData struct, in particular, there is a local copy
  // of the user objective and bound scaling options values, and
  // records of resulting extreme data values that prevent the user
  // objective and bound scaling from being applied.
  initialiseUserScaleData(this->options_, data);
  // Determine whether user scaling yields excessively large cost,
  // Hessian values, column/row bounds or matrix values. If not,
  // then apply the user scaling to the model...
  if (this->userScaleModel(data) == HighsStatus::kError)
    return HighsStatus::kError;
  // ... and the solution
  this->userScaleSolution(data);
  // Indicate that the scaling has been applied
  data.applied = true;
  return HighsStatus::kOk;
}

HighsStatus Highs::userUnscale(HighsUserScaleData& data) {
  if (!data.applied) return HighsStatus::kOk;
  // Unscale the incumbent model and solution
  HighsStatus status = HighsStatus::kOk;
  // Flip the scaling sign
  data.user_objective_scale *= -1;
  data.user_bound_scale *= -1;
  HighsStatus unscale_status = this->userScaleModel(data);
  if (unscale_status == HighsStatus::kError) {
    highsLogUser(
        this->options_.log_options, HighsLogType::kError,
        "Unexpected error removing user scaling from the incumbent model\n");
    assert(unscale_status != HighsStatus::kError);
  }
  const bool update_kkt = true;
  unscale_status = this->userScaleSolution(data, update_kkt);
  highsLogUser(this->options_.log_options, HighsLogType::kInfo,
               "After solving the user-scaled model, the unscaled solution "
               "has objective value %.12g\n",
               this->info_.objective_function_value);
  if (model_status_ == HighsModelStatus::kOptimal &&
      unscale_status != HighsStatus::kOk) {
    // KKT errors in the unscaled optimal solution, so log a warning and
    // return
    highsLogUser(
        this->options_.log_options, HighsLogType::kWarning,
        "User scaled problem solved to optimality, but unscaled solution "
        "does not satisfy feasibility and optimality tolerances\n");
    status = HighsStatus::kWarning;
  }
  return status;
}

HighsStatus Highs::userScaleModel(HighsUserScaleData& data) {
  // Consider applying user objective and bound scaling to the model
  // by first identifying whether it causes any errors due to creating
  // extreme data values...
  userScaleLp(this->model_.lp_, data, false);
  userScaleHessian(this->model_.hessian_, data, false);
  HighsStatus return_status = userScaleStatus(this->options_.log_options, data);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  // ... and, if not, actually apply the scaling
  userScaleLp(this->model_.lp_, data);
  userScaleHessian(this->model_.hessian_, data);
  return return_status;
}

HighsStatus Highs::userScaleSolution(HighsUserScaleData& data,
                                     bool update_kkt) {
  HighsStatus return_status = HighsStatus::kOk;
  if (!data.user_objective_scale && !data.user_bound_scale)
    return HighsStatus::kOk;
  double objective_scale_value = std::pow(2, data.user_objective_scale);
  double bound_scale_value = std::pow(2, data.user_bound_scale);
  const HighsLp& lp = this->model_.lp_;
  const bool has_integrality = lp.integrality_.size();
  if (info_.primal_solution_status != kSolutionStatusNone) {
    if (data.user_bound_scale) {
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        if (has_integrality &&
            lp.integrality_[iCol] != HighsVarType::kContinuous)
          continue;
        this->solution_.col_value[iCol] *= bound_scale_value;
      }
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
        this->solution_.row_value[iRow] *= bound_scale_value;
    }
  }
  if (info_.dual_solution_status != kSolutionStatusNone) {
    if (data.user_objective_scale) {
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
        this->solution_.col_dual[iCol] *= objective_scale_value;
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
        this->solution_.row_dual[iRow] *= objective_scale_value;
    }
  }
  if (!update_kkt) return return_status;
  // In scaling the objective function value, have to consider the offset
  info_.objective_function_value *= (bound_scale_value * objective_scale_value);
  if (has_integrality)
    info_.mip_dual_bound *= (bound_scale_value * objective_scale_value);

  const double primal_feasibility_tolerance =
      options_.primal_feasibility_tolerance;
  // Remember to check primal feasibility using the tolerance used
  // when the model was solved
  if (model_.isMip() && !options_.solve_relaxation)
    options_.primal_feasibility_tolerance =
        options_.getMipPrimalFeasibilityTolerance();
  getKktFailures(options_, model_, solution_, basis_, info_);
  // Restore the primal feasibility tolerance
  options_.primal_feasibility_tolerance = primal_feasibility_tolerance;
  return reportKktFailures(model_.lp_, options_, info_,
                           "After removing user scaling")
             ? HighsStatus::kWarning
             : return_status;
}

void HighsIllConditioning::clear() { this->record.clear(); }

HighsStatus Highs::computeIllConditioning(
    HighsIllConditioning& ill_conditioning, const bool constraint,
    const HighsInt method, const double ill_conditioning_bound) {
  const double kZeroMultiplier = 1e-6;
  ill_conditioning.clear();
  HighsLp& incumbent_lp = this->model_.lp_;
  Highs conditioning;
  const bool dev_conditioning = false;
  conditioning.setOptionValue("output_flag", false);  // dev_conditioning);
  std::vector<HighsInt> basic_var;
  HighsLp& ill_conditioning_lp = conditioning.model_.lp_;
  // Form the ill-conditioning LP according to method
  if (method == 0) {
    formIllConditioningLp0(ill_conditioning_lp, basic_var, constraint);
  } else {
    formIllConditioningLp1(ill_conditioning_lp, basic_var, constraint,
                           ill_conditioning_bound);
    //    conditioning.writeModel("");
  }

  //  if (dev_conditioning) conditioning.writeModel("");

  assert(assessLp(ill_conditioning_lp, this->options_) == HighsStatus::kOk);
  // Solve the ill-conditioning analysis LP
  HighsStatus return_status = conditioning.run();
  HighsModelStatus model_status = conditioning.getModelStatus();
  const std::string type = constraint ? "Constraint" : "Column";
  const bool failed =
      return_status != HighsStatus::kOk ||
      (method == 0 && model_status != HighsModelStatus::kOptimal) ||
      (method == 1 && (model_status != HighsModelStatus::kOptimal &&
                       model_status != HighsModelStatus::kInfeasible));
  if (failed) {
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "\n%s view ill-conditioning analysis has failed\n",
                 type.c_str());
    return HighsStatus::kError;
  }
  if (method == 1 && model_status == HighsModelStatus::kInfeasible) {
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "\n%s view ill-conditioning bound of %g is insufficient for "
                 "analysis: try %g\n",
                 type.c_str(), ill_conditioning_bound,
                 1e1 * ill_conditioning_bound);
    return HighsStatus::kOk;
  }
  if (dev_conditioning) conditioning.writeSolution("", 1);
  // Extract and normalise the multipliers
  HighsSolution& solution = conditioning.solution_;
  double multiplier_norm = 0;
  for (HighsInt iRow = 0; iRow < incumbent_lp.num_row_; iRow++)
    multiplier_norm += std::fabs(solution.col_value[iRow]);
  assert(multiplier_norm > 0);
  const double ill_conditioning_measure =
      (method == 0 ? conditioning.getInfo().objective_function_value
                   : solution.row_value[conditioning.getNumRow() - 1]) /
      multiplier_norm;
  highsLogUser(
      options_.log_options, HighsLogType::kInfo,
      "\n%s view ill-conditioning analysis: 1-norm distance of basis matrix "
      "from singularity is estimated to be %g\n",
      type.c_str(), ill_conditioning_measure);
  std::vector<std::pair<double, HighsInt>> abs_list;
  for (HighsInt iRow = 0; iRow < incumbent_lp.num_row_; iRow++) {
    double abs_multiplier =
        std::fabs(solution.col_value[iRow]) / multiplier_norm;
    if (abs_multiplier <= kZeroMultiplier) continue;
    abs_list.push_back(std::make_pair(abs_multiplier, iRow));
  }
  std::sort(abs_list.begin(), abs_list.end());
  // Report on ill-conditioning multipliers
  std::stringstream ss;
  const bool has_row_names =
      HighsInt(incumbent_lp.row_names_.size()) == incumbent_lp.num_row_;
  const bool has_col_names =
      HighsInt(incumbent_lp.col_names_.size()) == incumbent_lp.num_col_;
  const double coefficient_zero_tolerance = 1e-8;
  auto printCoefficient = [&](const double multiplier, const bool first) {
    if (std::fabs(multiplier) < coefficient_zero_tolerance) {
      ss << "+ 0";
    } else if (std::fabs(multiplier - 1) < coefficient_zero_tolerance) {
      std::string str = first ? "" : "+ ";
      ss << str;
    } else if (std::fabs(multiplier + 1) < coefficient_zero_tolerance) {
      std::string str = first ? "-" : "- ";
      ss << str;
    } else if (multiplier < 0) {
      std::string str = first ? "-" : "- ";
      ss << str << -multiplier << " ";
    } else {
      std::string str = first ? "" : "+ ";
      ss << str << multiplier << " ";
    }
  };

  for (HighsInt iX = int(abs_list.size()) - 1; iX >= 0; iX--) {
    HighsInt iRow = abs_list[iX].second;
    HighsIllConditioningRecord record;
    record.index = iRow;
    record.multiplier = solution.col_value[iRow] / multiplier_norm;
    ill_conditioning.record.push_back(record);
  }
  HighsSparseMatrix& incumbent_matrix = incumbent_lp.a_matrix_;
  if (constraint) {
    HighsInt num_nz;
    std::vector<HighsInt> index(incumbent_lp.num_col_);
    std::vector<double> value(incumbent_lp.num_col_);
    HighsInt* p_index = index.data();
    double* p_value = value.data();
    for (HighsInt iX = 0; iX < HighsInt(ill_conditioning.record.size()); iX++) {
      ss.str(std::string());
      bool newline = false;
      HighsInt iRow = ill_conditioning.record[iX].index;
      double multiplier = ill_conditioning.record[iX].multiplier;
      // Extract the row corresponding to this constraint
      num_nz = 0;
      incumbent_matrix.getRow(iRow, num_nz, p_index, p_value);
      std::string row_name = has_row_names ? incumbent_lp.row_names_[iRow]
                                           : "R" + std::to_string(iRow);
      ss << "(Mu=" << multiplier << ")" << row_name << ": ";
      const double lower = incumbent_lp.row_lower_[iRow];
      const double upper = incumbent_lp.row_upper_[iRow];
      if (lower > -kHighsInf && lower != upper)
        ss << incumbent_lp.row_lower_[iRow] << " <= ";
      for (HighsInt iEl = 0; iEl < num_nz; iEl++) {
        if (newline) {
          ss << "  ";
          newline = false;
        }
        HighsInt iCol = index[iEl];
        printCoefficient(value[iEl], iEl == 0);
        std::string col_name = has_col_names ? incumbent_lp.col_names_[iCol]
                                             : "C" + std::to_string(iCol);
        ss << col_name << " ";
        HighsInt length_ss = ss.str().length();
        if (length_ss > 72 && iEl < num_nz - 1) {
          highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                       ss.str().c_str());
          ss.str(std::string());
          newline = true;
        }
      }
      if (upper < kHighsInf) {
        if (lower == upper) {
          ss << "= " << upper;
        } else {
          ss << "<= " << upper;
        }
      }
      if (ss.str().length())
        highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                     ss.str().c_str());
    }
  } else {
    for (const auto& rec : ill_conditioning.record) {
      ss.str(std::string());
      bool newline = false;
      double multiplier = rec.multiplier;
      HighsInt iCol = basic_var[rec.index];
      if (iCol < incumbent_lp.num_col_) {
        std::string col_name = has_col_names ? incumbent_lp.col_names_[iCol]
                                             : "C" + std::to_string(iCol);
        ss << "(Mu=" << multiplier << ")" << col_name << ": ";
        for (HighsInt iEl = incumbent_matrix.start_[iCol];
             iEl < incumbent_matrix.start_[iCol + 1]; iEl++) {
          if (newline) {
            ss << "  ";
            newline = false;
          } else {
            if (iEl > incumbent_matrix.start_[iCol]) ss << " | ";
          }
          HighsInt iRow = incumbent_matrix.index_[iEl];
          printCoefficient(incumbent_matrix.value_[iEl], true);
          std::string row_name = has_row_names ? incumbent_lp.row_names_[iRow]
                                               : "R" + std::to_string(iRow);
          ss << row_name;
          HighsInt length_ss = ss.str().length();
          if (length_ss > 72 && iEl < incumbent_matrix.start_[iCol + 1] - 1) {
            ss << " | ";
            highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                         ss.str().c_str());
            ss.str(std::string());
            newline = true;
          }
        }
      } else {
        HighsInt iRow = iCol - incumbent_lp.num_col_;
        std::string col_name = has_row_names
                                   ? "Slack_" + incumbent_lp.row_names_[iRow]
                                   : "Slack_R" + std::to_string(iRow);
        ss << "(Mu=" << multiplier << ")" << col_name << ": ";
      }
      if (ss.str().length())
        highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                     ss.str().c_str());
    }
  }
  return HighsStatus::kOk;
}

void Highs::formIllConditioningLp0(HighsLp& ill_conditioning_lp,
                                   std::vector<HighsInt>& basic_var,
                                   const bool constraint) {
  HighsLp& incumbent_lp = this->model_.lp_;
  // Conditioning LP minimizes the infeasibilities of
  //
  // [B^T]y = [0]; y free - for constraint view
  // [e^T]    [1]
  //
  // [ B ]y = [0]; y free - for column view
  // [e^T]    [1]
  //
  ill_conditioning_lp.num_row_ = incumbent_lp.num_row_ + 1;
  for (HighsInt iRow = 0; iRow < incumbent_lp.num_row_; iRow++) {
    ill_conditioning_lp.row_lower_.push_back(0);
    ill_conditioning_lp.row_upper_.push_back(0);
  }
  ill_conditioning_lp.row_lower_.push_back(1);
  ill_conditioning_lp.row_upper_.push_back(1);
  HighsSparseMatrix& incumbent_matrix = incumbent_lp.a_matrix_;
  incumbent_matrix.ensureColwise();
  HighsSparseMatrix& ill_conditioning_matrix = ill_conditioning_lp.a_matrix_;
  ill_conditioning_matrix.num_row_ = ill_conditioning_lp.num_row_;
  // Form the basis matrix and
  //
  // * For constraint view, add the column e, and transpose the
  // * resulting matrix
  //
  // * For column view, add a unit entry to each column
  //
  const HighsInt ill_conditioning_lp_e_row = ill_conditioning_lp.num_row_ - 1;
  for (HighsInt iCol = 0; iCol < incumbent_lp.num_col_; iCol++) {
    if (this->basis_.col_status[iCol] != HighsBasisStatus::kBasic) continue;
    // Basic column goes into conditioning LP, possibly with unit
    // coefficient for constraint e^Ty=1
    basic_var.push_back(iCol);
    ill_conditioning_lp.col_cost_.push_back(0);
    ill_conditioning_lp.col_lower_.push_back(-kHighsInf);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    for (HighsInt iEl = incumbent_matrix.start_[iCol];
         iEl < incumbent_matrix.start_[iCol + 1]; iEl++) {
      ill_conditioning_matrix.index_.push_back(incumbent_matrix.index_[iEl]);
      ill_conditioning_matrix.value_.push_back(incumbent_matrix.value_[iEl]);
    }
    if (!constraint) {
      ill_conditioning_matrix.index_.push_back(ill_conditioning_lp_e_row);
      ill_conditioning_matrix.value_.push_back(1.0);
    }
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
  }
  for (HighsInt iRow = 0; iRow < incumbent_lp.num_row_; iRow++) {
    if (this->basis_.row_status[iRow] != HighsBasisStatus::kBasic) continue;
    // Basic slack goes into conditioning LP
    basic_var.push_back(incumbent_lp.num_col_ + iRow);
    ill_conditioning_lp.col_cost_.push_back(0);
    ill_conditioning_lp.col_lower_.push_back(-kHighsInf);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    ill_conditioning_matrix.index_.push_back(iRow);
    ill_conditioning_matrix.value_.push_back(-1.0);
    if (!constraint) {
      ill_conditioning_matrix.index_.push_back(ill_conditioning_lp_e_row);
      ill_conditioning_matrix.value_.push_back(1.0);
    }
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
  }
  if (constraint) {
    // Add the column e, and transpose the resulting matrix
    for (HighsInt iRow = 0; iRow < incumbent_lp.num_row_; iRow++) {
      ill_conditioning_matrix.index_.push_back(iRow);
      ill_conditioning_matrix.value_.push_back(1.0);
    }
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    ill_conditioning_matrix.num_row_ = incumbent_lp.num_row_;
    ill_conditioning_matrix.num_col_ = incumbent_lp.num_row_ + 1;
    ill_conditioning_matrix.ensureRowwise();
    ill_conditioning_matrix.format_ = MatrixFormat::kColwise;
  }
  // Now add the variables to measure the infeasibilities
  for (HighsInt iRow = 0; iRow < incumbent_lp.num_row_; iRow++) {
    // Adding x_+ with cost 1
    ill_conditioning_lp.col_cost_.push_back(1);
    ill_conditioning_lp.col_lower_.push_back(0);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    ill_conditioning_matrix.index_.push_back(iRow);
    ill_conditioning_matrix.value_.push_back(1.0);
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    // Subtracting x_- with cost 1
    ill_conditioning_lp.col_cost_.push_back(1);
    ill_conditioning_lp.col_lower_.push_back(0);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    ill_conditioning_matrix.index_.push_back(iRow);
    ill_conditioning_matrix.value_.push_back(-1.0);
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
  }
  ill_conditioning_lp.num_col_ = 3 * incumbent_lp.num_row_;
  ill_conditioning_matrix.num_col_ = ill_conditioning_lp.num_col_;
  ill_conditioning_matrix.num_row_ = ill_conditioning_lp.num_row_;
}

void Highs::formIllConditioningLp1(HighsLp& ill_conditioning_lp,
                                   std::vector<HighsInt>& basic_var,
                                   const bool constraint,
                                   const double ill_conditioning_bound) {
  HighsLp& incumbent_lp = this->model_.lp_;
  const HighsInt incumbent_num_row = incumbent_lp.num_row_;
  //
  // Using notation from Klotz14
  //
  // For constraint view, conditioning LP minimizes the
  // infeasibilities of c7
  //
  // c4: B^Ty         -   s +   t   = 0
  // c1:    y - u + w               = 0
  // c7:        u + w               = 0
  // c6: e^Ty                       = 1
  // c5:               e^Ts + e^Tt <= eps
  // y free; u, w, s, t >= 0
  //
  // Column view uses B rather than B^T
  //
  // Set up offsets
  //
  const HighsInt c4_offset = 0;
  const HighsInt c1_offset = incumbent_num_row;
  const HighsInt c7_offset = 2 * incumbent_num_row;
  const HighsInt c6_offset = 3 * incumbent_num_row;
  const HighsInt c5_offset = 3 * incumbent_num_row + 1;
  for (HighsInt iRow = 0; iRow < c6_offset; iRow++) {
    ill_conditioning_lp.row_lower_.push_back(0);
    ill_conditioning_lp.row_upper_.push_back(0);
  }
  HighsSparseMatrix& incumbent_matrix = incumbent_lp.a_matrix_;
  incumbent_matrix.ensureColwise();
  HighsSparseMatrix& ill_conditioning_matrix = ill_conditioning_lp.a_matrix_;
  // Form the basis matrix and
  //
  // * For constraint view, add the identity matrix and vector of
  // * ones, and transpose the resulting matrix
  //
  // * For column view, add an identity matrix column and unit entry
  // * below each column
  //
  ill_conditioning_lp.num_col_ = 0;
  for (HighsInt iCol = 0; iCol < incumbent_lp.num_col_; iCol++) {
    if (this->basis_.col_status[iCol] != HighsBasisStatus::kBasic) continue;
    // Basic column goes into ill-conditioning LP, possibly with
    // identity matrix column for constraint y - u + w = 0 and unit
    // entry for e^Ty = 1
    basic_var.push_back(iCol);
    ill_conditioning_lp.col_names_.push_back(
        "y_" + std::to_string(ill_conditioning_lp.num_col_));
    ill_conditioning_lp.col_cost_.push_back(0);
    ill_conditioning_lp.col_lower_.push_back(-kHighsInf);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    for (HighsInt iEl = incumbent_matrix.start_[iCol];
         iEl < incumbent_matrix.start_[iCol + 1]; iEl++) {
      ill_conditioning_matrix.index_.push_back(incumbent_matrix.index_[iEl]);
      ill_conditioning_matrix.value_.push_back(incumbent_matrix.value_[iEl]);
    }
    if (!constraint) {
      // Add identity matrix column for constraint y - u + w = 0
      ill_conditioning_matrix.index_.push_back(c1_offset +
                                               ill_conditioning_lp.num_col_);
      ill_conditioning_matrix.value_.push_back(1.0);
      // Add unit entry for e^Ty = 1
      ill_conditioning_matrix.index_.push_back(c6_offset);
      ill_conditioning_matrix.value_.push_back(1.0);
    }
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    ill_conditioning_lp.num_col_++;
  }

  for (HighsInt iRow = 0; iRow < incumbent_num_row; iRow++) {
    if (this->basis_.row_status[iRow] != HighsBasisStatus::kBasic) continue;
    // Basic slack goes into conditioning LP
    basic_var.push_back(incumbent_lp.num_col_ + iRow);
    ill_conditioning_lp.col_names_.push_back(
        "y_" + std::to_string(ill_conditioning_lp.num_col_));
    ill_conditioning_lp.col_cost_.push_back(0);
    ill_conditioning_lp.col_lower_.push_back(-kHighsInf);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    ill_conditioning_matrix.index_.push_back(iRow);
    ill_conditioning_matrix.value_.push_back(-1.0);
    if (!constraint) {
      // Add identity matrix column for constraint y - u + w = 0
      ill_conditioning_matrix.index_.push_back(c1_offset +
                                               ill_conditioning_lp.num_col_);
      ill_conditioning_matrix.value_.push_back(1.0);
      // Add unit entry for e^Ty = 1
      ill_conditioning_matrix.index_.push_back(c6_offset);
      ill_conditioning_matrix.value_.push_back(1.0);
    }
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    ill_conditioning_lp.num_col_++;
  }
  assert(ill_conditioning_lp.num_col_ == incumbent_num_row);
  if (constraint) {
    // Add the identity matrix for constraint y - u + w = 0
    for (HighsInt iRow = 0; iRow < incumbent_num_row; iRow++) {
      ill_conditioning_matrix.index_.push_back(iRow);
      ill_conditioning_matrix.value_.push_back(1.0);
      ill_conditioning_matrix.start_.push_back(
          HighsInt(ill_conditioning_matrix.index_.size()));
    }
    // Add the square zero matrix of c7
    for (HighsInt iRow = 0; iRow < incumbent_num_row; iRow++)
      ill_conditioning_matrix.start_.push_back(
          HighsInt(ill_conditioning_matrix.index_.size()));
    // Add the vector of ones for e^Ty = 1
    for (HighsInt iRow = 0; iRow < incumbent_num_row; iRow++) {
      ill_conditioning_matrix.index_.push_back(iRow);
      ill_conditioning_matrix.value_.push_back(1.0);
    }
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));

    // Transpose the resulting matrix
    ill_conditioning_matrix.num_col_ = c6_offset + 1;
    ill_conditioning_matrix.num_row_ = incumbent_num_row;
    ill_conditioning_matrix.ensureRowwise();
    ill_conditioning_matrix.format_ = MatrixFormat::kColwise;
    ill_conditioning_matrix.num_col_ = incumbent_num_row;
    ill_conditioning_matrix.num_row_ = c6_offset + 1;
  }

  assert(ill_conditioning_lp.num_col_ == incumbent_num_row);
  ill_conditioning_lp.num_row_ = 3 * incumbent_num_row + 2;

  // Now add the variables u and w
  for (HighsInt iRow = 0; iRow < incumbent_num_row; iRow++) {
    // Adding u with cost 0
    ill_conditioning_lp.col_names_.push_back("u_" + std::to_string(iRow));
    ill_conditioning_lp.col_cost_.push_back(0);
    ill_conditioning_lp.col_lower_.push_back(0);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    // Contribution to c1: y - u + w = 0
    ill_conditioning_matrix.index_.push_back(c1_offset + iRow);
    ill_conditioning_matrix.value_.push_back(-1.0);
    // Contribution to c7: u + w = 0
    ill_conditioning_matrix.index_.push_back(c7_offset + iRow);
    ill_conditioning_matrix.value_.push_back(1.0);
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    ill_conditioning_lp.num_col_++;
    // Adding w with cost 0
    ill_conditioning_lp.col_names_.push_back("w_" + std::to_string(iRow));
    ill_conditioning_lp.col_cost_.push_back(0);
    ill_conditioning_lp.col_lower_.push_back(0);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    // Contribution to c1: y - u + w = 0
    ill_conditioning_matrix.index_.push_back(c1_offset + iRow);
    ill_conditioning_matrix.value_.push_back(1.0);
    // Contribution to c7: u + w = 0
    ill_conditioning_matrix.index_.push_back(c7_offset + iRow);
    ill_conditioning_matrix.value_.push_back(1.0);
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    ill_conditioning_lp.num_col_++;
  }
  // Now add the variables s and t
  for (HighsInt iRow = 0; iRow < incumbent_num_row; iRow++) {
    // Adding s with cost 0
    ill_conditioning_lp.col_names_.push_back("s_" + std::to_string(iRow));
    ill_conditioning_lp.col_cost_.push_back(0);
    ill_conditioning_lp.col_lower_.push_back(0);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    // Contribution to c4: B^Ty - s + t = 0
    ill_conditioning_matrix.index_.push_back(c4_offset + iRow);
    ill_conditioning_matrix.value_.push_back(-1.0);
    // Contribution to c5: e^Ts + e^Tt <= eps
    ill_conditioning_matrix.index_.push_back(c5_offset);
    ill_conditioning_matrix.value_.push_back(1.0);
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    ill_conditioning_lp.num_col_++;
    // Adding t with cost 0
    ill_conditioning_lp.col_names_.push_back("t_" + std::to_string(iRow));
    ill_conditioning_lp.col_cost_.push_back(0);
    ill_conditioning_lp.col_lower_.push_back(0);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    // Contribution to c4: B^Ty - s + t = 0
    ill_conditioning_matrix.index_.push_back(c4_offset + iRow);
    ill_conditioning_matrix.value_.push_back(1.0);
    // Contribution to c5: e^Ts + e^Tt <= eps
    ill_conditioning_matrix.index_.push_back(c5_offset);
    ill_conditioning_matrix.value_.push_back(1.0);
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    ill_conditioning_lp.num_col_++;
  }
  // Add the bounds for c6: e^Ty = 1
  ill_conditioning_lp.row_lower_.push_back(1);
  ill_conditioning_lp.row_upper_.push_back(1);
  // Add the bounds for c5: e^Ts + e^Tt <= eps
  assert(ill_conditioning_bound > 0);
  ill_conditioning_lp.row_lower_.push_back(-kHighsInf);
  ill_conditioning_lp.row_upper_.push_back(ill_conditioning_bound);
  assert(HighsInt(ill_conditioning_lp.row_lower_.size()) ==
         ill_conditioning_lp.num_row_);
  assert(HighsInt(ill_conditioning_lp.row_upper_.size()) ==
         ill_conditioning_lp.num_row_);

  // Now add the variables to measure the infeasibilities in
  //
  // c7: u + w = r^+ - r^-
  for (HighsInt iRow = 0; iRow < incumbent_num_row; iRow++) {
    // Adding r^+ with cost 1
    ill_conditioning_lp.col_names_.push_back("IfsPlus_" + std::to_string(iRow));
    ill_conditioning_lp.col_cost_.push_back(1);
    ill_conditioning_lp.col_lower_.push_back(0);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    ill_conditioning_matrix.index_.push_back(c7_offset + iRow);
    ill_conditioning_matrix.value_.push_back(-1.0);
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    ill_conditioning_lp.num_col_++;
    // Adding r^- with cost 1
    ill_conditioning_lp.col_names_.push_back("IfsMinus_" +
                                             std::to_string(iRow));
    ill_conditioning_lp.col_cost_.push_back(1);
    ill_conditioning_lp.col_lower_.push_back(0);
    ill_conditioning_lp.col_upper_.push_back(kHighsInf);
    ill_conditioning_matrix.index_.push_back(c7_offset + iRow);
    ill_conditioning_matrix.value_.push_back(1.0);
    ill_conditioning_matrix.start_.push_back(
        HighsInt(ill_conditioning_matrix.index_.size()));
    ill_conditioning_lp.num_col_++;
  }
  assert(ill_conditioning_lp.num_col_ == 7 * incumbent_num_row);
  assert(ill_conditioning_lp.num_row_ == 3 * incumbent_num_row + 2);
  ill_conditioning_matrix.num_col_ = ill_conditioning_lp.num_col_;
  ill_conditioning_matrix.num_row_ = ill_conditioning_lp.num_row_;
}

// ToDo: this method shouldn't change the incumbent model irrevocably,
// since this loses user data
bool Highs::processedBoundsFindsInfeasibility() {
  const HighsLogOptions& log_options = this->options_.log_options;
  HighsLp& lp = this->model_.lp_;

  HighsInt num_true_infeasible_bound = 0;
  HighsInt num_ok_infeasible_bound = 0;
  const bool has_integrality = lp.integrality_.size() > 0;
  bool performed_inward_integer_rounding = false;
  // Lambda for assessing infeasible bounds
  auto infeasibleBoundOk = [&](const std::string& type, const HighsInt iX,
                               double& lower, double& upper) {
    double range = upper - lower;
    // Should only be called if lower > upper, so range < 0
    assert(range < 0);
    if (range >= 0) return true;
    if (range > -this->options_.primal_feasibility_tolerance) {
      // Infeasibility is less than feasibility tolerance, so fix
      // bounds at lower (upper) if lower (upper) is an integer - and
      // both can't be integer, otherwise the range <= -1 - otherwise
      // fix at 0.5 * (lower + upper)
      num_ok_infeasible_bound++;
      bool report = num_ok_infeasible_bound <= 10;
      bool integer_lower = lower == std::floor(lower + 0.5);
      bool integer_upper = upper == std::floor(upper + 0.5);
      assert(!integer_lower || !integer_upper);
      if (integer_lower) {
        if (report)
          highsLogUser(log_options, HighsLogType::kInfo,
                       "%s %d bounds [%g, %g] have infeasibility = %g so set "
                       "upper bound to %g\n",
                       type.c_str(), int(iX), lower, upper, range, lower);
        upper = lower;
      } else if (integer_upper) {
        if (report)
          highsLogUser(log_options, HighsLogType::kInfo,
                       "%s %d bounds [%g, %g] have infeasibility = %g so set "
                       "lower bound to %g\n",
                       type.c_str(), int(iX), lower, upper, range, upper);
        lower = upper;
      } else {
        double mid = 0.5 * (lower + upper);
        if (report)
          highsLogUser(log_options, HighsLogType::kInfo,
                       "%s %d bounds [%g, %g] have infeasibility = %g so set "
                       "both bounds to %g\n",
                       type.c_str(), int(iX), lower, upper, range, mid);
        lower = mid;
        upper = mid;
      }
      return true;
    }
    // Infeasibility is greater than feasibility tolerance, so report
    // this (up to 10 times)
    num_true_infeasible_bound++;
    if (num_true_infeasible_bound <= 10)
      highsLogUser(
          log_options, HighsLogType::kInfo,
          "%s %d bounds [%g, %g] have excessive infeasibility = %g%s\n",
          type.c_str(), int(iX), lower, upper, range,
          performed_inward_integer_rounding ? " after inward integer rounding"
                                            : "");
    return false;
  };

  const bool perform_inward_integer_rounding = !this->options_.solve_relaxation;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    performed_inward_integer_rounding = false;
    double lower = lp.col_lower_[iCol];
    double upper = lp.col_upper_[iCol];
    if (has_integrality) {
      // Semi-variables cannot have inconsistent bounds
      if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
          lp.integrality_[iCol] == HighsVarType::kSemiInteger)
        continue;
      if (perform_inward_integer_rounding &&
          lp.integrality_[iCol] == HighsVarType::kInteger) {
        // Assess bounds after inward integer rounding
        double feastol = this->options_.mip_feasibility_tolerance;
        double integer_lower = std::ceil(lower - feastol);
        double integer_upper = std::floor(upper + feastol);
        assert(integer_lower >= lower - feastol);
        assert(integer_upper <= upper + feastol);
        performed_inward_integer_rounding =
            integer_lower > lower || integer_upper < upper;
        lower = integer_lower;
        upper = integer_upper;
      }
    }
    //
    if (lower > upper) {
      if (infeasibleBoundOk("Column", iCol, lower, upper)) {
        // Bound infeasibility is OK (less than the tolerance), so can
        // change the model data
        lp.col_lower_[iCol] = lower;
        lp.col_upper_[iCol] = upper;
      }
    }
    // Note that any inward integer rounding can't be used to change
    // the model data, since it may be a significant change and make
    // the relaxation infeasible when previously it was feasible. In
    // particular, when inward integer rounding leads to inconsistent
    // bounds being propagated to the relaxation, this can prevent a
    // dual ray from being constructed
  }
  // Set performed_inward_integer_rounding = false since it's still
  // checked in infeasibleBoundOk, and may be true from last pass with
  // integer variable
  performed_inward_integer_rounding = false;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (lp.row_lower_[iRow] > lp.row_upper_[iRow])
      infeasibleBoundOk("Row", iRow, lp.row_lower_[iRow], lp.row_upper_[iRow]);
  }
  if (num_ok_infeasible_bound > 0)
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Model has %d small inconsistent bound(s): rectified\n",
                 int(num_ok_infeasible_bound));
  if (num_true_infeasible_bound > 0)
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Model has %d significant inconsistent bound(s): infeasible\n",
                 int(num_true_infeasible_bound));
  return num_true_infeasible_bound > 0;
}

bool Highs::validLinearObjective(const HighsLinearObjective& linear_objective,
                                 const HighsInt iObj) const {
  HighsInt linear_objective_coefficients_size =
      linear_objective.coefficients.size();
  if (linear_objective_coefficients_size != this->model_.lp_.num_col_) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Coefficient vector for linear objective %s has size %d != %d = "
        "lp.num_col_\n",
        iObj >= 0 ? std::to_string(iObj).c_str() : "",
        int(linear_objective_coefficients_size),
        int(this->model_.lp_.num_col_));
    return false;
  }
  if (!options_.blend_multi_objectives &&
      hasRepeatedLinearObjectivePriorities(&linear_objective)) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Repeated priorities for lexicographic optimization is illegal\n");
    return false;
  }
  return true;
}

bool Highs::hasRepeatedLinearObjectivePriorities(
    const HighsLinearObjective* linear_objective) const {
  // Look for repeated values in the linear objective priorities, also
  // comparing linear_objective if it's not a null pointer. Cost is
  // O(n^2), but who will have more than O(1) linear objectives!
  HighsInt num_linear_objective = this->multi_linear_objective_.size();
  if (num_linear_objective <= 0 ||
      (num_linear_objective <= 1 && !linear_objective))
    return false;
  for (HighsInt iObj0 = 0; iObj0 < num_linear_objective; iObj0++) {
    HighsInt priority0 = this->multi_linear_objective_[iObj0].priority;
    for (HighsInt iObj1 = iObj0 + 1; iObj1 < num_linear_objective; iObj1++) {
      HighsInt priority1 = this->multi_linear_objective_[iObj1].priority;
      if (priority1 == priority0) return true;
    }
    if (linear_objective) {
      if (linear_objective->priority == priority0) return true;
    }
  }
  return false;
}

static bool comparison(std::pair<HighsInt, HighsInt> x1,
                       std::pair<HighsInt, HighsInt> x2) {
  return x1.first > x2.first;
}

HighsStatus Highs::returnFromLexicographicOptimization(
    HighsStatus return_status, HighsInt original_lp_num_row) {
  // Save model_status_ and info_ since they are cleared by calling
  // deleteRows
  HighsModelStatus model_status = this->model_status_;
  HighsInfo info = this->info_;
  HighsInt num_linear_objective = this->multi_linear_objective_.size();
  if (num_linear_objective > 1) {
    this->deleteRows(original_lp_num_row, this->model_.lp_.num_row_ - 1);
    // Recover model_status_ and info_, and then account for lack of basis or
    // dual solution
    this->model_status_ = model_status;
    this->info_ = info;
    info_.objective_function_value = 0;
    info_.basis_validity = kBasisValidityInvalid;
    info_.invalidateDualKkt();
    this->solution_.value_valid = true;
    this->model_.lp_.col_cost_.assign(this->model_.lp_.num_col_, 0);
  }
  return return_status;
}

HighsStatus Highs::multiobjectiveSolve() {
  const HighsInt coeff_logging_size_limit = 10;
  HighsInt num_linear_objective = this->multi_linear_objective_.size();

  assert(num_linear_objective > 0);
  HighsLp& lp = this->model_.lp_;
  for (HighsInt iObj = 0; iObj < num_linear_objective; iObj++) {
    HighsLinearObjective& multi_linear_objective =
        this->multi_linear_objective_[iObj];
    if (multi_linear_objective.coefficients.size() !=
        static_cast<size_t>(lp.num_col_)) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Multiple linear objective coefficient vector %d has size "
                   "incompatible with model\n",
                   int(iObj));
      return HighsStatus::kError;
    }
  }

  std::unique_ptr<std::stringstream> multi_objective_log;
  highsLogUser(options_.log_options, HighsLogType::kInfo,
               "Solving with %d multiple linear objectives, %s\n",
               int(num_linear_objective),
               this->options_.blend_multi_objectives
                   ? "blending objectives by weight"
                   : "using lexicographic optimization by priority");
  highsLogUser(
      options_.log_options, HighsLogType::kInfo,
      "Ix      weight      offset     abs_tol     rel_tol    priority%s\n",
      lp.num_col_ < coeff_logging_size_limit ? "   coefficients" : "");
  for (HighsInt iObj = 0; iObj < num_linear_objective; iObj++) {
    HighsLinearObjective& linear_objective =
        this->multi_linear_objective_[iObj];
    multi_objective_log =
        std::unique_ptr<std::stringstream>(new std::stringstream());
    *multi_objective_log << highsFormatToString(
        "%2d %11.6g %11.6g %11.6g %11.6g %11d  ", int(iObj),
        linear_objective.weight, linear_objective.offset,
        linear_objective.abs_tolerance, linear_objective.rel_tolerance,
        int(linear_objective.priority));
    if (lp.num_col_ < coeff_logging_size_limit) {
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
        *multi_objective_log << highsFormatToString(
            "%s c_{%1d} = %g", iCol == 0 ? "" : ",", int(iCol),
            linear_objective.coefficients[iCol]);
    }
    *multi_objective_log << "\n";
    highsLogUser(options_.log_options, HighsLogType::kInfo, "%s",
                 multi_objective_log->str().c_str());
  }
  // Solving with a different objective, but don't call
  // this->clearSolver() since this loses the current solution - that
  // may have been provided by the user (#2419). Just clear the dual
  // data.
  //
  this->clearSolverDualData();
  if (this->options_.blend_multi_objectives) {
    // Objectives are blended by weight and minimized
    lp.offset_ = 0;
    lp.col_cost_.assign(lp.num_col_, 0);
    for (HighsInt iObj = 0; iObj < num_linear_objective; iObj++) {
      HighsLinearObjective& multi_linear_objective =
          this->multi_linear_objective_[iObj];
      lp.offset_ +=
          multi_linear_objective.weight * multi_linear_objective.offset;
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
        lp.col_cost_[iCol] += multi_linear_objective.weight *
                              multi_linear_objective.coefficients[iCol];
    }
    lp.sense_ = ObjSense::kMinimize;

    multi_objective_log =
        std::unique_ptr<std::stringstream>(new std::stringstream());
    *multi_objective_log << highsFormatToString(
        "Solving with blended objective");
    if (lp.num_col_ < coeff_logging_size_limit) {
      *multi_objective_log << highsFormatToString(
          ": %s %g", lp.sense_ == ObjSense::kMinimize ? "min" : "max",
          lp.offset_);
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
        *multi_objective_log << highsFormatToString(
            " + (%g) x[%d]", lp.col_cost_[iCol], int(iCol));
    }
    *multi_objective_log << "\n";
    highsLogUser(options_.log_options, HighsLogType::kInfo, "%s",
                 multi_objective_log->str().c_str());
    return this->optimizeModel();
  }

  // Objectives are applied lexicographically
  if (model_.isQp() && num_linear_objective > 1) {
    // Lexicographic optimization with a single linear objective is
    // trivially standard optimization, so is OK
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Cannot perform non-trivial lexicographic optimization for QP\n");
    return HighsStatus::kError;
  }
  // Check whether there are repeated linear objective priorities
  if (hasRepeatedLinearObjectivePriorities()) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Repeated priorities for lexicographic optimization is illegal\n");
    return HighsStatus::kError;
  }
  std::vector<std::pair<HighsInt, HighsInt>> priority_objective;
  priority_objective.reserve(num_linear_objective);

  for (HighsInt iObj = 0; iObj < num_linear_objective; iObj++)
    priority_objective.push_back(
        std::make_pair(this->multi_linear_objective_[iObj].priority, iObj));
  std::sort(priority_objective.begin(), priority_objective.end(), comparison);
  // Clear LP objective
  lp.offset_ = 0;
  lp.col_cost_.assign(lp.num_col_, 0);
  const HighsInt original_lp_num_row = lp.num_row_;
  std::vector<HighsInt> index(lp.num_col_);
  std::vector<double> value(lp.num_col_);
  // Use the solution of one MIP to provide an integer feasible
  // solution of the next
  HighsSolution solution;
  for (HighsInt iIx = 0; iIx < num_linear_objective; iIx++) {
    HighsInt priority = priority_objective[iIx].first;
    HighsInt iObj = priority_objective[iIx].second;
    // Use this objective
    HighsLinearObjective& linear_objective =
        this->multi_linear_objective_[iObj];
    lp.offset_ = linear_objective.offset;
    lp.col_cost_ = linear_objective.coefficients;
    lp.sense_ =
        linear_objective.weight > 0 ? ObjSense::kMinimize : ObjSense::kMaximize;
    if (lp.isMip() && solution.value_valid) {
      HighsStatus set_solution_status = this->setSolution(solution);
      if (set_solution_status == HighsStatus::kError) {
        highsLogUser(options_.log_options, HighsLogType::kError,
                     "Failure to use one MIP to provide an integer feasible "
                     "solution of the next\n");
        return returnFromLexicographicOptimization(HighsStatus::kError,
                                                   original_lp_num_row);
      }
      bool valid, integral, feasible;
      HighsStatus assess_primal_solution =
          assessPrimalSolution(valid, integral, feasible);
      if (!valid || !integral || !feasible) {
        highsLogUser(options_.log_options, HighsLogType::kWarning,
                     "Failure to use one MIP to provide an integer feasible "
                     "solution of the next: "
                     "status is valid = %s, integral = %s, feasible = %s\n",
                     highsBoolToString(valid).c_str(),
                     highsBoolToString(integral).c_str(),
                     highsBoolToString(feasible).c_str());
      }
    }
    multi_objective_log =
        std::unique_ptr<std::stringstream>(new std::stringstream());
    *multi_objective_log << highsFormatToString("Solving with objective %d",
                                                int(iObj));
    if (lp.num_col_ < coeff_logging_size_limit) {
      *multi_objective_log << highsFormatToString(
          ": %s %g", lp.sense_ == ObjSense::kMinimize ? "min" : "max",
          lp.offset_);
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
        *multi_objective_log << highsFormatToString(
            " + (%g) x[%d]", lp.col_cost_[iCol], int(iCol));
    }
    *multi_objective_log << "\n";
    highsLogUser(options_.log_options, HighsLogType::kInfo, "%s",
                 multi_objective_log->str().c_str());
    HighsStatus optimize_model_status = this->optimizeModel();
    if (optimize_model_status == HighsStatus::kError)
      return returnFromLexicographicOptimization(HighsStatus::kError,
                                                 original_lp_num_row);
    if (model_status_ != HighsModelStatus::kOptimal) {
      highsLogUser(options_.log_options, HighsLogType::kWarning,
                   "After priority %d solve, model status is %s\n",
                   int(priority), modelStatusToString(model_status_).c_str());
      return returnFromLexicographicOptimization(HighsStatus::kWarning,
                                                 original_lp_num_row);
    }
    if (iIx == num_linear_objective - 1) break;
    if (lp.isMip()) {
      // Save the solution to provide an integer feasible solution of
      // the next MIP
      solution.col_value = this->solution_.col_value;
      solution.value_valid = true;
    }
    // Add the constraint
    HighsInt nnz = 0;
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.col_cost_[iCol]) {
        index[nnz] = iCol;
        value[nnz] = lp.col_cost_[iCol];
        nnz++;
      }
    }
    double objective = info_.objective_function_value;
    HighsStatus add_row_status;
    double lower_bound = -kHighsInf;
    double upper_bound = kHighsInf;
    if (lp.sense_ == ObjSense::kMinimize) {
      // Minimizing, so set a greater upper bound than the objective
      if (linear_objective.abs_tolerance >= 0)
        upper_bound = objective + linear_objective.abs_tolerance;
      if (linear_objective.rel_tolerance >= 0) {
        if (objective >= 0) {
          // Guarantees objective of at least (1+t).f^*
          //
          // so ((1+t).f^*-f^*)/f^* = t
          upper_bound = std::min(
              objective * (1.0 + linear_objective.rel_tolerance), upper_bound);
        } else if (objective < 0) {
          // Guarantees objective of at least (1-t).f^*
          //
          // so ((1-t).f^*-f^*)/f^* = -t
          upper_bound = std::min(
              objective * (1.0 - linear_objective.rel_tolerance), upper_bound);
        }
      }
      upper_bound -= lp.offset_;
    } else {
      // Maximizing, so set a lesser lower bound than the objective
      if (linear_objective.abs_tolerance >= 0)
        lower_bound = objective - linear_objective.abs_tolerance;
      if (linear_objective.rel_tolerance >= 0) {
        if (objective >= 0) {
          // Guarantees objective of at most (1-t).f^*
          //
          // so ((1-t).f^*-f^*)/f^* = -t
          lower_bound = std::max(
              objective * (1.0 - linear_objective.rel_tolerance), lower_bound);
        } else if (objective < 0) {
          // Guarantees objective of at least (1+t).f^*
          //
          // so ((1+t).f^*-f^*)/f^* = t
          lower_bound = std::max(
              objective * (1.0 + linear_objective.rel_tolerance), lower_bound);
        }
      }
      lower_bound -= lp.offset_;
    }
    if (lower_bound == -kHighsInf && upper_bound == kHighsInf)
      highsLogUser(options_.log_options, HighsLogType::kWarning,
                   "After priority %d solve, no objective constraint due to "
                   "absolute tolerance being %g < 0,"
                   " and relative tolerance being %g < 0\n",
                   int(priority), linear_objective.abs_tolerance,
                   linear_objective.rel_tolerance);
    multi_objective_log =
        std::unique_ptr<std::stringstream>(new std::stringstream());
    *multi_objective_log << highsFormatToString(
        "Add constraint for objective %d: ", int(iObj));
    if (nnz < coeff_logging_size_limit) {
      *multi_objective_log << highsFormatToString("%g <= ", lower_bound);
      for (HighsInt iEl = 0; iEl < nnz; iEl++)
        *multi_objective_log << highsFormatToString(
            "%s(%g) x[%d]", iEl > 0 ? " + " : "", value[iEl], int(index[iEl]));
      *multi_objective_log << highsFormatToString(" <= %g\n", upper_bound);
    } else {
      *multi_objective_log << highsFormatToString("Bounds [%g, %g]\n",
                                                  lower_bound, upper_bound);
    }
    highsLogUser(options_.log_options, HighsLogType::kInfo, "%s",
                 multi_objective_log->str().c_str());
    add_row_status =
        this->addRow(lower_bound, upper_bound, nnz, index.data(), value.data());
    assert(add_row_status == HighsStatus::kOk);
  }
  return returnFromLexicographicOptimization(HighsStatus::kOk,
                                             original_lp_num_row);
}

bool Highs::tryPdlpCleanup(HighsInt& pdlp_cleanup_iteration_limit,
                           const HighsInfo& presolved_lp_info) const {
  // Primal/dual infeasibilities/residuals can be magnified in
  // postsolve after PDLP, and IPX without crossover can fail,
  // both leading to model_status_ == HighsModelStatus::kUnknown.
  //
  // If the primal/dual infeasibilities/residuals are too large, then it's not
  // worth it, so measure this
  //
  const double tolerance_margin = 1e2;
  bool no_cleanup = false;
  double max_relative_violation = 0;
  // Lambda for updating no_cleanup and max_relative_violation
  auto noCleanup = [&](const std::string& kkt_name, const double kkt_error,
                       const double kkt_tolerance) {
    double use_kkt_tolerance =
        this->options_.kkt_tolerance != kDefaultKktTolerance
            ? this->options_.kkt_tolerance
            : kkt_tolerance;
    double relative_violation = kkt_error / use_kkt_tolerance;
    if (relative_violation > tolerance_margin)
      printf(
          "KKT measure (%11.4g, %11.4g) gives relative violation of %11.4g "
          "for "
          "%s\n",
          kkt_error, use_kkt_tolerance, relative_violation, kkt_name.c_str());
    max_relative_violation =
        std::max(relative_violation, max_relative_violation);
    no_cleanup = max_relative_violation > tolerance_margin;
  };
  noCleanup("Max relative primal infeasibility",
            this->info_.max_relative_primal_infeasibility,
            this->options_.primal_feasibility_tolerance);
  noCleanup("Max relative dual infeasibility",
            this->info_.max_relative_dual_infeasibility,
            this->options_.dual_feasibility_tolerance);
  noCleanup("Max relative primal residual error",
            this->info_.max_relative_primal_residual_error,
            this->options_.primal_residual_tolerance);
  noCleanup("Max relative dual residual error",
            this->info_.max_relative_dual_residual_error,
            this->options_.dual_residual_tolerance);
  noCleanup("Primal-dual objective error",
            this->info_.primal_dual_objective_error,
            this->options_.optimality_tolerance);
  if (no_cleanup) {
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "No PDLP cleanup due to KKT errors exceeding tolerances by a "
                 "max factor = %g > %g = allowed margin\n",
                 max_relative_violation, tolerance_margin);
    return false;
  }
  //
  // Force PDLP to be used with an iteration limit
  if (presolved_lp_info.pdlp_iteration_count > 0) {
    // PDLP was used, so allow 10% of the iterations to clean up
    HighsInt ten_percent_pdlp_iteration_count =
        presolved_lp_info.pdlp_iteration_count / 10;
    pdlp_cleanup_iteration_limit =
        std::max(HighsInt(10000), ten_percent_pdlp_iteration_count);
  } else {
    // IPX without crossover was used, so can only guess what PDLP iteration
    // limit to use
    pdlp_cleanup_iteration_limit = 1000;
  }
  return true;
}

void HighsLinearObjective::clear() {
  this->weight = 0.0;
  this->offset = 0.0;
  this->coefficients.clear();
  this->abs_tolerance = 0.0;
  this->rel_tolerance = 0.0;
  this->priority = 0;
}

void HighsProfiling::initialize(HighsTimer& timer_, const bool sub_solver,
                                const bool mip) {
  // NB this->multi_threaded is set externally:
  //
  // * true in HighsProfiling::clear()
  //
  // * false in initializeSingleThreadedProfiling
  //
  // Hence don't call this->clear() here!
  this->timer = &timer_;
  this->num_profiling_clock_ = kToPresolveSolvePostsolve;
  this->name.assign(this->num_profiling_clock_, "");
  this->name[kPresolveTime] = "Presolve";
  this->name[kSolveTime] = "Solve";
  this->name[kPostsolveTime] = "Postsolve";
  // Now add clocks if performing subsolver profiling
  this->sub_solver_ = sub_solver;
  if (this->sub_solver_) {
    this->num_profiling_clock_ = kToSubSolver;
    this->name.resize(this->num_profiling_clock_);
    this->name[kSubSolverDuSimplexBasis] = "Du simplex (basis)";
    this->name[kSubSolverDuSimplexNoBasis] = "Du simplex (no basis)";
    this->name[kSubSolverPrSimplexBasis] = "Pr simplex (basis)";
    this->name[kSubSolverPrSimplexNoBasis] = "Pr simplex (no basis)";
    this->name[kSubSolverHipo] = "HiPO";
    this->name[kSubSolverIpx] = "IPX";
    this->name[kSubSolverHipoAc] = "HiPO (AC)";
    this->name[kSubSolverIpxAc] = "IPX (AC)";
    this->name[kSubSolverPdlp] = "PDLP";
    this->name[kSubSolverQpAsm] = "QP ASM";
    this->name[kSubSolverMip] = "MIP";
    this->name[kSubSolverSubMip] = "Sub-MIP";
  }
  // Now add clocks if also performing MIP profiling
  // Cannot perform MIP profiling without subsolver profiling
  if (mip) assert(sub_solver);
  this->mip_ = sub_solver && mip;
  if (this->mip_) {
    this->num_profiling_clock_ = kToMipClock;
    this->name.resize(this->num_profiling_clock_);
    initialiseMipProfilingNames(this->name);
  }
  HighsProfilingRecord thread_record;
  thread_record.num_call.assign(this->num_profiling_clock_, 0);
  thread_record.run_time.assign(this->num_profiling_clock_, 0);
  thread_record.start_time.assign(this->num_profiling_clock_, 1);
  HighsInt num_thread = this->numThread();
  assert(num_thread > 0);
  this->submip.assign(num_thread, false);
  this->record.assign(num_thread, thread_record);
  this->submip_record.assign(num_thread, thread_record);
  this->initialized = true;
}

void HighsProfiling::clear() {
  this->timer = nullptr;
  this->multi_threaded = true;
  this->model_name_ = "";
  this->sub_solver_ = false;
  this->mip_ = false;
  this->num_profiling_clock_ = -1;
  this->name.clear();
  this->submip.clear();
  this->record.clear();
  this->submip_record.clear();
  this->initialized = false;
}

HighsInt HighsProfiling::numThread() {
  return this->multi_threaded ? highs::parallel::num_threads() : 1;
}

HighsInt HighsProfiling::myThread() {
  return this->multi_threaded ? highs::parallel::thread_num() : 0;
}

void HighsProfiling::setSubMip(const bool submip) {
  this->submip[this->myThread()] = submip;
}

bool HighsProfiling::isSubMip() { return this->submip[this->myThread()]; }

// Gets the MIP or sub-MIP record according to record_type which, by
// default is kChooseRecord so this->submip is used
HighsProfilingRecord* HighsProfiling::getHighsProfilingRecord(
    const HighsInt record_type) {
  HighsInt thread = this->myThread();
  if (record_type == kSubMipRecord ||
      (record_type == kChooseRecord && this->submip[thread]))
    return &this->submip_record[thread];
  return &this->record[thread];
}

void HighsProfiling::start(const HighsInt profiling_clock, const bool restart) {
  assert(profiling_clock >= 0);
  if (profiling_clock >= this->num_profiling_clock_) return;
  // For a sub-MIP, don't start the clock for anything but a
  // sub-solver
  if (this->isSubMip() && profiling_clock >= kToSubSolver) return;
  // Start timing sub-solver profiling_clock
  HighsInt thread = this->myThread();
  HighsProfilingRecord* thread_record = this->getHighsProfilingRecord();
  double time_start = timer->read();

  if (profiling_clock == kMipClockSubMipSolve) {
    printf("HighsProfiling::start SubMipSolve on thread %2d with submip = %s\n",
           int(thread), this->submip[thread] ? "T" : "F");
  }
  const bool clock_running =
      std::signbit(thread_record->start_time[profiling_clock]);
  if (clock_running && profiling_clock != kSubSolverHipoAc &&
      profiling_clock != kSubSolverIpxAc) {
    // Check for starting a clock that's currently running - except
    // for analytic centre calculation that may by stopped by
    // terminating the task
    printf(
        "HighsProfiling: clock running for thread %d when starting clock %d "
        "(%s) and subMip = %s\n",
        int(thread), int(profiling_clock), this->name[profiling_clock].c_str(),
        this->submip[thread] ? "T" : "F");
    assert(!clock_running);
  }
  thread_record->start_time[profiling_clock] = -time_start;
  if (restart) thread_record->num_call[profiling_clock]--;
}

void HighsProfiling::stop(const HighsInt profiling_clock) {
  assert(profiling_clock >= 0);
  if (profiling_clock >= this->num_profiling_clock_) return;
  // For a sub-MIP, don't start the clock for anything but a
  // sub-solver
  if (this->isSubMip() && profiling_clock >= kToSubSolver) return;
  HighsInt thread = this->myThread();
  HighsProfilingRecord* thread_record = this->getHighsProfilingRecord();
  double time_stop = timer->read();
  double time_start = thread_record->start_time[profiling_clock];
  const bool clock_running = std::signbit(time_start);
  if (!clock_running) {
    // Check for stopping a clock that's currently stopped
    printf(
        "HighsProfiling: clock not running for thread %d when stopping clock "
        "%d (%s) \n",
        int(thread), int(profiling_clock), this->name[profiling_clock].c_str());
    assert(clock_running);
  } else {
    thread_record->num_call[profiling_clock]++;
    thread_record->run_time[profiling_clock] += (time_stop + time_start);
  }
  thread_record->start_time[profiling_clock] = time_stop;
}

double HighsProfiling::read(const HighsInt profiling_clock,
                            const HighsInt record_type) {
  assert(profiling_clock >= 0);
  if (profiling_clock >= this->num_profiling_clock_) return -kHighsInf;
  HighsProfilingRecord* thread_record =
      this->getHighsProfilingRecord(record_type);
  // If the clock is running, work out current running time
  const double current_running_time =
      this->running(profiling_clock, record_type)
          ? thread_record->start_time[profiling_clock] + timer->read()
          : 0;
  return thread_record->run_time[profiling_clock] + current_running_time;
}

HighsInt HighsProfiling::numCall(const HighsInt profiling_clock,
                                 const HighsInt record_type) {
  assert(profiling_clock >= 0);
  if (profiling_clock >= this->num_profiling_clock_) return -kHighsIInf;
  HighsProfilingRecord* thread_record =
      this->getHighsProfilingRecord(record_type);
  const HighsInt this_call_counts =
      this->running(profiling_clock, record_type) ? 1 : 0;
  return thread_record->num_call[profiling_clock] + this_call_counts;
}

bool HighsProfiling::running(const HighsInt profiling_clock,
                             const HighsInt record_type) {
  assert(profiling_clock >= 0);
  if (profiling_clock >= this->num_profiling_clock_) return false;
  HighsProfilingRecord* thread_record =
      this->getHighsProfilingRecord(record_type);
  double time_start = thread_record->start_time[profiling_clock];
  const bool clock_running = std::signbit(time_start);
  return clock_running;
}

void HighsProfiling::solveCall(const std::string& model, const bool submip) {
  const bool printing = false;
  if (this->num_profiling_clock_ <= kToPresolveSolvePostsolve) return;
  const bool local_submip_ok = this->isSubMip() == submip;
  HighsInt thread = this->myThread();
  if (!local_submip_ok) {
    printf("Solving %3s for %4sMIP on thread %d with isSubMip() = %4sMIP\n",
           model.c_str(), submip ? "sub-" : "", int(thread),
           isSubMip() ? "sub-" : "");
  }
  assert(local_submip_ok);
  if (thread != 0 || submip) {
    if (model == "MIP") {
      if (printing)
        printf("Solving MIP for %4sMIP on thread %d\n", submip ? "sub-" : "",
               int(thread));
    } else {
      if (printing)
        printf("Solving %3s for %4sMIP on thread %d\n", model.c_str(),
               submip ? "sub-" : "", int(thread));
    }
  }
}

// HighsInt HighsProfiling::getSepaClockIndex(const std::string& name) {
// assert(1==4);  return 0;}

void Highs::reportProfiling() const {
  if (!this->profiling_->sub_solver_) return;
  HighsInt num_thread = this->profiling_->numThread();
  double mip_time = 0;
  double max_sumip_time = 0;
  const std::vector<HighsProfilingRecord>& record = this->profiling_->record;
  const std::vector<HighsProfilingRecord>& submip_record =
      this->profiling_->submip_record;
  for (HighsInt thread_num = 0; thread_num < num_thread; thread_num++) {
    mip_time = std::max(record[thread_num].run_time[kSubSolverMip], mip_time);
    max_sumip_time =
        std::max(record[thread_num].run_time[kSubSolverSubMip], max_sumip_time);
  }

  std::vector<HighsInt> used_thread;
  for (HighsInt thread_num = 0; thread_num < num_thread; thread_num++) {
    bool used = false;
    for (HighsInt Ix = kFromSubSolver; Ix < kToSubSolver; Ix++)
      if (record[thread_num].num_call[Ix]) used = true;
    for (HighsInt Ix = kFromSubSolver; Ix < kToSubSolver; Ix++)
      if (submip_record[thread_num].num_call[Ix]) used = true;
    if (!used) continue;
    used_thread.push_back(thread_num);
  }
  const double num_threads_used = used_thread.size();
  std::stringstream ss;
  std::vector<HighsBool> mip_used_sub_solver(kToSubSolver, false);
  std::vector<HighsBool> submip_used_sub_solver(kToSubSolver, false);
  const HighsInt to_k = max_sumip_time > 0 ? 2 : 1;
  const std::vector<std::string>& name = this->profiling_->name;
  double sum_sum_mip_sub_solve_time = 0;
  for (HighsInt k = 0; k < to_k; k++) {
    if (k == 0) {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "\nMIP sub-solver profiling: number of threads used = %d\n",
                   int(num_threads_used));
    } else {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "\nSub-MIP sub-solver profiling\n");
    }
    for (HighsInt thread_ix = 0; thread_ix < HighsInt(num_threads_used);
         thread_ix++) {
      HighsInt thread_num = used_thread[thread_ix];
      double ideal_time =
          k == 0
              ? mip_time
              : this->profiling_->record[thread_num].run_time[kSubSolverSubMip];
      if (ideal_time <= 0) continue;
      const std::vector<HighsProfilingRecord>& record =
          k == 0 ? this->profiling_->record : this->profiling_->submip_record;
      std::vector<HighsBool>& used_sub_solver =
          k == 0 ? mip_used_sub_solver : submip_used_sub_solver;
      const std::vector<HighsInt>& num_call = record[thread_num].num_call;
      const std::vector<double>& run_time = record[thread_num].run_time;
      ss.str(std::string());
      ss << highsFormatToString(
          "\nThread %d\n"
          "Solver                    Calls    Time       "
          "Time/call",
          int(thread_num));
      if (ideal_time > 0) ss << (k == 0 ? "      MIP%" : "  Sub-MIP%");
      highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                   ss.str().c_str());
      double sum_mip_sub_solve_time = 0;
      for (HighsInt Ix = kFromSubSolver; Ix < kToSubSolver; Ix++) {
        double pct = 0;
        if (!num_call[Ix]) continue;
        used_sub_solver[Ix] = true;
        ss.str(std::string());
        ss << highsFormatToString("%-21s %9d %11.4e %11.4e", name[Ix].c_str(),
                                  int(num_call[Ix]), run_time[Ix],
                                  run_time[Ix] / (1.0 * num_call[Ix]));
        if (ideal_time > 0 && Ix != kSubSolverMip) {
          sum_mip_sub_solve_time += run_time[Ix];
          ss << highsFormatToString("     %5.1f",
                                    1e2 * run_time[Ix] / ideal_time);
        }
        highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                     ss.str().c_str());
      }
      sum_sum_mip_sub_solve_time += sum_mip_sub_solve_time;
      if (ideal_time > 0 && sum_mip_sub_solve_time > 0)
        highsLogUser(
            options_.log_options, HighsLogType::kInfo,
            "TOTAL                           %11.4e                 %5.1f\n",
            sum_mip_sub_solve_time, 1e2 * sum_mip_sub_solve_time / ideal_time);
    }
  }
  if (mip_time <= 0) return;
  if (sum_sum_mip_sub_solve_time <= 0) return;
  // Lambda for horizontal rule
  auto hrule = [&]() {
    ss.str(std::string());
    ss << "=====================";
    for (HighsInt thread_ix = 0; thread_ix < HighsInt(num_threads_used);
         thread_ix++)
      ss << "======";
    highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                 ss.str().c_str());
  };
  // Determine the sub-solver percentage breakdown over all threads
  // when solving MIPs
  highsLogUser(options_.log_options, HighsLogType::kInfo,
               "\nPercent (sub-)MIP time by thread\n");
  for (HighsInt k = 0; k < to_k; k++) {
    if (k == 0) {
      ss.str(std::string());
      ss << highsFormatToString("\nMIP sub-solver       ");
    } else {
      if (max_sumip_time <= 0) continue;
      ss.str(std::string());
      ss << highsFormatToString("\nSub-MIP sub-solver   ");
    }
    for (HighsInt thread_ix = 0; thread_ix < HighsInt(num_threads_used);
         thread_ix++) {
      HighsInt thread_num = used_thread[thread_ix];
      ss << highsFormatToString("%6d", int(thread_num));
    }
    highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                 ss.str().c_str());
    std::vector<HighsBool>& used_sub_solver =
        k == 0 ? mip_used_sub_solver : submip_used_sub_solver;
    const std::vector<HighsProfilingRecord>& record =
        k == 0 ? this->profiling_->record : this->profiling_->submip_record;
    std::vector<double> totalPct(num_threads_used, 0);
    for (HighsInt Ix = kFromSubSolver + 1; Ix < kToSubSolver; Ix++) {
      if (!used_sub_solver[Ix]) continue;
      ss.str(std::string());
      ss << highsFormatToString("%-21s", name[Ix].c_str());
      for (HighsInt thread_ix = 0; thread_ix < HighsInt(num_threads_used);
           thread_ix++) {
        HighsInt thread_num = used_thread[thread_ix];
        double ideal_time = k == 0 ? mip_time
                                   : this->profiling_->record[thread_num]
                                         .run_time[kSubSolverSubMip];
        HighsInt num_call = record[thread_num].num_call[Ix];
        double run_time = record[thread_num].run_time[Ix];
        if (num_call && ideal_time > 0) {
          double pct = 1e2 * run_time / ideal_time;
          totalPct[thread_ix] += pct;
          ss << highsFormatToString(" %5.1f", pct);
        } else {
          ss << "      ";
        }
      }
      highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                   ss.str().c_str());
    }
    hrule();
    ss.str(std::string());
    ss << "Total                ";
    for (HighsInt thread_ix = 0; thread_ix < HighsInt(num_threads_used);
         thread_ix++) {
      if (totalPct[thread_ix]) {
        ss << highsFormatToString(" %5.1f", totalPct[thread_ix]);
      } else {
        ss << "      ";
      }
    }
    highsLogUser(options_.log_options, HighsLogType::kInfo, "%s\n",
                 ss.str().c_str());
    hrule();
  }
}

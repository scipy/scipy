/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkk.cpp
 * @brief
 */
#include "simplex/HEkk.h"

#include "lp_data/HighsLpSolverObject.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolutionDebug.h"
#include "parallel/HighsParallel.h"
#include "simplex/HEkkDual.h"
#include "simplex/HEkkPrimal.h"
#include "simplex/HSimplexDebug.h"
#include "simplex/HSimplexReport.h"
#include "simplex/SimplexTimer.h"

using std::fabs;
using std::max;
using std::min;

// using std::cout;
// using std::endl;

void HEkk::clear() {
  // Clears Ekk entirely. Clears all associated pointers, data scalars
  // and vectors, and the status values.
  this->clearEkkLp();
  this->clearEkkDualize();
  this->clearEkkData();
  this->clearEkkDualEdgeWeightData();
  this->clearEkkPointers();
  this->basis_.clear();
  this->simplex_nla_.clear();
  this->clearEkkAllStatus();
  this->clearRayRecords();
}

void HEkk::clearEkkAllStatus() {
  // Clears the Ekk status entirely. Functionally junks all
  // information relating to the simplex solve, but doesn't clear the
  // associated data scalars and vectors.
  HighsSimplexStatus& status = this->status_;
  status.initialised_for_new_lp = false;
  status.initialised_for_solve = false;
  this->clearNlaStatus();
  this->clearEkkDataStatus();
}

void HEkk::clearEkkDataStatus() {
  // Just clears the Ekk status values associated with Ekk-specific
  // data: doesn't clear "initialised_for_new_lp", "initialised_for_solve" or
  // NLA status
  HighsSimplexStatus& status = this->status_;
  status.has_ar_matrix = false;
  status.has_dual_steepest_edge_weights = false;
  status.has_fresh_rebuild = false;
  status.has_dual_objective_value = false;
  status.has_primal_objective_value = false;
}

void HEkk::clearNlaStatus() {
  // Clears Ekk status values associated with NLA. Functionally junks
  // NLA, but doesn't clear the associated data scalars and vectors
  HighsSimplexStatus& status = this->status_;
  status.has_basis = false;
  status.has_nla = false;
  clearNlaInvertStatus();
}

void HEkk::clearNlaInvertStatus() {
  this->status_.has_invert = false;
  this->status_.has_fresh_invert = false;
}

void HEkk::clearRayRecords() {
  this->dual_ray_record_.clear();
  this->primal_ray_record_.clear();
}

void HEkk::clearEkkPointers() {
  this->callback_ = nullptr;
  this->options_ = nullptr;
  this->timer_ = nullptr;
}

void HEkk::clearEkkLp() {
  this->lp_.clear();
  lp_name_ = "";
}

void HEkk::clearEkkDualize() {
  this->original_col_cost_.clear();
  this->original_col_lower_.clear();
  this->original_col_upper_.clear();
  this->original_row_lower_.clear();
  this->original_row_upper_.clear();
  this->upper_bound_col_.clear();
  this->upper_bound_row_.clear();
}

void HEkk::clearEkkDualEdgeWeightData() {
  this->dual_edge_weight_.clear();
  this->scattered_dual_edge_weight_.clear();
}

void HEkk::clearEkkData() {
  // Clears Ekk-specific data scalars and vectors. Doesn't clear
  // status, as this is done elsewhere
  //
  // Doesn't clear the LP, simplex NLA or simplex basis as part of
  // clearing Ekk data, so that the simplex basis and HFactor instance
  // are maintained
  //
  // analysis_; No clear yet

  this->clearEkkDataInfo();
  model_status_ = HighsModelStatus::kNotset;
  // random_; Has no data

  this->simplex_in_scaled_space_ = false;
  this->ar_matrix_.clear();
  this->scaled_a_matrix_.clear();

  this->cost_scale_ = 1;
  this->iteration_count_ = 0;
  this->dual_simplex_cleanup_level_ = 0;
  this->dual_simplex_phase1_cleanup_level_ = 0;

  this->previous_iteration_cycling_detected = -kHighsIInf;

  this->solve_bailout_ = false;
  this->called_return_from_solve_ = false;
  this->exit_algorithm_ = SimplexAlgorithm::kNone;
  this->return_primal_solution_status_ = 0;
  this->return_dual_solution_status_ = 0;

  this->proof_index_.clear();
  this->proof_value_.clear();

  this->clearRayRecords();

  this->build_synthetic_tick_ = 0.0;
  this->total_synthetic_tick_ = 0.0;

  // Clear values used for debugging
  this->debug_solve_call_num_ = 0;
  this->debug_basis_id_ = 0;
  this->time_report_ = false;
  this->debug_initial_build_synthetic_tick_ = 0;
  this->debug_solve_report_ = false;
  this->debug_iteration_report_ = false;
  this->debug_basis_report_ = false;
  this->debug_dual_feasible = false;
  this->debug_max_relative_dual_steepest_edge_weight_error = 0;

  clearBadBasisChange();

  this->primal_phase1_dual_.clear();
}

void HEkk::clearEkkDataInfo() {
  HighsSimplexInfo& info = this->info_;
  info.workCost_.clear();
  info.workDual_.clear();
  info.workShift_.clear();
  info.workLower_.clear();
  info.workUpper_.clear();
  info.workRange_.clear();
  info.workValue_.clear();
  info.workLowerShift_.clear();
  info.workUpperShift_.clear();
  info.baseLower_.clear();
  info.baseUpper_.clear();
  info.baseValue_.clear();
  info.numTotRandomValue_.clear();
  info.numTotPermutation_.clear();
  info.numColPermutation_.clear();
  info.devex_index_.clear();
  info.pivot_.clear();
  info.index_chosen_.clear();
  info.phase1_backtracking_test_done = false;
  info.phase2_backtracking_test_done = false;
  info.backtracking_ = false;
  info.valid_backtracking_basis_ = false;
  info.backtracking_basis_.clear();
  info.backtracking_basis_costs_shifted_ = false;
  info.backtracking_basis_costs_perturbed_ = false;
  info.backtracking_basis_bounds_shifted_ = false;
  info.backtracking_basis_bounds_perturbed_ = false;
  info.backtracking_basis_workShift_.clear();
  info.backtracking_basis_workLowerShift_.clear();
  info.backtracking_basis_workUpperShift_.clear();
  info.backtracking_basis_edge_weight_.clear();
  info.simplex_strategy = 0;
  info.dual_edge_weight_strategy = 0;
  info.primal_edge_weight_strategy = 0;
  info.price_strategy = 0;
  info.dual_simplex_cost_perturbation_multiplier = 1;
  info.primal_simplex_phase1_cost_perturbation_multiplier = 1;
  info.primal_simplex_bound_perturbation_multiplier = 1;
  info.allow_dual_steepest_edge_to_devex_switch = 0;
  info.dual_steepest_edge_weight_log_error_threshold = 0;
  info.run_quiet = false;
  info.store_squared_primal_infeasibility = false;
  info.report_simplex_inner_clock = false;
  info.report_simplex_outer_clock = false;
  info.report_simplex_phases_clock = false;
  info.report_HFactor_clock = false;
  info.analyse_lp = false;
  info.analyse_iterations = false;
  info.analyse_invert_form = false;

  info.allow_cost_shifting = true;
  info.allow_cost_perturbation = true;
  info.allow_bound_perturbation = true;
  info.costs_shifted = false;
  info.costs_perturbed = false;
  info.bounds_shifted = false;
  info.bounds_perturbed = false;

  info.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
  info.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
  info.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  info.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
  info.dual_phase1_iteration_count = 0;
  info.dual_phase2_iteration_count = 0;
  info.primal_phase1_iteration_count = 0;
  info.primal_phase2_iteration_count = 0;
  info.primal_bound_swap = 0;
  info.min_concurrency = 1;
  info.num_concurrency = 1;
  info.max_concurrency = kSimplexConcurrencyLimit;
  info.multi_iteration = 0;
  info.update_count = 0;
  info.dual_objective_value = 0;
  info.primal_objective_value = 0;
  info.updated_dual_objective_value = 0;
  info.updated_primal_objective_value = 0;
  info.num_basic_logicals = 0;
}

void HEkk::clearEkkControlInfo() {
  HighsSimplexInfo& info = this->info_;
  info.control_iteration_count0 = 0;
  info.col_aq_density = 0.0;
  info.row_ep_density = 0.0;
  info.row_ap_density = 0.0;
  info.row_DSE_density = 0.0;
  info.col_steepest_edge_density = 0.0;
  info.col_basic_feasibility_change_density = 0.0;
  info.row_basic_feasibility_change_density = 0.0;
  info.col_BFRT_density = 0.0;
  info.primal_col_density = 0.0;
  info.dual_col_density = 0.0;
  info.costly_DSE_frequency = 0;
  info.num_costly_DSE_iteration = 0;
  info.costly_DSE_measure = 0;
  info.average_log_low_DSE_weight_error = 0;
  info.average_log_high_DSE_weight_error = 0;
}

void HEkk::clearEkkNlaInfo() {
  HighsSimplexInfo& info = this->info_;
  info.factor_pivot_threshold = 0;
  info.update_limit = 0;
}

void HEkk::invalidate() {
  this->status_.initialised_for_new_lp = false;
  assert(!this->status_.is_dualized);
  assert(!this->status_.is_permuted);
  this->status_.initialised_for_solve = false;
  this->invalidateBasisMatrix();
  this->simplex_stats_.initialise();
}

void HEkk::invalidateBasisMatrix() {
  // When the constraint matrix changes - dimensions or just (basic)
  // values, the simplex NLA becomes invalid, the simplex basis is
  // no longer valid, and
  this->status_.has_nla = false;
  invalidateBasis();
}

void HEkk::invalidateBasis() {
  // Invalidate the basis of the simplex LP, and all its other
  // basis-related properties
  this->status_.has_basis = false;
  this->invalidateBasisArtifacts();
}

void HEkk::invalidateBasisArtifacts() {
  // Invalidate the artifacts of the basis of the simplex LP
  this->status_.has_ar_matrix = false;
  this->status_.has_dual_steepest_edge_weights = false;
  this->status_.has_invert = false;
  this->status_.has_fresh_invert = false;
  this->status_.has_fresh_rebuild = false;
  this->status_.has_dual_objective_value = false;
  this->status_.has_primal_objective_value = false;
  this->clearRayRecords();
}

void HEkk::updateStatus(LpAction action) {
  assert(!this->status_.is_dualized);
  assert(!this->status_.is_permuted);
  switch (action) {
    case LpAction::kScale:
      this->invalidateBasisMatrix();
      break;
    case LpAction::kNewCosts:
      this->status_.has_fresh_rebuild = false;
      this->status_.has_dual_objective_value = false;
      this->status_.has_primal_objective_value = false;
      break;
    case LpAction::kNewBounds:
      this->status_.has_fresh_rebuild = false;
      this->status_.has_dual_objective_value = false;
      this->status_.has_primal_objective_value = false;
      break;
    case LpAction::kNewBasis:
      this->invalidateBasis();
      break;
    case LpAction::kNewCols:
      this->clear();
      //    this->invalidateBasisArtifacts();
      break;
    case LpAction::kNewRows:
      if (kExtendInvertWhenAddingRows) {
        // Just clear Ekk data
        this->clearEkkData();
      } else {
        // Clear everything
        this->clear();
      }
      //    this->invalidateBasisArtifacts();
      break;
    case LpAction::kDelCols:
      this->clear();
      //    this->invalidateBasis();
      break;
    case LpAction::kDelNonbasicCols:
      this->clear();
      //    this->invalidateBasis();
      break;
    case LpAction::kDelRows:
      this->clear();
      //   this->invalidateBasis();
      break;
    case LpAction::kDelRowsBasisOk:
      assert(1 == 0);
      //      info.lp_ = true;
      break;
    case LpAction::kScaledCol:
      this->invalidateBasisMatrix();
      break;
    case LpAction::kScaledRow:
      this->invalidateBasisMatrix();
      break;
    case LpAction::kBacktracking:
      this->status_.has_ar_matrix = false;
      this->status_.has_fresh_rebuild = false;
      this->status_.has_dual_objective_value = false;
      this->status_.has_primal_objective_value = false;
      break;
    default:
      break;
  }
}

void HEkk::setNlaPointersForLpAndScale(const HighsLp& lp) {
  assert(status_.has_nla);
  simplex_nla_.setLpAndScalePointers(&lp);
}

void HEkk::setNlaPointersForTrans(const HighsLp& lp) {
  assert(status_.has_nla);
  assert(status_.has_basis);
  simplex_nla_.setLpAndScalePointers(&lp);
  simplex_nla_.basic_index_ = basis_.basicIndex_.data();
}

void HEkk::setNlaRefactorInfo() {
  simplex_nla_.factor_.refactor_info_ = this->hot_start_.refactor_info;
  simplex_nla_.factor_.refactor_info_.use = true;
}

void HEkk::btran(HVector& rhs, const double expected_density) {
  assert(status_.has_nla);
  simplex_nla_.btran(rhs, expected_density);
}

void HEkk::ftran(HVector& rhs, const double expected_density) {
  assert(status_.has_nla);
  simplex_nla_.ftran(rhs, expected_density);
}

void HEkk::moveLp(HighsLpSolverObject& solver_object) {
  // Move the incumbent LP to EKK
  HighsLp& incumbent_lp = solver_object.lp_;
  this->lp_ = std::move(incumbent_lp);
  incumbent_lp.is_moved_ = true;
  //
  // Invalidate the row-wise matrix
  this->status_.has_ar_matrix = false;
  //
  // The simplex algorithm runs in the same space as the LP that has
  // just been moved in. This is a scaled space if the LP is scaled.
  this->simplex_in_scaled_space_ = this->lp_.is_scaled_;
  //
  // Update other EKK pointers. Currently just pointers to the
  // HighsOptions and HighsTimer members of the Highs class that are
  // communicated by reference via the HighsLpSolverObject instance.
  this->setPointers(&solver_object.callback_, &solver_object.options_,
                    &solver_object.timer_);
  // Initialise Ekk if this has not been done. Ekk isn't initialised
  // if moveLp hasn't been called for this instance of HiGHS, or if
  // the Ekk instance is junked due to removing rows from the LP
  this->initialiseEkk();
}

void HEkk::setPointers(HighsCallback* callback, HighsOptions* options,
                       HighsTimer* timer) {
  this->callback_ = callback;
  this->options_ = options;
  this->timer_ = timer;
  this->analysis_.timer_ = this->timer_;
}

HighsSparseMatrix* HEkk::getScaledAMatrixPointer() {
  // Return a pointer to either the constraint matrix or a scaled copy
  // (that is a member of the HEkk class), with the latter returned if
  // the LP has scaling factors but is unscaled.
  HighsSparseMatrix* local_scaled_a_matrix = &(this->lp_.a_matrix_);
  if (this->lp_.scale_.has_scaling && !this->lp_.is_scaled_) {
    scaled_a_matrix_ = this->lp_.a_matrix_;
    scaled_a_matrix_.applyScale(this->lp_.scale_);
    local_scaled_a_matrix = &scaled_a_matrix_;
  }
  return local_scaled_a_matrix;
}

HighsStatus HEkk::dualize() {
  assert(lp_.a_matrix_.isColwise());
  original_num_col_ = lp_.num_col_;
  original_num_row_ = lp_.num_row_;
  original_num_nz_ = lp_.numNz();
  original_offset_ = lp_.offset_;
  original_col_cost_ = lp_.col_cost_;
  original_col_lower_ = lp_.col_lower_;
  original_col_upper_ = lp_.col_upper_;
  original_row_lower_ = lp_.row_lower_;
  original_row_upper_ = lp_.row_upper_;
  // Reserve space for simple dual
  lp_.col_cost_.reserve(original_num_row_);
  lp_.col_lower_.reserve(original_num_row_);
  lp_.col_upper_.reserve(original_num_row_);
  lp_.row_lower_.reserve(original_num_col_);
  lp_.row_upper_.reserve(original_num_col_);
  // Invalidate the original data
  lp_.col_cost_.resize(0);
  lp_.col_lower_.resize(0);
  lp_.col_upper_.resize(0);
  lp_.row_lower_.resize(0);
  lp_.row_upper_.resize(0);
  // The bulk of the constraint matrix of the dual LP is the transpose
  // of the primal constraint matrix. This is obtained row-wise by
  // copying the matrix and flipping the dimensions
  HighsSparseMatrix dual_matrix = lp_.a_matrix_;
  dual_matrix.num_row_ = original_num_col_;
  dual_matrix.num_col_ = original_num_row_;
  dual_matrix.format_ = MatrixFormat::kRowwise;
  // The primal_bound_value vector accumulates the values of the
  // finite bounds on variables - or zero for a free variable - used
  // later to compute the offset and shift for the costs. Many of
  // these components will be zero - all for the case x>=0 - so
  // maintain a list of the nonzeros and corresponding indices. Don't
  // reserve space since they may not be needed
  vector<double> primal_bound_value;
  vector<HighsInt> primal_bound_index;
  const double inf = kHighsInf;
  for (HighsInt iCol = 0; iCol < original_num_col_; iCol++) {
    const double cost = original_col_cost_[iCol];
    const double lower = original_col_lower_[iCol];
    const double upper = original_col_upper_[iCol];
    double primal_bound = inf;
    double row_lower = inf;
    double row_upper = -inf;
    if (lower == upper) {
      // Fixed
      primal_bound = lower;
      // Dual activity a^Ty is free, implying dual for primal column
      // (slack for dual row) is free
      row_lower = -inf;
      row_upper = inf;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        //
        // Treat as lower
        primal_bound = lower;
        // Dual activity a^Ty is bounded above by cost, implying dual
        // for primal column (slack for dual row) is non-negative
        row_lower = -inf;
        row_upper = cost;
        // Treat upper bound as additional constraint
        upper_bound_col_.push_back(iCol);
      } else {
        // Lower (since upper bound is infinite)
        primal_bound = lower;
        // Dual activity a^Ty is bounded above by cost, implying dual
        // for primal column (slack for dual row) is non-negative
        row_lower = -inf;
        row_upper = cost;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      primal_bound = upper;
      // Dual activity a^Ty is bounded below by cost, implying dual
      // for primal column (slack for dual row) is non-positive
      row_lower = cost;
      row_upper = inf;
    } else {
      // FREE
      //
      // Dual activity a^Ty is fixed by cost, implying dual for primal
      // column (slack for dual row) is fixed at zero
      primal_bound = 0;
      row_lower = cost;
      row_upper = cost;
    }
    assert(row_lower < inf);
    assert(row_upper > -inf);
    assert(primal_bound < inf);
    lp_.row_lower_.push_back(row_lower);
    lp_.row_upper_.push_back(row_upper);
    if (primal_bound) {
      primal_bound_value.push_back(primal_bound);
      primal_bound_index.push_back(iCol);
    }
  }
  for (HighsInt iRow = 0; iRow < original_num_row_; iRow++) {
    double lower = original_row_lower_[iRow];
    double upper = original_row_upper_[iRow];
    double col_cost = inf;
    double col_lower = inf;
    double col_upper = -inf;
    if (lower == upper) {
      // Equality constraint
      //
      // Dual variable has primal RHS as cost and is free
      col_cost = lower;
      col_lower = -inf;
      col_upper = inf;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        //
        // Treat as lower
        col_cost = lower;
        col_lower = 0;
        col_upper = inf;
        // Treat upper bound as additional constraint
        upper_bound_row_.push_back(iRow);
      } else {
        // Lower (since upper bound is infinite)
        col_cost = lower;
        col_lower = 0;
        col_upper = inf;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      col_cost = upper;
      col_lower = -inf;
      col_upper = 0;
    } else {
      // FREE
      // Shouldn't get free rows, but handle them anyway
      col_cost = 0;
      col_lower = 0;
      col_upper = 0;
    }
    assert(col_lower < inf);
    assert(col_upper > -inf);
    assert(col_cost < inf);
    lp_.col_cost_.push_back(col_cost);
    lp_.col_lower_.push_back(col_lower);
    lp_.col_upper_.push_back(col_upper);
  }
  vector<HighsInt>& start = lp_.a_matrix_.start_;
  vector<HighsInt>& index = lp_.a_matrix_.index_;
  vector<double>& value = lp_.a_matrix_.value_;
  // Boxed variables and constraints yield extra columns in the dual LP
  HighsSparseMatrix extra_columns;
  extra_columns.ensureColwise();
  extra_columns.num_row_ = original_num_col_;
  HighsInt num_upper_bound_col = upper_bound_col_.size();
  HighsInt num_upper_bound_row = upper_bound_row_.size();
  double one = 1;
  for (HighsInt iX = 0; iX < num_upper_bound_col; iX++) {
    HighsInt iCol = upper_bound_col_[iX];
    const double upper = original_col_upper_[iCol];
    extra_columns.addVec(1, &iCol, &one);
    lp_.col_cost_.push_back(upper);
    lp_.col_lower_.push_back(-inf);
    lp_.col_upper_.push_back(0);
  }

  if (num_upper_bound_row) {
    // Need to identify the submatrix of constraint matrix rows
    // corresponding to those with a row index in
    // upper_bound_row_. When identifying numbers of entries in each
    // row of submatrix, use indirection to get corresponding row
    // index, with a dummy row for rows not in the submatrix.
    HighsInt dummy_row = num_upper_bound_row;
    vector<HighsInt> indirection;
    vector<HighsInt> count;
    indirection.assign(original_num_row_, dummy_row);
    count.assign(num_upper_bound_row + 1, 0);
    HighsInt extra_iRow = 0;
    for (HighsInt iX = 0; iX < num_upper_bound_row; iX++) {
      HighsInt iRow = upper_bound_row_[iX];
      indirection[iRow] = extra_iRow++;
      double upper = original_row_upper_[iRow];
      lp_.col_cost_.push_back(upper);
      lp_.col_lower_.push_back(-inf);
      lp_.col_upper_.push_back(0);
    }
    for (HighsInt iEl = 0; iEl < original_num_nz_; iEl++)
      count[indirection[index[iEl]]]++;
    extra_columns.start_.resize(num_upper_bound_col + num_upper_bound_row + 1);
    for (HighsInt iRow = 0; iRow < num_upper_bound_row; iRow++) {
      extra_columns.start_[num_upper_bound_col + iRow + 1] =
          extra_columns.start_[num_upper_bound_col + iRow] + count[iRow];
      count[iRow] = extra_columns.start_[num_upper_bound_col + iRow];
    }
    HighsInt extra_columns_num_nz =
        extra_columns.start_[num_upper_bound_col + num_upper_bound_row];
    extra_columns.index_.resize(extra_columns_num_nz);
    extra_columns.value_.resize(extra_columns_num_nz);
    for (HighsInt iCol = 0; iCol < original_num_col_; iCol++) {
      for (HighsInt iEl = start[iCol]; iEl < start[iCol + 1]; iEl++) {
        HighsInt iRow = indirection[index[iEl]];
        if (iRow < num_upper_bound_row) {
          HighsInt extra_columns_iEl = count[iRow];
          assert(extra_columns_iEl < extra_columns_num_nz);
          extra_columns.index_[extra_columns_iEl] = iCol;
          extra_columns.value_[extra_columns_iEl] = value[iEl];
          count[iRow]++;
        }
      }
    }
    extra_columns.num_col_ += num_upper_bound_row;
  }
  // Incorporate the cost shift by subtracting A*primal_bound from the
  // cost vector; compute the objective offset
  double delta_offset = 0;
  for (size_t iX = 0; iX < primal_bound_index.size(); iX++) {
    HighsInt iCol = primal_bound_index[iX];
    double multiplier = primal_bound_value[iX];
    delta_offset += multiplier * original_col_cost_[iCol];
    for (HighsInt iEl = start[iCol]; iEl < start[iCol + 1]; iEl++)
      lp_.col_cost_[index[iEl]] -= multiplier * value[iEl];
  }
  if (extra_columns.num_col_) {
    // Incorporate the cost shift by subtracting
    // extra_columns*primal_bound from the cost vector for the extra
    // dual variables
    //
    // Have to scatter the packed primal bound values into a
    // full-length vector
    //
    // ToDo Make this more efficient?
    vector<double> primal_bound;
    primal_bound.assign(original_num_col_, 0);
    for (size_t iX = 0; iX < primal_bound_index.size(); iX++)
      primal_bound[primal_bound_index[iX]] = primal_bound_value[iX];

    for (HighsInt iCol = 0; iCol < extra_columns.num_col_; iCol++) {
      double cost = lp_.col_cost_[original_num_row_ + iCol];
      for (HighsInt iEl = extra_columns.start_[iCol];
           iEl < extra_columns.start_[iCol + 1]; iEl++)
        cost -=
            primal_bound[extra_columns.index_[iEl]] * extra_columns.value_[iEl];
      lp_.col_cost_[original_num_row_ + iCol] = cost;
    }
  }
  lp_.offset_ += delta_offset;
  // Copy the row-wise dual LP constraint matrix and transpose it.
  lp_.a_matrix_ = dual_matrix.getColwise();
  assert(lp_.a_matrix_.isColwise());
  // lp_.a_matrix_.ensureColwise(); is now redundant

  // Add the extra columns to the dual LP constraint matrix
  lp_.a_matrix_.addCols(extra_columns);

  HighsInt dual_num_col =
      original_num_row_ + num_upper_bound_col + num_upper_bound_row;
  HighsInt dual_num_row = original_num_col_;
  assert(dual_num_col == (int)lp_.col_cost_.size());
  assert(lp_.a_matrix_.num_col_ == dual_num_col);
  const bool ignore_scaling = true;
  if (!ignore_scaling) {
    // Flip any scale factors
    if (lp_.scale_.has_scaling) {
      std::vector<double> temp_scale = lp_.scale_.row;
      lp_.scale_.row = lp_.scale_.col;
      lp_.scale_.col = temp_scale;
      lp_.scale_.num_col = dual_num_col;
      lp_.scale_.num_row = dual_num_row;
    }
  }
  // Change optimization sense
  if (lp_.sense_ == ObjSense::kMinimize) {
    lp_.sense_ = ObjSense::kMaximize;
  } else {
    lp_.sense_ = ObjSense::kMinimize;
  }
  // Flip LP dimensions
  lp_.num_col_ = dual_num_col;
  lp_.num_row_ = dual_num_row;
  status_.is_dualized = true;
  status_.has_basis = false;
  status_.has_ar_matrix = false;
  status_.has_nla = false;
  highsLogUser(options_->log_options, HighsLogType::kInfo,
               "Solving dual LP with %d columns", (int)dual_num_col);
  if (num_upper_bound_col + num_upper_bound_row) {
    highsLogUser(options_->log_options, HighsLogType::kInfo, " [%d extra from",
                 (int)dual_num_col - original_num_row_);
    if (num_upper_bound_col)
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   " %d boxed variable(s)", (int)num_upper_bound_col);
    if (num_upper_bound_col && num_upper_bound_row)
      highsLogUser(options_->log_options, HighsLogType::kInfo, " and");
    if (num_upper_bound_row)
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   " %d boxed constraint(s)", (int)num_upper_bound_row);
    highsLogUser(options_->log_options, HighsLogType::kInfo, "]");
  }
  highsLogUser(options_->log_options, HighsLogType::kInfo, " and %d rows\n",
               (int)dual_num_row);
  //  reportLp(options_->log_options, lp_, HighsLogType::kVerbose);
  return HighsStatus::kOk;
}

HighsStatus HEkk::undualize() {
  if (!this->status_.is_dualized) return HighsStatus::kOk;
  HighsInt dual_num_col = lp_.num_col_;
  HighsInt primal_num_tot = original_num_col_ + original_num_row_;
  // These two aren't used (yet)
  // vector<double>& dual_work_dual = info_.workDual_;
  // vector<double>& primal_work_value = info_.workValue_;
  // Take copies of the nonbasic information for the dual LP, since
  // its values will be over-written in constructing the corresponding
  // data for the primal problem
  vector<int8_t> dual_nonbasic_flag = basis_.nonbasicFlag_;
  vector<HighsInt>& primal_basic_index = basis_.basicIndex_;
  vector<int8_t>& primal_nonbasic_flag = basis_.nonbasicFlag_;
  vector<int8_t>& primal_nonbasic_move = basis_.nonbasicMove_;
  basis_.nonbasicFlag_.assign(primal_num_tot, kIllegalFlagValue);
  basis_.nonbasicMove_.assign(primal_num_tot, kIllegalMoveValue);
  basis_.basicIndex_.resize(0);
  // The number of dual rows is the number of primal columns, so all
  // dual basic variables are nonbasic in the primal problem.
  //
  // If there are extra dual columns due to upper bounds on boxed
  // primal variables/constraints, there will be an excess of dual
  // nonbasic variables for the required primal basic variables.
  //
  // For each pair of dual variables associated with a boxed primal
  // variable/constraint:
  //
  // * If one is basic then it yields a nonbasic primal variable, at
  //   a bound given by the basic dual
  //
  // * If both are nonbasic, then they yield a basic primal variable
  //
  // Keep track of the dual column added to handle upper bounds on
  // boxed variables/constraints
  HighsInt upper_bound_col = original_num_row_;
  for (HighsInt iCol = 0; iCol < original_num_col_; iCol++) {
    const double lower = original_col_lower_[iCol];
    const double upper = original_col_upper_[iCol];
    int8_t move = kIllegalMoveValue;
    HighsInt dual_variable = dual_num_col + iCol;
    bool dual_basic = dual_nonbasic_flag[dual_variable] == kNonbasicFlagFalse;
    if (lower == upper) {
      // Fixed
      if (dual_basic) move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (dual_basic) {
          // Primal variable is nonbasic at its lower bound
          move = kNonbasicMoveUp;
        } else {
          // Look at the corresponding dual variable for the upper bound
          dual_variable = upper_bound_col;
          dual_basic = dual_nonbasic_flag[dual_variable] == kNonbasicFlagFalse;
          if (dual_basic) {
            // Primal variable is nonbasic at its upper bound
            move = kNonbasicMoveDn;
          }
        }
        upper_bound_col++;
      } else {
        // Lower (since upper bound is infinite)
        if (dual_basic) move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      if (dual_basic) move = kNonbasicMoveDn;
    } else {
      // FREE
      //
      // Dual activity a^Ty is fixed by cost, implying dual for primal
      // column (slack for dual row) is fixed at zero
      if (dual_basic) {
        assert(4 == 0);
        move = kNonbasicMoveZe;
      }
    }
    if (dual_basic) {
      // Primal nonbasic column from basic dual row
      assert(move != kIllegalMoveValue);
      primal_nonbasic_flag[iCol] = kNonbasicFlagTrue;
      primal_nonbasic_move[iCol] = move;
    } else {
      // Primal basic column from nonbasic dual row
      primal_basic_index.push_back(iCol);
      primal_nonbasic_flag[iCol] = kNonbasicFlagFalse;
      primal_nonbasic_move[iCol] = 0;
    }
  }
  for (HighsInt iRow = 0; iRow < original_num_row_; iRow++) {
    double lower = original_row_lower_[iRow];
    double upper = original_row_upper_[iRow];
    int8_t move = kIllegalMoveValue;
    HighsInt dual_variable = iRow;
    bool dual_basic = dual_nonbasic_flag[dual_variable] == kNonbasicFlagFalse;
    if (lower == upper) {
      // Equality constraint
      //
      // Dual variable has primal RHS as cost and is free
      if (dual_basic) {
        move = kNonbasicMoveZe;
      }
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (dual_basic) {
          // Primal variable is nonbasic at its lower bound
          move = kNonbasicMoveDn;
        } else {
          // Look at the corresponding dual variable for the upper bound
          dual_variable = upper_bound_col;
          dual_basic = dual_nonbasic_flag[dual_variable] == kNonbasicFlagFalse;
          if (dual_basic) {
            // Primal variable is nonbasic at its upper bound
            move = kNonbasicMoveUp;
          }
        }
        upper_bound_col++;
      } else {
        // Lower (since upper bound is infinite)
        if (dual_basic) {
          move = kNonbasicMoveDn;
        }
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      if (dual_basic) {
        move = kNonbasicMoveUp;
      }
    } else {
      // FREE
      if (dual_basic) {
        assert(14 == 0);
        move = kNonbasicMoveZe;
      }
      // Shouldn't get free rows, but handle them anyway
    }
    if (dual_basic) {
      // Primal nonbasic column from basic dual row
      assert(move != kIllegalMoveValue);
      primal_nonbasic_flag[original_num_col_ + iRow] = kNonbasicFlagTrue;
      primal_nonbasic_move[original_num_col_ + iRow] = move;
    } else {
      // Primal basic column from nonbasic dual row
      primal_basic_index.push_back(original_num_col_ + iRow);
      primal_nonbasic_flag[original_num_col_ + iRow] = kNonbasicFlagFalse;
      primal_nonbasic_move[original_num_col_ + iRow] = 0;
    }
  }
  const bool ignore_scaling = true;
  if (!ignore_scaling) {
    // Flip any scale factors
    if (lp_.scale_.has_scaling) {
      std::vector<double> temp_scale = lp_.scale_.row;
      lp_.scale_.row = lp_.scale_.col;
      lp_.scale_.col = temp_scale;
      lp_.scale_.col.resize(original_num_col_);
      lp_.scale_.row.resize(original_num_row_);
      lp_.scale_.num_col = original_num_col_;
      lp_.scale_.num_row = original_num_row_;
    }
  }
  // Change optimization sense
  if (lp_.sense_ == ObjSense::kMinimize) {
    lp_.sense_ = ObjSense::kMaximize;
  } else {
    lp_.sense_ = ObjSense::kMinimize;
  }
  // Flip LP dimensions
  lp_.num_col_ = original_num_col_;
  lp_.num_row_ = original_num_row_;
  // Restore the original offset
  lp_.offset_ = original_offset_;
  // Copy back the costs and bounds
  lp_.col_cost_ = original_col_cost_;
  lp_.col_lower_ = original_col_lower_;
  lp_.col_upper_ = original_col_upper_;
  lp_.row_lower_ = original_row_lower_;
  lp_.row_upper_ = original_row_upper_;
  // The primal constraint matrix is available row-wise as the first
  // original_num_row_ vectors of the dual constraint matrix
  HighsSparseMatrix primal_matrix;
  primal_matrix.start_.resize(original_num_row_ + 1);
  primal_matrix.index_.resize(original_num_nz_);
  primal_matrix.value_.resize(original_num_nz_);

  for (HighsInt iCol = 0; iCol < original_num_row_ + 1; iCol++)
    primal_matrix.start_[iCol] = lp_.a_matrix_.start_[iCol];
  for (HighsInt iEl = 0; iEl < original_num_nz_; iEl++) {
    primal_matrix.index_[iEl] = lp_.a_matrix_.index_[iEl];
    primal_matrix.value_[iEl] = lp_.a_matrix_.value_[iEl];
  }
  primal_matrix.num_col_ = original_num_col_;
  primal_matrix.num_row_ = original_num_row_;
  primal_matrix.format_ = MatrixFormat::kRowwise;
  // Copy the row-wise primal LP constraint matrix and transpose it.
  lp_.a_matrix_ = primal_matrix.getColwise();
  assert(lp_.a_matrix_.isColwise());
  // lp_.a_matrix_.ensureColwise(); is now redundant

  // Some sanity checks
  assert(lp_.num_col_ == original_num_col_);
  assert(lp_.num_row_ == original_num_row_);
  assert(lp_.numNz() == original_num_nz_);
  HighsInt num_basic_variables = primal_basic_index.size();
  bool num_basic_variables_ok = num_basic_variables == original_num_row_;
  if (!num_basic_variables_ok)
    printf("HEkk::undualize: Have %d basic variables, not %d\n",
           (int)num_basic_variables, (int)original_num_row_);
  assert(num_basic_variables_ok);

  // Clear the data retained when solving dual LP
  clearEkkDualize();
  status_.is_dualized = false;
  // Now solve with this basis. Should just be a case of reinverting
  // and re-solving for optimal primal and dual values, but
  // numerically marginal LPs will need clean-up
  status_.has_basis = true;
  status_.has_ar_matrix = false;
  status_.has_nla = false;
  status_.has_invert = false;
  HighsInt primal_solve_iteration_count = -iteration_count_;
  HighsStatus return_status = solve();
  primal_solve_iteration_count += iteration_count_;
  //  if (primal_solve_iteration_count)
  highsLogUser(options_->log_options, HighsLogType::kInfo,
               "Solving the primal LP (%s) using the optimal basis of its dual "
               "required %d simplex iterations\n",
               lp_.model_name_.c_str(), (int)primal_solve_iteration_count);
  return return_status;
}

HighsStatus HEkk::permute() {
  assert(1 == 0);
  return HighsStatus::kError;
}

HighsStatus HEkk::unpermute() {
  if (!this->status_.is_permuted) return HighsStatus::kOk;
  assert(1 == 0);
  return HighsStatus::kError;
}

HighsStatus HEkk::solve(const bool force_phase2) {
  debugInitialise();

  initialiseAnalysis();
  initialiseControl();

  if (analysis_.analyse_simplex_time)
    analysis_.simplexTimerStart(SimplexTotalClock);
  dual_simplex_cleanup_level_ = 0;
  dual_simplex_phase1_cleanup_level_ = 0;

  previous_iteration_cycling_detected = -kHighsIInf;

  initialiseForSolve();

  const HighsDebugStatus simplex_nla_status =
      simplex_nla_.debugCheckData("Before HEkk::solve()");
  const bool simplex_nla_ok = simplex_nla_status == HighsDebugStatus::kOk;
  if (!simplex_nla_ok) {
    highsLogUser(options_->log_options, HighsLogType::kError,
                 "Error in simplex NLA data\n");
    assert(simplex_nla_ok);
    return returnFromEkkSolve(HighsStatus::kError);
  }

  const bool report_initial_basis = false;
  if (report_initial_basis) debugReportInitialBasis();

  assert(status_.has_basis);
  assert(status_.has_invert);
  assert(status_.initialised_for_solve);
  if (model_status_ == HighsModelStatus::kOptimal)
    return returnFromEkkSolve(HighsStatus::kOk);

  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  std::string algorithm_name;

  // Indicate that dual and primal rays are not known
  this->clearRayRecords();

  // Allow primal and dual perturbations in case a block on them is
  // hanging over from a previous call
  info_.allow_cost_shifting = true;
  info_.allow_cost_perturbation = true;
  info_.allow_bound_perturbation = true;

  chooseSimplexStrategyThreads(*options_, info_);
  HighsInt& simplex_strategy = info_.simplex_strategy;

  // Initial solve according to strategy
  if (simplex_strategy == kSimplexStrategyPrimal) {
    algorithm_name = "primal";
    reportSimplexPhaseIterations(options_->log_options, iteration_count_, info_,
                                 true);
    highsLogUser(options_->log_options, HighsLogType::kInfo, "Using %s\n",
                 simplexStrategyToString(kSimplexStrategyPrimal).c_str());
    HEkkPrimal primal_solver(*this);
    call_status = primal_solver.solve(force_phase2);
    assert(called_return_from_solve_);
    return_status = interpretCallStatus(options_->log_options, call_status,
                                        return_status, "HEkkPrimal::solve");
  } else {
    algorithm_name = "dual";
    reportSimplexPhaseIterations(options_->log_options, iteration_count_, info_,
                                 true);
    // Solve, depending on the particular strategy
    if (simplex_strategy == kSimplexStrategyDualTasks) {
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   "Using %s with concurrency of %d\n",
                   simplexStrategyToString(kSimplexStrategyDualTasks).c_str(),
                   int(info_.num_concurrency));
    } else if (simplex_strategy == kSimplexStrategyDualMulti) {
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   "Using %s with concurrency of %d\n",
                   simplexStrategyToString(kSimplexStrategyDualMulti).c_str(),
                   int(info_.num_concurrency));
    } else {
      std::stringstream ss;
      ss.str(std::string());
      assert(simplex_strategy == kSimplexStrategyChoose ||
             simplex_strategy == kSimplexStrategyDual);
      ss << highsFormatToString(
          "Choosing to use %s%s",
          options_->parallel == kHighsOnString ? "serial " : "",
          simplexStrategyToString(kSimplexStrategyDual).c_str());
      if (options_->parallel == kHighsOnString) {
        ss << highsFormatToString(" despite setting parallel = on\n");
        ss << highsFormatToString(
            "To force parallel dual simplex, set simplex_strategy to "
            "kSimplexStrategyDualTasks = %d or (preferably) "
            "kSimplexStrategyDualMulti = %d\n",
            kSimplexStrategyDualTasks, kSimplexStrategyDualMulti);
      }
      highsLogUser(options_->log_options, HighsLogType::kInfo, "%s\n",
                   ss.str().c_str());

      highsLogUser(options_->log_options, HighsLogType::kInfo, "Using %s\n",
                   simplexStrategyToString(kSimplexStrategyDual).c_str());
    }
    HEkkDual dual_solver(*this);
    call_status = dual_solver.solve(force_phase2);
    assert(called_return_from_solve_);
    return_status = interpretCallStatus(options_->log_options, call_status,
                                        return_status, "HEkkDual::solve");

    // Dual simplex solver may set model_status to be
    // kUnboundedOrInfeasible, and Highs::run() may not allow that to
    // be returned, so use primal simplex to distinguish
    if (model_status_ == HighsModelStatus::kUnboundedOrInfeasible &&
        !options_->allow_unbounded_or_infeasible) {
      HEkkPrimal primal_solver(*this);
      call_status = primal_solver.solve();
      assert(called_return_from_solve_);
      return_status = interpretCallStatus(options_->log_options, call_status,
                                          return_status, "HEkkPrimal::solve");
    }
  }

  reportSimplexPhaseIterations(options_->log_options, iteration_count_, info_);
  if (return_status == HighsStatus::kError)
    return returnFromEkkSolve(return_status);
  highsLogDev(options_->log_options, HighsLogType::kInfo,
              "%s simplex solver returns %" HIGHSINT_FORMAT
              " primal and %" HIGHSINT_FORMAT
              " dual infeasibilities: "
              "Status %s\n",
              algorithm_name.c_str(), info_.num_primal_infeasibilities,
              info_.num_dual_infeasibilities,
              utilModelStatusToString(model_status_).c_str());
  // Can model_status_ = HighsModelStatus::kNotset be returned?
  assert(model_status_ != HighsModelStatus::kNotset);

  if (analysis_.analyse_simplex_summary_data) analysis_.summaryReport();
  if (analysis_.analyse_factor_data) analysis_.reportInvertFormData();
  if (analysis_.analyse_factor_time) analysis_.reportFactorTimer();
  return returnFromEkkSolve(return_status);
}

HighsStatus HEkk::setBasis() {
  // Set up nonbasicFlag and basicIndex for a logical basis
  const HighsInt num_col = lp_.num_col_;
  const HighsInt num_row = lp_.num_row_;

  basis_.setup(num_col, num_row);
  basis_.debug_origin_name = "HEkk::setBasis - logical";

  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    basis_.nonbasicFlag_[iCol] = kNonbasicFlagTrue;
    double lower = lp_.col_lower_[iCol];
    double upper = lp_.col_upper_[iCol];
    int8_t move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed. Set to bound of LP that is closer to
        // zero
        if (move == kIllegalMoveValue) {
          if (fabs(lower) < fabs(upper)) {
            move = kNonbasicMoveUp;
          } else {
            move = kNonbasicMoveDn;
          }
        }
      } else {
        // Lower (since upper bound is infinite)
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = kNonbasicMoveDn;
    } else {
      // FREE
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis_.nonbasicMove_[iCol] = move;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    basis_.nonbasicFlag_[iVar] = kNonbasicFlagFalse;
    HighsHashHelpers::sparse_combine(basis_.hash, iVar);
    basis_.basicIndex_[iRow] = iVar;
  }
  info_.num_basic_logicals = num_row;
  status_.has_basis = true;
  return HighsStatus::kOk;
}

HighsStatus HEkk::setBasis(const HighsBasis& highs_basis) {
  // Shouldn't have to check the incoming basis since this is an
  // internal call, but it may be a basis that's set up internally
  // with errors :-) ...
  //
  if (kDebugMipNodeDualFeasible) {
    // The basis should be dual feasible unless it was alien
    debug_dual_feasible = !highs_basis.was_alien;
  } else {
    assert(!debug_dual_feasible);
  }
  HighsOptions& options = *options_;
  if (debugHighsBasisConsistent(options, lp_, highs_basis) ==
      HighsDebugStatus::kLogicalError) {
    highsLogDev(options_->log_options, HighsLogType::kError,
                "Supposed to be a Highs basis, but not valid\n");
    return HighsStatus::kError;
  }
  HighsInt num_col = lp_.num_col_;
  HighsInt num_row = lp_.num_row_;
  // Set up the basis in case it has not yet been done for this LP
  basis_.setup(num_col, num_row);
  basis_.debug_id = highs_basis.debug_id;
  basis_.debug_update_count = highs_basis.debug_update_count;
  basis_.debug_origin_name = highs_basis.debug_origin_name;
  assert(basis_.debug_origin_name != "");
  HighsInt num_basic_variables = 0;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    HighsInt iVar = iCol;
    const double lower = lp_.col_lower_[iCol];
    const double upper = lp_.col_upper_[iCol];
    if (highs_basis.col_status[iCol] == HighsBasisStatus::kBasic) {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagFalse;
      basis_.nonbasicMove_[iVar] = 0;
      basis_.basicIndex_[num_basic_variables++] = iVar;
      HighsHashHelpers::sparse_combine(basis_.hash, iVar);
    } else {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagTrue;
      if (lower == upper) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      } else if (highs_basis.col_status[iCol] == HighsBasisStatus::kLower) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveUp;
      } else if (highs_basis.col_status[iCol] == HighsBasisStatus::kUpper) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveDn;
      } else {
        assert(highs_basis.col_status[iCol] == HighsBasisStatus::kZero);
        basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      }
    }
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    const double lower = lp_.row_lower_[iRow];
    const double upper = lp_.row_upper_[iRow];
    if (highs_basis.row_status[iRow] == HighsBasisStatus::kBasic) {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagFalse;
      basis_.nonbasicMove_[iVar] = 0;
      basis_.basicIndex_[num_basic_variables++] = iVar;
      HighsHashHelpers::sparse_combine(basis_.hash, iVar);
    } else {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagTrue;
      if (lower == upper) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      } else if (highs_basis.row_status[iRow] == HighsBasisStatus::kLower) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveDn;
      } else if (highs_basis.row_status[iRow] == HighsBasisStatus::kUpper) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveUp;
      } else {
        assert(highs_basis.row_status[iRow] == HighsBasisStatus::kZero);
        basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      }
    }
  }
  status_.has_basis = true;
  return HighsStatus::kOk;
}

void HEkk::addCols(const HighsLp& lp,
                   const HighsSparseMatrix& scaled_a_matrix) {
  // Should be extendSimplexLpRandomVectors
  //  if (valid_simplex_basis)
  //    appendBasicRowsToBasis(simplex_lp, simplex_basis, XnumNewRow);
  //  ekk_instance_.updateStatus(LpAction::kNewRows);
  //  if (valid_simplex_lp) {
  //    simplex_lp.num_row_ += XnumNewRow;
  //    ekk_instance_.initialiseSimplexLpRandomVectors();
  //  }
  //  if (valid_simplex_lp)
  //    assert(ekk_instance_.lp_.dimensionsOk("addCols - simplex"));
  if (this->status_.has_nla) this->simplex_nla_.addCols(&lp);
  this->updateStatus(LpAction::kNewCols);
}

void HEkk::addRows(const HighsLp& lp,
                   const HighsSparseMatrix& scaled_ar_matrix) {
  // Should be extendSimplexLpRandomVectors
  //  if (valid_simplex_basis)
  //    appendBasicRowsToBasis(simplex_lp, simplex_basis, XnumNewRow);
  //  ekk_instance_.updateStatus(LpAction::kNewRows);
  //  if (valid_simplex_lp) {
  //    simplex_lp.num_row_ += XnumNewRow;
  //    ekk_instance_.initialiseSimplexLpRandomVectors();
  //  }
  //  if (valid_simplex_lp)
  //    assert(ekk_instance_.lp_.dimensionsOk("addRows - simplex"));
  if (kExtendInvertWhenAddingRows && this->status_.has_nla) {
    this->simplex_nla_.addRows(&lp, basis_.basicIndex_.data(),
                               &scaled_ar_matrix);
    setNlaPointersForTrans(lp);
    this->debugNlaCheckInvert("HEkk::addRows - on entry",
                              kHighsDebugLevelExpensive + 1);
  }
  // Update the number of rows in the simplex LP so that it's
  // consistent with simplex basis information
  this->lp_.num_row_ = lp.num_row_;
  this->updateStatus(LpAction::kNewRows);
}

void HEkk::deleteCols(const HighsIndexCollection& index_collection) {
  this->updateStatus(LpAction::kDelCols);
}
void HEkk::deleteRows(const HighsIndexCollection& index_collection) {
  this->updateStatus(LpAction::kDelRows);
}

void HEkk::unscaleSimplex(const HighsLp& incumbent_lp) {
  if (!this->simplex_in_scaled_space_) return;
  assert(incumbent_lp.scale_.has_scaling);
  const HighsInt num_col = incumbent_lp.num_col_;
  const HighsInt num_row = incumbent_lp.num_row_;
  const vector<double>& col_scale = incumbent_lp.scale_.col;
  const vector<double>& row_scale = incumbent_lp.scale_.row;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    const HighsInt iVar = iCol;
    const double factor = col_scale[iCol];
    this->info_.workCost_[iVar] /= factor;
    this->info_.workDual_[iVar] /= factor;
    this->info_.workShift_[iVar] /= factor;
    this->info_.workLower_[iVar] *= factor;
    this->info_.workUpper_[iVar] *= factor;
    this->info_.workRange_[iVar] *= factor;
    this->info_.workValue_[iVar] *= factor;
    this->info_.workLowerShift_[iVar] *= factor;
    this->info_.workUpperShift_[iVar] *= factor;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    const HighsInt iVar = num_col + iRow;
    const double factor = row_scale[iRow];
    this->info_.workCost_[iVar] *= factor;
    this->info_.workDual_[iVar] *= factor;
    this->info_.workShift_[iVar] *= factor;
    this->info_.workLower_[iVar] /= factor;
    this->info_.workUpper_[iVar] /= factor;
    this->info_.workRange_[iVar] /= factor;
    this->info_.workValue_[iVar] /= factor;
    this->info_.workLowerShift_[iVar] /= factor;
    this->info_.workUpperShift_[iVar] /= factor;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    double factor;
    const HighsInt iVar = this->basis_.basicIndex_[iRow];
    if (iVar < num_col) {
      factor = col_scale[iVar];
    } else {
      factor = 1.0 / row_scale[iVar - num_col];
    }
    this->info_.baseLower_[iRow] *= factor;
    this->info_.baseUpper_[iRow] *= factor;
    this->info_.baseValue_[iRow] *= factor;
  }
  this->simplex_in_scaled_space_ = false;
}

HighsSolution HEkk::getSolution() {
  HighsSolution solution;
  // Scatter the basic primal values
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++)
    info_.workValue_[basis_.basicIndex_[iRow]] = info_.baseValue_[iRow];
  // Zero the basic dual values
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++)
    info_.workDual_[basis_.basicIndex_[iRow]] = 0;

  // Now we can get the solution
  solution.col_value.resize(lp_.num_col_);
  solution.col_dual.resize(lp_.num_col_);
  solution.row_value.resize(lp_.num_row_);
  solution.row_dual.resize(lp_.num_row_);

  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    solution.col_value[iCol] = info_.workValue_[iCol];
    solution.col_dual[iCol] = (HighsInt)lp_.sense_ * info_.workDual_[iCol];
  }
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    solution.row_value[iRow] = -info_.workValue_[lp_.num_col_ + iRow];
    // @FlipRowDual negate RHS
    solution.row_dual[iRow] =
        -(HighsInt)lp_.sense_ * info_.workDual_[lp_.num_col_ + iRow];
  }
  solution.value_valid = true;
  solution.dual_valid = true;
  return solution;
}

HighsBasis HEkk::getHighsBasis(HighsLp& use_lp) const {
  HighsInt num_col = use_lp.num_col_;
  HighsInt num_row = use_lp.num_row_;
  HighsBasis highs_basis;
  highs_basis.col_status.resize(num_col);
  highs_basis.row_status.resize(num_row);
  assert(status_.has_basis);
  highs_basis.valid = false;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    HighsInt iVar = iCol;
    const double lower = use_lp.col_lower_[iCol];
    const double upper = use_lp.col_upper_[iCol];
    HighsBasisStatus basis_status = HighsBasisStatus::kNonbasic;
    if (!basis_.nonbasicFlag_[iVar]) {
      basis_status = HighsBasisStatus::kBasic;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveUp) {
      basis_status = HighsBasisStatus::kLower;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveDn) {
      basis_status = HighsBasisStatus::kUpper;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveZe) {
      if (lower == upper) {
        const double dual = (HighsInt)lp_.sense_ * info_.workDual_[iCol];
        basis_status =
            dual >= 0 ? HighsBasisStatus::kLower : HighsBasisStatus::kUpper;
      } else {
        basis_status = HighsBasisStatus::kZero;
      }
    }
    highs_basis.col_status[iCol] = basis_status;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    const double lower = use_lp.row_lower_[iRow];
    const double upper = use_lp.row_upper_[iRow];
    HighsBasisStatus basis_status = HighsBasisStatus::kNonbasic;
    if (!basis_.nonbasicFlag_[iVar]) {
      basis_status = HighsBasisStatus::kBasic;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveUp) {
      basis_status = HighsBasisStatus::kUpper;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveDn) {
      basis_status = HighsBasisStatus::kLower;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveZe) {
      if (lower == upper) {
        const double dual = (HighsInt)lp_.sense_ * info_.workDual_[iVar];
        basis_status =
            dual >= 0 ? HighsBasisStatus::kLower : HighsBasisStatus::kUpper;
      } else {
        basis_status = HighsBasisStatus::kZero;
      }
    }
    highs_basis.row_status[iRow] = basis_status;
  }
  highs_basis.valid = true;
  highs_basis.alien = false;
  highs_basis.useful = true;
  highs_basis.was_alien = false;
  highs_basis.debug_id =
      (HighsInt)(build_synthetic_tick_ + total_synthetic_tick_);
  highs_basis.debug_update_count = info_.update_count;
  highs_basis.debug_origin_name = basis_.debug_origin_name;
  return highs_basis;
}

HighsStatus HEkk::initialiseSimplexLpBasisAndFactor(
    const bool only_from_known_basis) {
  // This is normally called only from HEkk::initialiseForSolve, with
  // only_from_known_basis = false.
  //
  // In this case, if there is no existing simplex basis then a
  // logical basis is set up. Otherwise the existing simplex basis is
  // factorized, with logicals introduced to handle rank deficiency.
  //
  // It is also called from HighsSolution's
  // formSimplexLpBasisAndFactor, with only_from_known_basis = true or
  // false. In both cases a simplex basis is known.
  //
  // Calls with only_from_known_basis = true originate from
  // Highs::getBasicVariablesInterface, and are made when
  // Highs::getBasicVariables has been called and there is no
  // factorization of the current basis matrix. No rank deficiency
  // handling is permitted so, if it occurs, an error must be
  // returned.
  //
  // Calls with only_from_known_basis = false originate from
  // Highs::setBasis. The simplex basis should be non-singular since
  // it's come from a HighsBasis that is either non-alien, or from an
  // alien HighsBasis that's been checked/completed. However, it's
  // conceivable that a singularity could occur, and it's fine to
  // accommodate it.
  //
  // If only_from_known_basis is true, then there should be a simplex
  // basis to use
  if (only_from_known_basis) assert(status_.has_basis);
  // If there is no simplex basis, set up a logical basis
  if (!status_.has_basis) setBasis();
  // The simplex NLA operates in the scaled space if the LP has
  // scaling factors. If they exist but haven't been applied, then the
  // simplex NLA needs a separate, scaled constraint matrix. Thus
  // getScaledAMatrixPointer() returns a pointer to either the
  // constraint matrix or a scaled copy (that is a member of the HEkk
  // class), with the latter returned if the LP has scaling factors
  // but is unscaled.
  //
  HighsSparseMatrix* local_scaled_a_matrix = getScaledAMatrixPointer();
  //
  // If simplex NLA is set up, pass the pointers that it uses. It
  // deduces any scaling factors that it must use by inspecting
  // whether the LP has scaling factors, and whether it is scaled.
  //
  // If simplex NLA is not set up, then it will be done if
  //
  if (this->status_.has_nla) {
    assert(lpFactorRowCompatible());
    this->simplex_nla_.setPointers(
        &(this->lp_), local_scaled_a_matrix, this->basis_.basicIndex_.data(),
        this->options_, this->timer_, &(this->analysis_));
  } else {
    // todo @ Julian: this fails on glass4
    assert(info_.factor_pivot_threshold >= options_->factor_pivot_threshold);
    simplex_nla_.setup(
        &(this->lp_),                     //&lp_,
        this->basis_.basicIndex_.data(),  // basis_.basicIndex_.data(),
        this->options_,                   // options_,
        this->timer_,                     // timer_,
        &(this->analysis_),               //&analysis_,
        local_scaled_a_matrix, this->info_.factor_pivot_threshold);
    status_.has_nla = true;
  }

  if (!status_.has_invert) {
    const HighsInt rank_deficiency = computeFactor();
    if (rank_deficiency) {
      // Basis is rank deficient
      highsLogDev(
          options_->log_options, HighsLogType::kInfo,
          "HEkk::initialiseSimplexLpBasisAndFactor (%s) Rank_deficiency %d: Id "
          "= "
          "%d; UpdateCount = %d\n",
          basis_.debug_origin_name.c_str(), (int)rank_deficiency,
          (int)basis_.debug_id, (int)basis_.debug_update_count);
      if (only_from_known_basis) {
        // If only this basis should be used, then return error
        highsLogDev(options_->log_options, HighsLogType::kError,
                    "Supposed to be a full-rank basis, but incorrect\n");
        return HighsStatus::kError;
      }
      // Account for rank deficiency by correcting nonbasicFlag
      handleRankDeficiency();
      this->updateStatus(LpAction::kNewBasis);
      setNonbasicMove();
      status_.has_basis = true;
      status_.has_invert = true;
      status_.has_fresh_invert = true;
    }
    // Record the synthetic clock for INVERT, and zero it for UPDATE
    resetSyntheticClock();
  }
  assert(status_.has_invert);
  return HighsStatus::kOk;
}

void HEkk::handleRankDeficiency() {
  HFactor& factor = simplex_nla_.factor_;
  HighsInt rank_deficiency = factor.rank_deficiency;
  vector<HighsInt>& row_with_no_pivot = factor.row_with_no_pivot;
  vector<HighsInt>& var_with_no_pivot = factor.var_with_no_pivot;
  for (HighsInt k = 0; k < rank_deficiency; k++) {
    HighsInt row_in = row_with_no_pivot[k];
    HighsInt variable_in = lp_.num_col_ + row_in;
    HighsInt variable_out = var_with_no_pivot[k];
    basis_.nonbasicFlag_[variable_in] = kNonbasicFlagFalse;
    basis_.nonbasicFlag_[variable_out] = kNonbasicFlagTrue;
    HighsInt row_out = row_with_no_pivot[k];
    assert(basis_.basicIndex_[row_out] == variable_in);
    highsLogDev(
        options_->log_options, HighsLogType::kInfo,
        "HEkk::handleRankDeficiency: %4d: Basic row of leaving variable (%4d "
        "is %s %4d) is "
        "%4d; Entering logical = %4d is variable %d)\n",
        (int)k, (int)variable_out,
        variable_out < lp_.num_col_ ? " column" : "logical",
        variable_out < lp_.num_col_ ? (int)variable_out
                                    : (int)(variable_out - lp_.num_col_),
        (int)row_out, (int)(row_in), (int)variable_in);
    // NB Parameters row_out, variable_in, variable_out, since
    // variable_in is the logical that must not come out to be
    // replaced by the structural variable_out
    addBadBasisChange(row_out, variable_in, variable_out,
                      BadBasisChangeReason::kSingular, true);
  }
  status_.has_ar_matrix = false;
}

// Private methods

void HEkk::initialiseEkk() {
  if (status_.initialised_for_new_lp) return;
  setSimplexOptions();
  initialiseControl();
  initialiseSimplexLpRandomVectors();
  simplex_nla_.clear();
  clearBadBasisChange();
  status_.initialised_for_new_lp = true;
}

bool HEkk::isUnconstrainedLp() const {
  bool is_unconstrained_lp = lp_.num_row_ <= 0;
  if (is_unconstrained_lp)
    highsLogDev(
        options_->log_options, HighsLogType::kError,
        "HEkkDual::solve called for LP with non-positive (%" HIGHSINT_FORMAT
        ") number of constraints\n",
        lp_.num_row_);
  assert(!is_unconstrained_lp);
  return is_unconstrained_lp;
}

void HEkk::initialiseForSolve() {
  const HighsStatus return_status = initialiseSimplexLpBasisAndFactor();
  assert(return_status == HighsStatus::kOk);
  assert(status_.has_basis);

  updateSimplexOptions();
  initialiseSimplexLpRandomVectors();
  initialisePartitionedRowwiseMatrix();  // Timed
  allocateWorkAndBaseArrays();
  initialiseCost(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown, false);
  initialiseBound(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown, false);
  initialiseNonbasicValueAndMove();
  computePrimal();                // Timed
  computeDual();                  // Timed
  computeSimplexInfeasible();     // Timed
  computeDualObjectiveValue();    // Timed
  computePrimalObjectiveValue();  // Timed
  status_.initialised_for_solve = true;

  bool primal_feasible = info_.num_primal_infeasibilities == 0;
  bool dual_feasible = info_.num_dual_infeasibilities == 0;
  visited_basis_.clear();
  visited_basis_.insert(basis_.hash);
  model_status_ = HighsModelStatus::kNotset;
  if (primal_feasible && dual_feasible)
    model_status_ = HighsModelStatus::kOptimal;
}

void HEkk::setSimplexOptions() {
  // Copy values of HighsOptions for the simplex solver
  // Currently most of these options are straight copies, but they
  // will become valuable when "choose" becomes a HiGHS strategy value
  // that will need converting into a specific simplex strategy value.
  //
  // NB simplex_strategy is set by chooseSimplexStrategyThreads in each call
  //
  info_.dual_edge_weight_strategy = options_->simplex_dual_edge_weight_strategy;
  info_.price_strategy = options_->simplex_price_strategy;
  info_.dual_simplex_cost_perturbation_multiplier =
      options_->dual_simplex_cost_perturbation_multiplier;
  info_.primal_simplex_bound_perturbation_multiplier =
      options_->primal_simplex_bound_perturbation_multiplier;
  info_.factor_pivot_threshold = options_->factor_pivot_threshold;
  info_.update_limit = options_->simplex_update_limit;
  random_.initialise(options_->random_seed);

  // Set values of internal options
  info_.store_squared_primal_infeasibility = true;
}

void HEkk::updateSimplexOptions() {
  // Update some simplex option values from HighsOptions when
  // (re-)solving an LP. Others aren't changed because better values
  // may have been learned due to solving this LP (possibly with some
  // modification) before.
  //
  // NB simplex_strategy is set by chooseSimplexStrategyThreads in each call
  //
  info_.dual_simplex_cost_perturbation_multiplier =
      options_->dual_simplex_cost_perturbation_multiplier;
  info_.primal_simplex_bound_perturbation_multiplier =
      options_->primal_simplex_bound_perturbation_multiplier;
}

void HEkk::initialiseSimplexLpRandomVectors() {
  const HighsInt num_col = lp_.num_col_;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  if (!num_tot) return;
  // Instantiate and (re-)initialise the random number generator
  //  HighsRandom random;
  HighsRandom& random = random_;
  //  random.initialise();

  if (num_col) {
    // Generate a random permutation of the column indices
    vector<HighsInt>& numColPermutation = info_.numColPermutation_;
    numColPermutation.resize(num_col);
    for (HighsInt i = 0; i < num_col; i++) numColPermutation[i] = i;
    random.shuffle(numColPermutation.data(), num_col);
  }

  // Re-initialise the random number generator and generate the
  // random vectors in the same order as hsol to maintain repeatable
  // performance
  // random.initialise();

  // Generate a random permutation of all the indices
  vector<HighsInt>& numTotPermutation = info_.numTotPermutation_;
  numTotPermutation.resize(num_tot);
  for (HighsInt i = 0; i < num_tot; i++) numTotPermutation[i] = i;
  random.shuffle(numTotPermutation.data(), num_tot);

  // Generate a vector of random reals
  info_.numTotRandomValue_.resize(num_tot);
  vector<double>& numTotRandomValue = info_.numTotRandomValue_;
  for (HighsInt i = 0; i < num_tot; i++) {
    numTotRandomValue[i] = random.fraction();
  }
}

void HEkk::chooseSimplexStrategyThreads(const HighsOptions& options,
                                        HighsSimplexInfo& info) {
  // Ensure that this is not called with an optimal basis
  assert(info.num_dual_infeasibilities > 0 ||
         info.num_primal_infeasibilities > 0);
  // Set the internal simplex strategy and number of threads for dual
  // simplex
  HighsInt& simplex_strategy = info.simplex_strategy;
  // By default, use the HighsOptions strategy. If this is
  // kSimplexStrategyChoose, then the strategy used will depend on
  // whether the current basis is primal feasible.
  simplex_strategy = options.simplex_strategy;
  if (simplex_strategy == kSimplexStrategyChoose) {
    // HiGHS is left to choose the simplex strategy
    if (info.num_primal_infeasibilities > 0) {
      // Not primal feasible, so use dual simplex
      simplex_strategy = kSimplexStrategyDual;
    } else {
      // Primal feasible. so use primal simplex
      simplex_strategy = kSimplexStrategyPrimal;
    }
  }
  // Set min/max_threads to correspond to serial code. They will be
  // set to other values if parallel options are used.
  info.min_concurrency = 1;
  info.max_concurrency = 1;
  // Record the min/max minimum concurrency in the options
  const HighsInt simplex_min_concurrency = options.simplex_min_concurrency;
  const HighsInt simplex_max_concurrency = options.simplex_max_concurrency;
  HighsInt max_threads = highs::parallel::num_threads();

  // Don't allow the parallel option to switch from
  // kSimplexStrategyDual to kSimplexStrategyDualMulti so the MIP
  // solver only ever uses serial simplex
  //
  //  if (options.parallel == kHighsOnString &&
  //      simplex_strategy == kSimplexStrategyDual) {
  //    // The parallel strategy is on and the simplex strategy is dual so use
  //    // PAMI if there are enough threads
  //    if (max_threads >= kDualMultiMinConcurrency)
  //      simplex_strategy = kSimplexStrategyDualMulti;
  //  }

  //
  // If parallel strategies are used, the minimum concurrency will be
  // set to be at least the minimum required for the strategy
  //
  // All this is independent of the number of threads available, since
  // code with multiple concurrency can be run in serial.

  if (simplex_strategy == kSimplexStrategyDualTasks) {
    info.min_concurrency =
        max(kDualTasksMinConcurrency, simplex_min_concurrency);
    info.max_concurrency = max(info.min_concurrency, simplex_max_concurrency);
  } else if (simplex_strategy == kSimplexStrategyDualMulti) {
    info.min_concurrency =
        max(kDualMultiMinConcurrency, simplex_min_concurrency);
    info.max_concurrency = max(info.min_concurrency, simplex_max_concurrency);
  }

  // Set the concurrency to be used to be the maximum number
  info.num_concurrency = info.max_concurrency;
  // Give a warning if the concurrency to be used is less than the
  // minimum concurrency allowed
  if (info.num_concurrency < simplex_min_concurrency) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Using concurrency of %" HIGHSINT_FORMAT
                 " for parallel strategy rather than "
                 "minimum number (%" HIGHSINT_FORMAT ") specified in options\n",
                 info.num_concurrency, simplex_min_concurrency);
  }
  // Give a warning if the concurrency to be used is more than the
  // maximum concurrency allowed
  if (info.num_concurrency > simplex_max_concurrency) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Using concurrency of %" HIGHSINT_FORMAT
                 " for parallel strategy rather than "
                 "maximum number (%" HIGHSINT_FORMAT ") specified in options\n",
                 info.num_concurrency, simplex_max_concurrency);
  }
  // Give a warning if the concurrency to be used is less than the
  // number of threads available
  if (info.num_concurrency > max_threads) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Number of threads available = %" HIGHSINT_FORMAT
                 " < %" HIGHSINT_FORMAT
                 " = Simplex concurrency to be used: Parallel performance may "
                 "be less than anticipated\n",
                 max_threads, info.num_concurrency);
  }
}

bool HEkk::getNonsingularInverse(const HighsInt solve_phase) {
  assert(status_.has_basis);
  const vector<HighsInt>& basicIndex = basis_.basicIndex_;
  // Take a copy of basicIndex from before INVERT to be used as the
  // saved ordering of basic variables - so reinvert will run
  // identically.
  const vector<HighsInt> basicIndex_before_compute_factor = basicIndex;
  // Save the number of updates performed in case it has to be used to determine
  // a limit
  const HighsInt simplex_update_count = info_.update_count;
  // Dual simplex edge weights are identified with rows, so must be
  // permuted according to INVERT. Scatter the edge weights so that,
  // after INVERT, they can be gathered according to the new
  // permutation of basicIndex
  analysis_.simplexTimerStart(PermWtClock);
  for (HighsInt i = 0; i < lp_.num_row_; i++)
    scattered_dual_edge_weight_[basicIndex[i]] = dual_edge_weight_[i];
  analysis_.simplexTimerStop(PermWtClock);

  // Call computeFactor to perform INVERT
  HighsInt rank_deficiency = computeFactor();
  if (rank_deficiency)
    highsLogDev(
        options_->log_options, HighsLogType::kInfo,
        "HEkk::getNonsingularInverse Rank_deficiency: solve %d (Iteration "
        "%d)\n",
        (int)debug_solve_call_num_, (int)iteration_count_);

  const bool artificial_rank_deficiency = false;  //  true;//
  if (artificial_rank_deficiency) {
    if (!info_.phase1_backtracking_test_done && solve_phase == kSolvePhase1) {
      // Claim rank deficiency to test backtracking
      // printf("Phase1 (Iter %" HIGHSINT_FORMAT
      //        ") Claiming rank deficiency to test backtracking\n",
      //        iteration_count_);
      rank_deficiency = 1;
      info_.phase1_backtracking_test_done = true;
    } else if (!info_.phase2_backtracking_test_done &&
               solve_phase == kSolvePhase2) {
      // Claim rank deficiency to test backtracking
      // printf("Phase2 (Iter %" HIGHSINT_FORMAT
      //        ") Claiming rank deficiency to test backtracking\n",
      //        iteration_count_);
      rank_deficiency = 1;
      info_.phase2_backtracking_test_done = true;
    }
  }
  if (rank_deficiency) {
    // Rank deficient basis, so backtrack to last full rank basis
    //
    // Get the last nonsingular basis - so long as there is one
    uint64_t deficient_hash = basis_.hash;
    if (!getBacktrackingBasis()) return false;
    // Record that backtracking is taking place
    info_.backtracking_ = true;
    visited_basis_.clear();
    visited_basis_.insert(basis_.hash);
    visited_basis_.insert(deficient_hash);
    this->updateStatus(LpAction::kBacktracking);
    HighsInt backtrack_rank_deficiency = computeFactor();
    // This basis has previously been inverted successfully, so it shouldn't be
    // singular
    if (backtrack_rank_deficiency) return false;
    // simplex update limit will be half of the number of updates
    // performed, so make sure that at least one update was performed
    if (simplex_update_count <= 1) return false;
    HighsInt use_simplex_update_limit = info_.update_limit;
    HighsInt new_simplex_update_limit = simplex_update_count / 2;
    info_.update_limit = new_simplex_update_limit;
    highsLogDev(options_->log_options, HighsLogType::kWarning,
                "Rank deficiency of %" HIGHSINT_FORMAT
                " after %" HIGHSINT_FORMAT
                " simplex updates, so "
                "backtracking: max updates reduced from %" HIGHSINT_FORMAT
                " to %" HIGHSINT_FORMAT "\n",
                rank_deficiency, simplex_update_count, use_simplex_update_limit,
                new_simplex_update_limit);
  } else {
    // Current basis is full rank so save it
    putBacktrackingBasis(basicIndex_before_compute_factor);
    // Indicate that backtracking is not taking place
    info_.backtracking_ = false;
    // Reset the update limit in case this is the first successful
    // inversion after backtracking
    info_.update_limit = options_->simplex_update_limit;
  }
  // Gather the edge weights according to the permutation of
  // basicIndex after INVERT
  analysis_.simplexTimerStart(PermWtClock);
  for (HighsInt i = 0; i < lp_.num_row_; i++)
    dual_edge_weight_[i] = scattered_dual_edge_weight_[basicIndex[i]];
  analysis_.simplexTimerStop(PermWtClock);
  return true;
}

bool HEkk::getBacktrackingBasis() {
  if (!info_.valid_backtracking_basis_) return false;
  basis_ = info_.backtracking_basis_;
  info_.costs_shifted = (info_.backtracking_basis_costs_shifted_ != 0);
  info_.costs_perturbed = (info_.backtracking_basis_costs_perturbed_ != 0);
  info_.bounds_shifted = (info_.backtracking_basis_bounds_shifted_ != 0);
  info_.bounds_perturbed = (info_.backtracking_basis_bounds_perturbed_ != 0);
  info_.workShift_ = info_.backtracking_basis_workShift_;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++)
    scattered_dual_edge_weight_[iVar] =
        info_.backtracking_basis_edge_weight_[iVar];
  return true;
}

void HEkk::putBacktrackingBasis() {
  const vector<HighsInt>& basicIndex = basis_.basicIndex_;
  analysis_.simplexTimerStart(PermWtClock);
  for (HighsInt i = 0; i < lp_.num_row_; i++)
    scattered_dual_edge_weight_[basicIndex[i]] = dual_edge_weight_[i];
  analysis_.simplexTimerStop(PermWtClock);
  putBacktrackingBasis(basicIndex);
}

void HEkk::putBacktrackingBasis(
    const vector<HighsInt>& basicIndex_before_compute_factor) {
  info_.valid_backtracking_basis_ = true;
  info_.backtracking_basis_ = basis_;
  info_.backtracking_basis_.basicIndex_ = basicIndex_before_compute_factor;
  info_.backtracking_basis_costs_shifted_ = info_.costs_shifted;
  info_.backtracking_basis_costs_perturbed_ = info_.costs_perturbed;
  info_.backtracking_basis_bounds_shifted_ = info_.bounds_shifted;
  info_.backtracking_basis_bounds_perturbed_ = info_.bounds_perturbed;
  info_.backtracking_basis_workShift_ = info_.workShift_;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++)
    info_.backtracking_basis_edge_weight_[iVar] =
        scattered_dual_edge_weight_[iVar];
}

void HEkk::computePrimalObjectiveValue() {
  analysis_.simplexTimerStart(ComputePrObjClock);
  info_.primal_objective_value = 0;
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    HighsInt iVar = basis_.basicIndex_[iRow];
    if (iVar < lp_.num_col_) {
      info_.primal_objective_value +=
          info_.baseValue_[iRow] * lp_.col_cost_[iVar];
    }
  }
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    if (basis_.nonbasicFlag_[iCol])
      info_.primal_objective_value +=
          info_.workValue_[iCol] * lp_.col_cost_[iCol];
  }
  info_.primal_objective_value *= cost_scale_;
  // Objective value calculation is done using primal values and
  // original costs so offset is vanilla
  info_.primal_objective_value += lp_.offset_;
  // Now have primal objective value
  status_.has_primal_objective_value = true;
  analysis_.simplexTimerStop(ComputePrObjClock);
}

void HEkk::computeDualObjectiveValue(const HighsInt phase) {
  analysis_.simplexTimerStart(ComputeDuObjClock);
  info_.dual_objective_value = 0;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    if (basis_.nonbasicFlag_[iCol])
      info_.dual_objective_value +=
          info_.workValue_[iCol] * info_.workDual_[iCol];
  }
  info_.dual_objective_value *= cost_scale_;
  if (phase != 1) {
    // In phase 1 the dual objective has no objective
    // shift. Otherwise, if minimizing the shift is added. If
    // maximizing, workCost (and hence workDual) are negated, so the
    // shift is subtracted. Hence the shift is added according to the
    // sign implied by sense_
    info_.dual_objective_value += ((HighsInt)lp_.sense_) * lp_.offset_;
  }
  // Now have dual objective value
  status_.has_dual_objective_value = true;
  analysis_.simplexTimerStop(ComputeDuObjClock);
}

bool HEkk::rebuildRefactor(HighsInt rebuild_reason) {
  // If no updates have been performed, then don't refactor!
  if (info_.update_count == 0) return false;
  // Otherwise, refactor by default
  bool refactor = true;
  double solution_error = 0;
  if (options_->no_unnecessary_rebuild_refactor) {
    // Consider whether not to refactor in rebuild
    //
    // Must have an INVERT just to consider this!
    assert(status_.has_invert);
    if (rebuild_reason == kRebuildReasonNo ||
        rebuild_reason == kRebuildReasonPossiblyOptimal ||
        rebuild_reason == kRebuildReasonPossiblyPhase1Feasible ||
        rebuild_reason == kRebuildReasonPossiblyPrimalUnbounded ||
        rebuild_reason == kRebuildReasonPossiblyDualUnbounded ||
        rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex) {
      // By default, don't refactor!
      refactor = false;
      // Possibly revise the decision based on accuracy when solving a
      // test system
      const double error_tolerance =
          options_->rebuild_refactor_solution_error_tolerance;
      if (error_tolerance > 0) {
        solution_error = factorSolveError();
        refactor = solution_error > error_tolerance;
      }
    }
  }
  const bool report_refactorization = false;
  if (report_refactorization) {
    const std::string logic = refactor ? "   " : "no   ";
    if (info_.update_count &&
        rebuild_reason != kRebuildReasonSyntheticClockSaysInvert)
      printf(
          "%srefactorization after %4d updates and solution error = %11.4g for "
          "rebuild reason = %s\n",
          logic.c_str(), (int)info_.update_count, solution_error,
          rebuildReason(rebuild_reason).c_str());
  }
  return refactor;
}

HighsInt HEkk::computeFactor() {
  assert(status_.has_nla);
  if (status_.has_fresh_invert) return 0;
  // Clear any bad basis changes
  clearBadBasisChange();
  highsAssert(lpFactorRowCompatible(),
              "HEkk::computeFactor: lpFactorRowCompatible");
  // Perform INVERT
  analysis_.simplexTimerStart(InvertClock);
  const HighsInt rank_deficiency = simplex_nla_.invert();
  analysis_.simplexTimerStop(InvertClock);
  //
  // Set up hot start information
  hot_start_.refactor_info = simplex_nla_.factor_.refactor_info_;
  hot_start_.nonbasicMove = basis_.nonbasicMove_;
  hot_start_.valid = true;

  if (analysis_.analyse_factor_data)
    analysis_.updateInvertFormData(simplex_nla_.factor_);

  HighsInt alt_debug_level = -1;
  if (rank_deficiency) alt_debug_level = kHighsDebugLevelCostly;
  debugNlaCheckInvert("HEkk::computeFactor - original", alt_debug_level);

  if (rank_deficiency) {
    // Have an invertible representation, but of B with column(s)
    // replacements due to singularity. So no (fresh) representation of
    // B^{-1}
    status_.has_invert = false;
    status_.has_fresh_invert = false;
  } else {
    // Now have a representation of B^{-1}, and it is fresh!
    status_.has_invert = true;
    status_.has_fresh_invert = true;
  }
  // Set the update count to zero since the corrected invertible
  // representation may be used for an initial basis. In any case the
  // number of updates shouldn't be positive
  info_.update_count = 0;

  simplex_stats_.num_invert++;
  return rank_deficiency;
}

void HEkk::computeDualSteepestEdgeWeights(const bool initial) {
  if (analysis_.analyse_simplex_time) {
    analysis_.simplexTimerStart(SimplexIzDseWtClock);
    analysis_.simplexTimerStart(DseIzClock);
  }
  const HighsInt num_row = lp_.num_row_;
  HVector row_ep;
  row_ep.setup(num_row);
  assert((HighsInt)dual_edge_weight_.size() >= num_row);
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    dual_edge_weight_[iRow] = computeDualSteepestEdgeWeight(iRow, row_ep);
  if (analysis_.analyse_simplex_time) {
    analysis_.simplexTimerStop(SimplexIzDseWtClock);
    analysis_.simplexTimerStop(DseIzClock);
    if (initial) {
      double IzDseWtTT = analysis_.simplexTimerRead(SimplexIzDseWtClock);
      highsLogDev(options_->log_options, HighsLogType::kDetailed,
                  "Computed %" HIGHSINT_FORMAT " initial DSE weights in %gs\n",
                  num_row, IzDseWtTT);
    }
  }
}

double HEkk::computeDualSteepestEdgeWeight(const HighsInt iRow,
                                           HVector& row_ep) {
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = iRow;
  row_ep.array[iRow] = 1;
  row_ep.packFlag = false;
  simplex_nla_.btranInScaledSpace(row_ep, info_.row_ep_density,
                                  analysis_.pointer_serial_factor_clocks);
  const double local_row_ep_density = (1.0 * row_ep.count) / lp_.num_row_;
  updateOperationResultDensity(local_row_ep_density, info_.row_ep_density);
  return row_ep.norm2();
}

// Update the DSE weights
void HEkk::updateDualSteepestEdgeWeights(
    const HighsInt row_out, const HighsInt variable_in, const HVector* column,
    const double new_pivotal_edge_weight, const double Kai,
    const double* dual_steepest_edge_array) {
  analysis_.simplexTimerStart(DseUpdateWeightClock);

  const HighsInt num_row = lp_.num_row_;
  const HighsInt column_count = column->count;
  const HighsInt* variable_index = column->index.data();
  const double* column_array = column->array.data();

  const double col_aq_scale = simplex_nla_.variableScaleFactor(variable_in);
  const double col_ap_scale = simplex_nla_.basicColScaleFactor(row_out);
  const double inv_col_ap_scale = 1.0 / col_ap_scale;

  const bool DSE_check = false;
  HVector alt_dual_steepest_edge_column;
  HVector alt_pivotal_column;
  if (DSE_check) {
    // Compute the DSE column otherwise to check
    alt_dual_steepest_edge_column.setup(num_row);
    alt_dual_steepest_edge_column.clear();
    alt_dual_steepest_edge_column.count = 1;
    alt_dual_steepest_edge_column.index[0] = row_out;
    alt_dual_steepest_edge_column.array[row_out] = 1;
    alt_dual_steepest_edge_column.packFlag = false;
    simplex_nla_.btranInScaledSpace(alt_dual_steepest_edge_column,
                                    info_.row_ep_density,
                                    analysis_.pointer_serial_factor_clocks);
    simplex_nla_.ftranInScaledSpace(alt_dual_steepest_edge_column,
                                    info_.row_DSE_density,
                                    analysis_.pointer_serial_factor_clocks);
    // Compute the pivotal column in the scaled space otherwise to check
    //
    // Need \bar{B}^{-1}(R.aq.cq) = \bar{B}^{-1}R.(cq.aq)
    //
    alt_pivotal_column.setup(num_row);
    alt_pivotal_column.clear();
    //
    // Determine cq, and apply it in forming RHS
    //
    lp_.a_matrix_.collectAj(alt_pivotal_column, variable_in, col_aq_scale);
    simplex_nla_.applyBasisMatrixRowScale(alt_pivotal_column);
    simplex_nla_.ftranInScaledSpace(alt_pivotal_column, info_.col_aq_density,
                                    analysis_.pointer_serial_factor_clocks);
    double max_dse_column_error = 0;
    double sum_dse_column_error = 0;
    HighsInt num_dse_column_error = 0;
    const double dse_column_value_tolerance = 1e-2;
    const double dse_column_error_tolerance = 1e-4;
    HighsInt DSE_array_count = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      const double dual_steepest_edge_array_value =
          dual_steepest_edge_array[iRow] * inv_col_ap_scale;
      if (dual_steepest_edge_array_value) DSE_array_count++;
      if (std::abs(dual_steepest_edge_array_value) >
              dse_column_value_tolerance ||
          std::abs(alt_dual_steepest_edge_column.array[iRow]) >
              dse_column_value_tolerance) {
        const double dse_column_error =
            std::abs(alt_dual_steepest_edge_column.array[iRow] -
                     dual_steepest_edge_array_value) /
            std::max(1.0, std::abs(dual_steepest_edge_array_value));
        sum_dse_column_error += dse_column_error;
        if (dse_column_error > dse_column_error_tolerance) {
          max_dse_column_error =
              std::max(dse_column_error, max_dse_column_error);
          num_dse_column_error++;
        }
      }
    }
    if (max_dse_column_error > dse_column_error_tolerance) {
      printf(
          "HEkk::updateDualSteepestEdgeWeights: Iter %2d has num / max / sum = "
          "%d / %g / %g DSE column errors exceeding = %g\n",
          (int)iteration_count_, (int)num_dse_column_error,
          max_dse_column_error, sum_dse_column_error,
          dse_column_error_tolerance);
      printf("DSE column count alt = %d; og = %d)\n",
             (int)alt_dual_steepest_edge_column.count, (int)DSE_array_count);
      for (HighsInt iRow = 0; iRow < num_row; iRow++) {
        const double dual_steepest_edge_array_value =
            dual_steepest_edge_array[iRow] * inv_col_ap_scale;
        if (alt_dual_steepest_edge_column.array[iRow] != 0 &&
            dual_steepest_edge_array_value != 0) {
          const double dse_column_error =
              std::abs(alt_dual_steepest_edge_column.array[iRow] -
                       dual_steepest_edge_array_value) /
              std::max(1.0, std::abs(dual_steepest_edge_array_value));
          if (dse_column_error > 1e-10)
            printf(
                "Row %4d: DSE column (alt = %11.4g; og = %11.4g) difference "
                "%10.4g\n",
                (int)iRow, alt_dual_steepest_edge_column.array[iRow],
                dual_steepest_edge_array_value, dse_column_error);
        }
      }
      fflush(stdout);
      assert(max_dse_column_error < dse_column_error_tolerance);
    }
  }

  if ((HighsInt)dual_edge_weight_.size() < num_row) {
    printf(
        "HEkk::updateDualSteepestEdgeWeights solve %d: "
        "dual_edge_weight_.size() = %d < %d\n",
        (int)debug_solve_call_num_, (int)dual_edge_weight_.size(),
        (int)num_row);
    fflush(stdout);
  }
  assert((HighsInt)dual_edge_weight_.size() >= num_row);
  HighsInt to_entry;
  const bool use_row_indices =
      simplex_nla_.sparseLoopStyle(column_count, num_row, to_entry);
  const bool convert_to_scaled_space = !simplex_in_scaled_space_;
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow = use_row_indices ? variable_index[iEntry] : iEntry;
    double aa_iRow = column_array[iRow];
    if (!aa_iRow) continue;
    double dual_steepest_edge_array_value = dual_steepest_edge_array[iRow];
    if (convert_to_scaled_space) {
      double basic_col_scale = simplex_nla_.basicColScaleFactor(iRow);
      aa_iRow /= basic_col_scale;
      aa_iRow *= col_aq_scale;
      dual_steepest_edge_array_value *= inv_col_ap_scale;
    }
    if (DSE_check) {
      const double pivotal_column_error =
          std::abs(aa_iRow - alt_pivotal_column.array[iRow]);
      if (pivotal_column_error > 1e-4) {
        printf(
            "HEkk::updateDualSteepestEdgeWeights Row %2d of pivotal column has "
            "error %10.4g\n",
            (int)iRow, pivotal_column_error);
        fflush(stdout);
      }
      assert(pivotal_column_error < 1e-4);
    }
    dual_edge_weight_[iRow] += aa_iRow * (new_pivotal_edge_weight * aa_iRow +
                                          Kai * dual_steepest_edge_array_value);
    dual_edge_weight_[iRow] =
        std::max(kMinDualSteepestEdgeWeight, dual_edge_weight_[iRow]);
  }
  analysis_.simplexTimerStop(DseUpdateWeightClock);
}

// Update the Devex weights
void HEkk::updateDualDevexWeights(const HVector* column,
                                  const double new_pivotal_edge_weight) {
  analysis_.simplexTimerStart(DevexUpdateWeightClock);

  const HighsInt num_row = lp_.num_row_;
  const HighsInt column_count = column->count;
  const HighsInt* variable_index = column->index.data();
  const double* column_array = column->array.data();

  if ((HighsInt)dual_edge_weight_.size() < num_row) {
    printf(
        "HEkk::updateDualDevexWeights solve %d: "
        "dual_edge_weight_.size() = %d < %d\n",
        (int)debug_solve_call_num_, (int)dual_edge_weight_.size(),
        (int)num_row);
    fflush(stdout);
  }
  assert((HighsInt)dual_edge_weight_.size() >= num_row);
  HighsInt to_entry;
  const bool use_row_indices =
      simplex_nla_.sparseLoopStyle(column_count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow = use_row_indices ? variable_index[iEntry] : iEntry;
    const double aa_iRow = column_array[iRow];
    dual_edge_weight_[iRow] = max(dual_edge_weight_[iRow],
                                  new_pivotal_edge_weight * aa_iRow * aa_iRow);
  }
  analysis_.simplexTimerStop(DevexUpdateWeightClock);
}

void HEkk::resetSyntheticClock() {
  this->build_synthetic_tick_ = this->simplex_nla_.build_synthetic_tick_;
  this->total_synthetic_tick_ = 0;
}

void HEkk::initialisePartitionedRowwiseMatrix() {
  if (status_.has_ar_matrix) return;
  analysis_.simplexTimerStart(matrixSetupClock);
  ar_matrix_.createRowwisePartitioned(lp_.a_matrix_,
                                      basis_.nonbasicFlag_.data());
  ar_matrix_.name_ = "partitioned row-wise constraint";
  // The partitioned row-wise matrix is generated from an accepted
  // constraint matrix, so it is OK so long as it's been constructed
  // correctly.
  assert(ar_matrix_.legalDimensions(options_->log_options));
  assert(ar_matrix_.trimDimensions(options_->log_options) == HighsStatus::kOk);
  // In particular, only its structure should be assessed, as its
  // values are inherited from the constraint matrix, and may contain
  // small values from presolve
  assert(ar_matrix_.assessStructure(*options_) == HighsStatus::kOk);
  assert(ar_matrix_.debugPartitionOk(basis_.nonbasicFlag_.data()));
  analysis_.simplexTimerStop(matrixSetupClock);
  status_.has_ar_matrix = true;
}

bool HEkk::lpFactorRowCompatible() const {
  return lpFactorRowCompatible(this->lp_.num_row_);
}

bool HEkk::lpFactorRowCompatible(const HighsInt expectedNumRow) const {
  // Check for LP-HFactor row compatibility
  const bool consistent_num_row =
      this->simplex_nla_.factor_.num_row == expectedNumRow;
  if (!consistent_num_row) {
    highsLogDev(options_->log_options, HighsLogType::kError,
                "HEkk::initialiseSimplexLpBasisAndFactor: LP(%6d, %6d) has "
                "factor_num_row = %d\n",
                (int)this->lp_.num_col_, expectedNumRow,
                (int)this->simplex_nla_.factor_.num_row);
  }
  return consistent_num_row;
}

std::string HEkk::simplexStrategyToString(
    const HighsInt simplex_strategy) const {
  assert(kSimplexStrategyMin <= simplex_strategy &&
         simplex_strategy <= kSimplexStrategyMax);
  if (simplex_strategy == kSimplexStrategyChoose)
    return "choose simplex solver";
  if (simplex_strategy == kSimplexStrategyDual) return "dual simplex solver";
  if (simplex_strategy == kSimplexStrategyDualPlain)
    return "serial dual simplex solver";
  if (simplex_strategy == kSimplexStrategyDualTasks)
    return "parallel dual simplex solver - SIP";
  if (simplex_strategy == kSimplexStrategyDualMulti)
    return "parallel dual simplex solver - PAMI";
  if (simplex_strategy == kSimplexStrategyPrimal)
    return "primal simplex solver";
  return "Unknown";
}

void HEkk::zeroBasicDuals() {
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++)
    info_.workDual_[basis_.basicIndex_[iRow]] = 0;
}

void HEkk::setNonbasicMove() {
  const bool have_solution = false;
  // Don't have a simplex basis since nonbasicMove is not set up.

  // Assign nonbasicMove using as much information as is available
  double lower;
  double upper;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  basis_.nonbasicMove_.resize(num_tot);

  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis_.nonbasicFlag_[iVar]) {
      // Basic variable
      basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      continue;
    }
    // Nonbasic variable
    if (iVar < lp_.num_col_) {
      lower = lp_.col_lower_[iVar];
      upper = lp_.col_upper_[iVar];
    } else {
      HighsInt iRow = iVar - lp_.num_col_;
      lower = -lp_.row_upper_[iRow];
      upper = -lp_.row_lower_[iRow];
    }
    int8_t move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        //
        // Determine the bound to set the value to according to, in order of
        // priority
        //
        // 1. Any solution value
        if (have_solution) {
          double midpoint = 0.5 * (lower + upper);
          double value = info_.workValue_[iVar];
          if (value < midpoint) {
            move = kNonbasicMoveUp;
          } else {
            move = kNonbasicMoveDn;
          }
        }
        // 2. Bound of original LP that is closer to zero
        if (move == kIllegalMoveValue) {
          if (fabs(lower) < fabs(upper)) {
            move = kNonbasicMoveUp;
          } else {
            move = kNonbasicMoveDn;
          }
        }
      } else {
        // Lower (since upper bound is infinite)
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = kNonbasicMoveDn;
    } else {
      // FREE
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis_.nonbasicMove_[iVar] = move;
  }
}

void HEkk::allocateWorkAndBaseArrays() {
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  info_.workCost_.resize(num_tot);
  info_.workDual_.resize(num_tot);
  info_.workShift_.resize(num_tot);

  info_.workLower_.resize(num_tot);
  info_.workUpper_.resize(num_tot);
  info_.workRange_.resize(num_tot);
  info_.workValue_.resize(num_tot);
  info_.workLowerShift_.resize(num_tot);
  info_.workUpperShift_.resize(num_tot);

  // Feel that it should be possible to resize this with in dual
  // solver, and only if Devex is being used, but a pointer to it
  // needs to be set up when constructing HDual
  info_.devex_index_.resize(num_tot);

  info_.baseLower_.resize(lp_.num_row_);
  info_.baseUpper_.resize(lp_.num_row_);
  info_.baseValue_.resize(lp_.num_row_);
}

void HEkk::initialiseLpColBound() {
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    info_.workLower_[iCol] = lp_.col_lower_[iCol];
    info_.workUpper_[iCol] = lp_.col_upper_[iCol];
    info_.workRange_[iCol] = info_.workUpper_[iCol] - info_.workLower_[iCol];
    info_.workLowerShift_[iCol] = 0;
    info_.workUpperShift_[iCol] = 0;
  }
}

void HEkk::initialiseLpRowBound() {
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    HighsInt iCol = lp_.num_col_ + iRow;
    info_.workLower_[iCol] = -lp_.row_upper_[iRow];
    info_.workUpper_[iCol] = -lp_.row_lower_[iRow];
    info_.workRange_[iCol] = info_.workUpper_[iCol] - info_.workLower_[iCol];
    info_.workLowerShift_[iCol] = 0;
    info_.workUpperShift_[iCol] = 0;
  }
}

void HEkk::initialiseCost(const SimplexAlgorithm algorithm,
                          const HighsInt solve_phase, const bool perturb) {
  // Copy the cost
  initialiseLpColCost();
  initialiseLpRowCost();
  info_.costs_shifted = false;
  info_.costs_perturbed = false;
  analysis_.net_num_single_cost_shift = 0;
  // Primal simplex costs are either from the LP or set specially in phase 1
  if (algorithm == SimplexAlgorithm::kPrimal) return;
  // Dual simplex costs are either from the LP or perturbed
  if (!perturb || info_.dual_simplex_cost_perturbation_multiplier == 0) return;
  // Perturb the original costs, scale down if is too big
  HighsInt num_original_nonzero_cost = 0;
  highsLogDev(options_->log_options, HighsLogType::kInfo,
              "Cost perturbation for %s\n", lp_.model_name_.c_str());
  // Of the following, min_abs_cost and sum_abs_cost are only used for
  // logging; only max_abs_cost determines perturbation
  double min_abs_cost = kHighsInf;
  double max_abs_cost = 0;
  double sum_abs_cost = 0;
  for (HighsInt i = 0; i < lp_.num_col_; i++) {
    const double abs_cost = fabs(info_.workCost_[i]);
    if (abs_cost) {
      num_original_nonzero_cost++;
      min_abs_cost = min(min_abs_cost, abs_cost);
    }
    sum_abs_cost += abs_cost;
    max_abs_cost = max(max_abs_cost, abs_cost);
  }
  const HighsInt pct0 = (100 * num_original_nonzero_cost) / lp_.num_col_;

  highsLogDev(options_->log_options, HighsLogType::kInfo,
              "   Initially have %" HIGHSINT_FORMAT
              " nonzero costs (%3" HIGHSINT_FORMAT "%%)",
              num_original_nonzero_cost, pct0);
  if (num_original_nonzero_cost) {
    highsLogDev(options_->log_options, HighsLogType::kInfo,
                " with min / average / max = %g / %g / %g\n", min_abs_cost,
                sum_abs_cost / num_original_nonzero_cost, max_abs_cost);
  } else {
    max_abs_cost = 1.0;
    highsLogDev(options_->log_options, HighsLogType::kInfo,
                " but perturb as if max cost was 1\n");
  }

  if (max_abs_cost > 100) {
    max_abs_cost = sqrt(sqrt(max_abs_cost));
    highsLogDev(
        options_->log_options, HighsLogType::kInfo,
        "   Large so set max_abs_cost = sqrt(sqrt(max_abs_cost)) = %g\n",
        max_abs_cost);
  }

  // If there are few boxed variables, we will just use simple perturbation
  double boxedRate = 0;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt i = 0; i < num_tot; i++)
    boxedRate += (info_.workRange_[i] < 1e30);
  boxedRate /= num_tot;
  if (boxedRate < 0.01) {
    max_abs_cost = min(max_abs_cost, 1.0);
    highsLogDev(options_->log_options, HighsLogType::kInfo,
                "   Small boxedRate (%g) so set max_abs_cost = "
                "min(max_abs_cost, 1.0) = "
                "%g\n",
                boxedRate, max_abs_cost);
  }
  // Determine the perturbation base
  cost_perturbation_max_abs_cost_ = max_abs_cost;
  cost_perturbation_base_ =
      info_.dual_simplex_cost_perturbation_multiplier * 5e-7 * max_abs_cost;
  highsLogDev(options_->log_options, HighsLogType::kInfo,
              "   Perturbation column base = %g\n", cost_perturbation_base_);

  // Now do the perturbation
  for (HighsInt i = 0; i < lp_.num_col_; i++) {
    double lower = lp_.col_lower_[i];
    double upper = lp_.col_upper_[i];
    double xpert = (1 + info_.numTotRandomValue_[i]) *
                   (fabs(info_.workCost_[i]) + 1) * cost_perturbation_base_;
    // const double previous_cost = info_.workCost_[i];
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free - no perturb
    } else if (upper >= kHighsInf) {  // Lower
      info_.workCost_[i] += xpert;
    } else if (lower <= -kHighsInf) {  // Upper
      info_.workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      info_.workCost_[i] += (info_.workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
  }
  const double row_cost_perturbation_base_ =
      info_.dual_simplex_cost_perturbation_multiplier * 1e-12;
  highsLogDev(options_->log_options, HighsLogType::kInfo,
              "   Perturbation row    base = %g\n",
              row_cost_perturbation_base_);
  for (HighsInt i = lp_.num_col_; i < num_tot; i++) {
    double perturbation2 =
        (0.5 - info_.numTotRandomValue_[i]) * row_cost_perturbation_base_;
    info_.workCost_[i] += perturbation2;
  }
  info_.costs_perturbed = true;
}

void HEkk::initialiseBound(const SimplexAlgorithm algorithm,
                           const HighsInt solve_phase, const bool perturb) {
  initialiseLpColBound();
  initialiseLpRowBound();
  info_.bounds_shifted = false;
  info_.bounds_perturbed = false;
  // Primal simplex bounds are either from the LP or perturbed
  if (algorithm == SimplexAlgorithm::kPrimal) {
    if (!perturb || info_.primal_simplex_bound_perturbation_multiplier == 0)
      return;
    // Perturb the bounds
    // Determine the smallest and largest finite lower/upper bounds
    HighsInt num_col = lp_.num_col_;
    HighsInt num_row = lp_.num_row_;
    HighsInt num_tot = num_col + num_row;
    double min_abs_lower = kHighsInf;
    double max_abs_lower = -1;
    double min_abs_upper = kHighsInf;
    double max_abs_upper = -1;
    for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
      double abs_lower = fabs(info_.workLower_[iVar]);
      double abs_upper = fabs(info_.workUpper_[iVar]);
      if (abs_lower && abs_lower < kHighsInf) {
        min_abs_lower = min(abs_lower, min_abs_lower);
        max_abs_lower = max(abs_lower, max_abs_lower);
      }
      if (abs_upper && abs_upper < kHighsInf) {
        min_abs_upper = min(abs_upper, min_abs_upper);
        max_abs_upper = max(abs_upper, max_abs_upper);
      }
    }
    // printf(
    //     "Nonzero finite lower bounds in [%9.4g, %9.4g]; upper bounds in "
    //     "[%9.4g, %9.4g]\n",
    //     min_abs_lower, max_abs_lower, min_abs_upper, max_abs_upper);

    const double base =
        info_.primal_simplex_bound_perturbation_multiplier * 5e-7;
    for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
      double lower = info_.workLower_[iVar];
      double upper = info_.workUpper_[iVar];
      const bool fixed = lower == upper;
      // Don't perturb bounds of nonbasic fixed variables as they stay nonbasic
      if (basis_.nonbasicFlag_[iVar] == kNonbasicFlagTrue && fixed) continue;
      double random_value = info_.numTotRandomValue_[iVar];
      if (lower > -kHighsInf) {
        if (lower < -1) {
          lower -= random_value * base * (-lower);
        } else if (lower < 1) {
          lower -= random_value * base;
        } else {
          lower -= random_value * base * lower;
        }
        info_.workLower_[iVar] = lower;
      }
      if (upper < kHighsInf) {
        if (upper < -1) {
          upper += random_value * base * (-upper);
        } else if (upper < 1) {
          upper += random_value * base;
        } else {
          upper += random_value * base * upper;
        }
        info_.workUpper_[iVar] = upper;
      }
      info_.workRange_[iVar] = info_.workUpper_[iVar] - info_.workLower_[iVar];
      if (basis_.nonbasicFlag_[iVar] == kNonbasicFlagFalse) continue;
      // Set values of nonbasic variables
      if (basis_.nonbasicMove_[iVar] > 0) {
        info_.workValue_[iVar] = lower;
      } else if (basis_.nonbasicMove_[iVar] < 0) {
        info_.workValue_[iVar] = upper;
      }
    }
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      HighsInt iVar = basis_.basicIndex_[iRow];
      info_.baseLower_[iRow] = info_.workLower_[iVar];
      info_.baseUpper_[iRow] = info_.workUpper_[iVar];
    }
    info_.bounds_perturbed = true;
    return;
  }
  // Dual simplex bounds are either from the LP or set to special values in
  // phase
  // 1
  assert(algorithm == SimplexAlgorithm::kDual);
  if (solve_phase == kSolvePhase2) return;

  // The dual objective is the sum of products of primal and dual
  // values for nonbasic variables. For dual simplex phase 1, the
  // primal bounds are set so that when the dual value is feasible, the
  // primal value is set to zero. Otherwise the value is +1/-1
  // according to the required sign of the dual, except for free
  // variables, where the bounds are [-1000, 1000]. Hence the dual
  // objective is the negation of the sum of infeasibilities, unless there are
  // free In Phase 1: change to dual phase 1 bound.
  const double inf = kHighsInf;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    if (info_.workLower_[iCol] == -inf && info_.workUpper_[iCol] == inf) {
      // Phase 1 bounds were not modified for free rows in hsol. Why?
      // To reduce the number of free variables on the assumption that
      // free rows should never be nonbasic? This is normally true
      // when starting from a logical basis or crash - unless
      // singularity makes a free logical nonbasic - but what about an
      // advanced basis start? It doesn't happen with the HiGHS MIP
      // solver. However, SCIP cut handling can lead to a nonbasic row
      // with nonzero dual being made free and, particularly if it's
      // (then) the only dual infeasibility, phase 1 fails to remove
      // it.
      //
      // if (iCol >= lp_.num_col_) continue;
      info_.workLower_[iCol] = -1000,
      info_.workUpper_[iCol] = 1000;  // FREE
    } else if (info_.workLower_[iCol] == -inf) {
      info_.workLower_[iCol] = -1,
      info_.workUpper_[iCol] = 0;  // UPPER
    } else if (info_.workUpper_[iCol] == inf) {
      info_.workLower_[iCol] = 0,
      info_.workUpper_[iCol] = 1;  // LOWER
    } else {
      info_.workLower_[iCol] = 0,
      info_.workUpper_[iCol] = 0;  // BOXED or FIXED
    }
    info_.workRange_[iCol] = info_.workUpper_[iCol] - info_.workLower_[iCol];
  }
}

void HEkk::initialiseLpColCost() {
  double cost_scale_factor = pow(2.0, options_->cost_scale_factor);
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    info_.workCost_[iCol] =
        (HighsInt)lp_.sense_ * cost_scale_factor * lp_.col_cost_[iCol];
    info_.workShift_[iCol] = 0;
  }
}

void HEkk::initialiseLpRowCost() {
  for (HighsInt iCol = lp_.num_col_; iCol < lp_.num_col_ + lp_.num_row_;
       iCol++) {
    info_.workCost_[iCol] = 0;
    info_.workShift_[iCol] = 0;
  }
}

void HEkk::initialiseNonbasicValueAndMove() {
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis_.nonbasicFlag_[iVar]) {
      // Basic variable
      basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      continue;
    }
    // Nonbasic variable
    const double lower = info_.workLower_[iVar];
    const double upper = info_.workUpper_[iVar];
    const int8_t original_move = basis_.nonbasicMove_[iVar];
    double value;
    int8_t move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      value = lower;
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (original_move == kNonbasicMoveUp) {
          // Set at lower
          value = lower;
          move = kNonbasicMoveUp;
        } else if (original_move == kNonbasicMoveDn) {
          // Set at upper
          value = upper;
          move = kNonbasicMoveDn;
        } else {
          // Invalid nonbasicMove: correct and set value at lower
          value = lower;
          move = kNonbasicMoveUp;
        }
      } else {
        // Lower
        value = lower;
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      value = upper;
      move = kNonbasicMoveDn;
    } else {
      // FREE
      value = 0;
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis_.nonbasicMove_[iVar] = move;
    info_.workValue_[iVar] = value;
  }
}

void HEkk::pivotColumnFtran(const HighsInt iCol, HVector& col_aq) {
  analysis_.simplexTimerStart(FtranClock);
  col_aq.clear();
  col_aq.packFlag = true;
  lp_.a_matrix_.collectAj(col_aq, iCol, 1);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordBefore(kSimplexNlaFtran, col_aq,
                                    info_.col_aq_density);
  simplex_nla_.ftran(col_aq, info_.col_aq_density,
                     analysis_.pointer_serial_factor_clocks);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaFtran, col_aq);
  HighsInt num_row = lp_.num_row_;
  const double local_col_aq_density = (double)col_aq.count / num_row;
  updateOperationResultDensity(local_col_aq_density, info_.col_aq_density);
  analysis_.simplexTimerStop(FtranClock);
}

void HEkk::unitBtran(const HighsInt iRow, HVector& row_ep) {
  analysis_.simplexTimerStart(BtranClock);
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = iRow;
  row_ep.array[iRow] = 1;
  row_ep.packFlag = true;
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordBefore(kSimplexNlaBtranEp, row_ep,
                                    info_.row_ep_density);
  simplex_nla_.btran(row_ep, info_.row_ep_density,
                     analysis_.pointer_serial_factor_clocks);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaBtranEp, row_ep);
  HighsInt num_row = lp_.num_row_;
  const double local_row_ep_density = (double)row_ep.count / num_row;
  updateOperationResultDensity(local_row_ep_density, info_.row_ep_density);
  analysis_.simplexTimerStop(BtranClock);
}

void HEkk::fullBtran(HVector& buffer) {
  // Performs BTRAN on the buffer supplied. Make sure that
  // buffer.count is large (>lp_.num_row_ to be sure) rather
  // than 0 if the indices of the RHS (and true value of buffer.count)
  // isn't known.
  analysis_.simplexTimerStart(BtranFullClock);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordBefore(kSimplexNlaBtranFull, buffer,
                                    info_.dual_col_density);
  simplex_nla_.btran(buffer, info_.dual_col_density,
                     analysis_.pointer_serial_factor_clocks);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaBtranFull, buffer);
  const double local_dual_col_density = (double)buffer.count / lp_.num_row_;
  updateOperationResultDensity(local_dual_col_density, info_.dual_col_density);
  analysis_.simplexTimerStop(BtranFullClock);
}

void HEkk::choosePriceTechnique(const HighsInt price_strategy,
                                const double row_ep_density,
                                bool& use_col_price,
                                bool& use_row_price_w_switch) const {
  // By default switch to column PRICE when pi_p has at least this
  // density
  const double density_for_column_price_switch = 0.75;
  use_col_price = (price_strategy == kSimplexPriceStrategyCol) ||
                  (price_strategy == kSimplexPriceStrategyRowSwitchColSwitch &&
                   row_ep_density > density_for_column_price_switch);
  use_row_price_w_switch =
      price_strategy == kSimplexPriceStrategyRowSwitch ||
      price_strategy == kSimplexPriceStrategyRowSwitchColSwitch;
}

void HEkk::tableauRowPrice(const bool quad_precision, const HVector& row_ep,
                           HVector& row_ap, const HighsInt debug_report) {
  analysis_.simplexTimerStart(PriceClock);
  const HighsInt solver_num_row = lp_.num_row_;
  const HighsInt solver_num_col = lp_.num_col_;
  const double local_density = 1.0 * row_ep.count / solver_num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  choosePriceTechnique(info_.price_strategy, local_density, use_col_price,
                       use_row_price_w_switch);
  if (analysis_.analyse_simplex_summary_data) {
    if (use_col_price) {
      const double expected_density = 1;
      analysis_.operationRecordBefore(kSimplexNlaPriceAp, row_ep,
                                      expected_density);
      analysis_.num_col_price++;
    } else if (use_row_price_w_switch) {
      analysis_.operationRecordBefore(kSimplexNlaPriceAp, row_ep,
                                      info_.row_ep_density);
      analysis_.num_row_price_with_switch++;
    } else {
      analysis_.operationRecordBefore(kSimplexNlaPriceAp, row_ep,
                                      info_.row_ep_density);
      analysis_.num_row_price++;
    }
  }
  row_ap.clear();
  if (use_col_price) {
    // Perform column-wise PRICE
    lp_.a_matrix_.priceByColumn(quad_precision, row_ap, row_ep, debug_report);
  } else if (use_row_price_w_switch) {
    // Perform hyper-sparse row-wise PRICE, but switch if the density of row_ap
    // becomes extreme
    const double switch_density = kHyperPriceDensity;
    ar_matrix_.priceByRowWithSwitch(quad_precision, row_ap, row_ep,
                                    info_.row_ap_density, 0, switch_density,
                                    debug_report);
  } else {
    // Perform hyper-sparse row-wise PRICE
    ar_matrix_.priceByRow(quad_precision, row_ap, row_ep, debug_report);
  }
  if (use_col_price) {
    // Column-wise PRICE computes components corresponding to basic
    // variables, so zero these by exploiting the fact that, for basic
    // variables, nonbasicFlag[*]=0
    const int8_t* nonbasicFlag = basis_.nonbasicFlag_.data();
    for (HighsInt iCol = 0; iCol < solver_num_col; iCol++)
      row_ap.array[iCol] *= nonbasicFlag[iCol];
  }
  // Update the record of average row_ap density
  const double local_row_ap_density = (double)row_ap.count / solver_num_col;
  updateOperationResultDensity(local_row_ap_density, info_.row_ap_density);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaPriceAp, row_ap);
  analysis_.simplexTimerStop(PriceClock);
}

void HEkk::fullPrice(const HVector& full_col, HVector& full_row) {
  analysis_.simplexTimerStart(PriceFullClock);
  full_row.clear();
  if (analysis_.analyse_simplex_summary_data) {
    const double expected_density = 1;
    analysis_.operationRecordBefore(kSimplexNlaPriceFull, full_col,
                                    expected_density);
  }
  const bool quad_precision = false;
  lp_.a_matrix_.priceByColumn(quad_precision, full_row, full_col);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaPriceFull, full_row);
  analysis_.simplexTimerStop(PriceFullClock);
}

void HEkk::computePrimal() {
  analysis_.simplexTimerStart(ComputePrimalClock);
  const HighsInt num_row = lp_.num_row_;
  const HighsInt num_col = lp_.num_col_;
  // Setup a local buffer for the values of basic variables
  HVector primal_col;
  primal_col.setup(num_row);
  primal_col.clear();
  for (HighsInt i = 0; i < num_col + num_row; i++) {
    if (basis_.nonbasicFlag_[i] && info_.workValue_[i] != 0) {
      lp_.a_matrix_.collectAj(primal_col, i, info_.workValue_[i]);
    }
  }
  // It's possible that the buffer has no nonzeros, so performing
  // FTRAN is unnecessary. Not much of a saving, but the zero density
  // looks odd in the analysis!
  if (primal_col.count) {
    simplex_nla_.ftran(primal_col, info_.primal_col_density,
                       analysis_.pointer_serial_factor_clocks);
    const double local_primal_col_density = (double)primal_col.count / num_row;
    updateOperationResultDensity(local_primal_col_density,
                                 info_.primal_col_density);
  }
  for (HighsInt i = 0; i < num_row; i++) {
    HighsInt iCol = basis_.basicIndex_[i];
    info_.baseValue_[i] = -primal_col.array[i];
    info_.baseLower_[i] = info_.workLower_[iCol];
    info_.baseUpper_[i] = info_.workUpper_[iCol];
  }
  // Indicate that the primal infeasibility information isn't known
  info_.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
  info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;

  analysis_.simplexTimerStop(ComputePrimalClock);
}

void HEkk::computeDual() {
  analysis_.simplexTimerStart(ComputeDualClock);
  // Create a local buffer for the pi vector
  HVector dual_col;
  dual_col.setup(lp_.num_row_);
  dual_col.clear();
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    const double value = info_.workCost_[basis_.basicIndex_[iRow]] +
                         info_.workShift_[basis_.basicIndex_[iRow]];
    if (value) {
      dual_col.index[dual_col.count++] = iRow;
      dual_col.array[iRow] = value;
    }
  }
  // If debugging, save the current duals
  const bool debug_compute_dual = false;
  if (debug_compute_dual) {
    debugComputeDual(true);
    debugSimplexDualInfeasible("(old duals)", true);
  }
  // Copy the costs in case the basic costs are all zero
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt i = 0; i < num_tot; i++)
    info_.workDual_[i] = info_.workCost_[i] + info_.workShift_[i];

  if (dual_col.count) {
    fullBtran(dual_col);
    // Create a local buffer for the values of reduced costs
    HVector dual_row;
    dual_row.setup(lp_.num_col_);
    fullPrice(dual_col, dual_row);
    for (HighsInt i = 0; i < lp_.num_col_; i++)
      info_.workDual_[i] -= dual_row.array[i];
    for (HighsInt i = lp_.num_col_; i < num_tot; i++)
      info_.workDual_[i] -= dual_col.array[i - lp_.num_col_];
    if (debug_compute_dual) {
      debugComputeDual();
      debugSimplexDualInfeasible("(new duals)", true);
    }
  }
  // Indicate that the dual infeasibility information isn't known
  info_.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;

  analysis_.simplexTimerStop(ComputeDualClock);
}

double HEkk::computeDualForTableauColumn(const HighsInt iVar,
                                         const HVector& tableau_column) const {
  const vector<double>& workCost = info_.workCost_;
  const vector<HighsInt>& basicIndex = basis_.basicIndex_;

  double dual = info_.workCost_[iVar];
  for (HighsInt i = 0; i < tableau_column.count; i++) {
    HighsInt iRow = tableau_column.index[i];
    dual -= tableau_column.array[iRow] * workCost[basicIndex[iRow]];
  }
  return dual;
}

bool HEkk::reinvertOnNumericalTrouble(
    const std::string& method_name, double& numerical_trouble_measure,
    const double alpha_from_col, const double alpha_from_row,
    const double numerical_trouble_tolerance) {
  double abs_alpha_from_col = fabs(alpha_from_col);
  double abs_alpha_from_row = fabs(alpha_from_row);
  double min_abs_alpha = min(abs_alpha_from_col, abs_alpha_from_row);
  double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  numerical_trouble_measure = abs_alpha_diff / min_abs_alpha;
  const HighsInt update_count = info_.update_count;
  // Reinvert if the relative difference is large enough, and updates have been
  // performed
  const bool numerical_trouble =
      numerical_trouble_measure > numerical_trouble_tolerance;
  const bool reinvert = numerical_trouble && update_count > 0;
  debugReportReinvertOnNumericalTrouble(method_name, numerical_trouble_measure,
                                        alpha_from_col, alpha_from_row,
                                        numerical_trouble_tolerance, reinvert);
  if (reinvert) {
    // Consider increasing the Markowitz multiplier
    const double current_pivot_threshold = info_.factor_pivot_threshold;
    double new_pivot_threshold = 0;
    if (current_pivot_threshold < kDefaultPivotThreshold) {
      // Threshold is below default value, so increase it
      new_pivot_threshold =
          min(current_pivot_threshold * kPivotThresholdChangeFactor,
              kDefaultPivotThreshold);
    } else if (current_pivot_threshold < kMaxPivotThreshold) {
      // Threshold is below max value, so increase it if few updates have been
      // performed
      if (update_count < 10)
        new_pivot_threshold =
            min(current_pivot_threshold * kPivotThresholdChangeFactor,
                kMaxPivotThreshold);
    }
    if (new_pivot_threshold) {
      highsLogUser(options_->log_options, HighsLogType::kWarning,
                   "   Increasing Markowitz threshold to %g\n",
                   new_pivot_threshold);
      info_.factor_pivot_threshold = new_pivot_threshold;
      simplex_nla_.setPivotThreshold(new_pivot_threshold);
    }
  }
  return reinvert;
}

// The major model updates. Factor calls factor_.update; Matrix
// calls matrix_.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void HEkk::transformForUpdate(HVector* column, HVector* row_ep,
                              const HighsInt variable_in, HighsInt* row_out) {
  simplex_nla_.transformForUpdate(column, row_ep, variable_in, *row_out);
}

void HEkk::flipBound(const HighsInt iCol) {
  const int8_t move = basis_.nonbasicMove_[iCol] = -basis_.nonbasicMove_[iCol];
  info_.workValue_[iCol] =
      move == 1 ? info_.workLower_[iCol] : info_.workUpper_[iCol];
}

void HEkk::updateFactor(HVector* column, HVector* row_ep, HighsInt* iRow,
                        HighsInt* hint) {
  analysis_.simplexTimerStart(UpdateFactorClock);
  simplex_nla_.update(column, row_ep, iRow, hint);
  // Now have a representation of B^{-1}, but it is not fresh
  status_.has_invert = true;
  if (info_.update_count >= info_.update_limit)
    *hint = kRebuildReasonUpdateLimitReached;

  // Determine whether to reinvert based on the synthetic clock
  bool reinvert_syntheticClock =
      this->total_synthetic_tick_ >= this->build_synthetic_tick_;
  const bool performed_min_updates =
      info_.update_count >= kSyntheticTickReinversionMinUpdateCount;
  if (reinvert_syntheticClock && performed_min_updates)
    *hint = kRebuildReasonSyntheticClockSaysInvert;
  analysis_.simplexTimerStop(UpdateFactorClock);
  // Use the next level down for the debug level, since the cost of
  // checking the INVERT every iteration is an order more expensive
  // than checking after factorization.
  HighsInt alt_debug_level = options_->highs_debug_level - 1;
  // Forced expensive debug for development work
  //  if (debug_solve_report_) alt_debug_level = kHighsDebugLevelExpensive;
  HighsDebugStatus debug_status =
      debugNlaCheckInvert("HEkk::updateFactor", alt_debug_level);
  if (debug_status == HighsDebugStatus::kError) {
    *hint = kRebuildReasonPossiblySingularBasis;
  }
}

void HEkk::updatePivots(const HighsInt variable_in, const HighsInt row_out,
                        const HighsInt move_out) {
  analysis_.simplexTimerStart(UpdatePivotsClock);
  HighsInt variable_out = basis_.basicIndex_[row_out];

  // update hash value of basis
  HighsHashHelpers::sparse_inverse_combine(basis_.hash, variable_out);
  HighsHashHelpers::sparse_combine(basis_.hash, variable_in);
  visited_basis_.insert(basis_.hash);

  // Incoming variable
  basis_.basicIndex_[row_out] = variable_in;
  basis_.nonbasicFlag_[variable_in] = 0;
  basis_.nonbasicMove_[variable_in] = 0;
  info_.baseLower_[row_out] = info_.workLower_[variable_in];
  info_.baseUpper_[row_out] = info_.workUpper_[variable_in];

  // Outgoing variable
  basis_.nonbasicFlag_[variable_out] = 1;
  if (info_.workLower_[variable_out] == info_.workUpper_[variable_out]) {
    info_.workValue_[variable_out] = info_.workLower_[variable_out];
    basis_.nonbasicMove_[variable_out] = 0;
  } else if (move_out == -1) {
    info_.workValue_[variable_out] = info_.workLower_[variable_out];
    basis_.nonbasicMove_[variable_out] = 1;
  } else {
    info_.workValue_[variable_out] = info_.workUpper_[variable_out];
    basis_.nonbasicMove_[variable_out] = -1;
  }
  // Update the dual objective value
  double nwValue = info_.workValue_[variable_out];
  double vrDual = info_.workDual_[variable_out];
  double dl_dual_objective_value = nwValue * vrDual;
  info_.updated_dual_objective_value += dl_dual_objective_value;
  info_.update_count++;
  // Update the number of basic logicals
  if (variable_out < lp_.num_col_) info_.num_basic_logicals++;
  if (variable_in < lp_.num_col_) info_.num_basic_logicals--;
  // No longer have a representation of B^{-1}, and certainly not
  // fresh!
  status_.has_invert = false;
  status_.has_fresh_invert = false;
  // Data are no longer fresh from rebuild
  status_.has_fresh_rebuild = false;
  analysis_.simplexTimerStop(UpdatePivotsClock);
}

bool HEkk::isBadBasisChange(const SimplexAlgorithm algorithm,
                            const HighsInt variable_in, const HighsInt row_out,
                            const HighsInt rebuild_reason) {
  if (rebuild_reason) return false;
  if (variable_in == -1 || row_out == -1) return false;
  uint64_t currhash = basis_.hash;
  HighsInt variable_out = basis_.basicIndex_[row_out];

  HighsHashHelpers::sparse_inverse_combine(currhash, variable_out);
  HighsHashHelpers::sparse_combine(currhash, variable_in);

  bool cycling_detected = false;
  const bool posible_cycling = visited_basis_.find(currhash) != nullptr;
  if (posible_cycling) {
    if (iteration_count_ == previous_iteration_cycling_detected + 1) {
      // Cycling detected on successive iterations suggests infinite cycling
      //      highsLogDev(options_->log_options, HighsLogType::kWarning,
      //		  "Cycling detected in %s simplex:");
      // printf("Cycling detected in %s simplex solve %d (Iteration %d)",
      //        algorithm == SimplexAlgorithm::kPrimal ? "primal" : "dual",
      //        (int)debug_solve_call_num_, (int)iteration_count_);
      cycling_detected = true;
    } else {
      previous_iteration_cycling_detected = iteration_count_;
    }
  }
  if (cycling_detected) {
    if (algorithm == SimplexAlgorithm::kDual) {
      analysis_.num_dual_cycling_detections++;
    } else {
      analysis_.num_primal_cycling_detections++;
    }
    highsLogDev(options_->log_options, HighsLogType::kWarning,
                " basis change (%d out; %d in) is bad\n", (int)variable_out,
                (int)variable_in);
    addBadBasisChange(row_out, variable_out, variable_in,
                      BadBasisChangeReason::kCycling, true);
    return true;
  } else {
    // Look to see whether this basis change is in the list of bad
    // ones
    for (auto& change : bad_basis_change_) {
      if (change.variable_out == variable_out &&
          change.variable_in == variable_in && change.row_out == row_out) {
        change.taboo = true;
        return true;
      }
    }
  }

  return false;
}

void HEkk::updateMatrix(const HighsInt variable_in,
                        const HighsInt variable_out) {
  analysis_.simplexTimerStart(UpdateMatrixClock);
  ar_matrix_.update(variable_in, variable_out, lp_.a_matrix_);
  //  assert(ar_matrix_.debugPartitionOk(basis_.nonbasicFlag_.data()));
  analysis_.simplexTimerStop(UpdateMatrixClock);
}

void HEkk::computeInfeasibilitiesForReporting(const SimplexAlgorithm algorithm,
                                              const HighsInt solve_phase) {
  if (algorithm == SimplexAlgorithm::kPrimal) {
    // Report the primal and dual infeasibilities
    computeSimplexInfeasible();
  } else {
    // Report the primal infeasibilities
    computeSimplexPrimalInfeasible();
    if (solve_phase == kSolvePhase1) {
      // In phase 1, report the simplex LP dual infeasibilities
      computeSimplexLpDualInfeasible();
    } else {
      // In phase 2, report the simplex dual infeasibilities
      computeSimplexDualInfeasible();
    }
  }
}

void HEkk::computeSimplexInfeasible() {
  computeSimplexPrimalInfeasible();
  computeSimplexDualInfeasible();
}

void HEkk::computeSimplexPrimalInfeasible() {
  // Computes num/max/sum of primal infeasibilities according to the
  // simplex bounds. This is used to determine optimality in dual
  // phase 1 and dual phase 2, albeit using different bounds in
  // workLower/Upper.
  analysis_.simplexTimerStart(ComputePrIfsClock);
  const double scaled_primal_feasibility_tolerance =
      options_->primal_feasibility_tolerance;
  HighsInt& num_primal_infeasibility = info_.num_primal_infeasibilities;
  double& max_primal_infeasibility = info_.max_primal_infeasibility;
  double& sum_primal_infeasibility = info_.sum_primal_infeasibilities;
  num_primal_infeasibility = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibility = 0;

  for (HighsInt i = 0; i < lp_.num_col_ + lp_.num_row_; i++) {
    if (basis_.nonbasicFlag_[i]) {
      // Nonbasic column
      double value = info_.workValue_[i];
      double lower = info_.workLower_[i];
      double upper = info_.workUpper_[i];
      // @primal_infeasibility calculation
      double primal_infeasibility = 0;
      if (value < lower - scaled_primal_feasibility_tolerance) {
        primal_infeasibility = lower - value;
      } else if (value > upper + scaled_primal_feasibility_tolerance) {
        primal_infeasibility = value - upper;
      }
      if (primal_infeasibility > 0) {
        if (primal_infeasibility > scaled_primal_feasibility_tolerance)
          num_primal_infeasibility++;
        max_primal_infeasibility =
            std::max(primal_infeasibility, max_primal_infeasibility);
        sum_primal_infeasibility += primal_infeasibility;
      }
    }
  }
  for (HighsInt i = 0; i < lp_.num_row_; i++) {
    // Basic variable
    double value = info_.baseValue_[i];
    double lower = info_.baseLower_[i];
    double upper = info_.baseUpper_[i];
    // @primal_infeasibility calculation
    double primal_infeasibility = 0;
    if (value < lower - scaled_primal_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + scaled_primal_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > 0) {
      if (primal_infeasibility > scaled_primal_feasibility_tolerance)
        num_primal_infeasibility++;
      max_primal_infeasibility =
          std::max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibility += primal_infeasibility;
    }
  }
  analysis_.simplexTimerStop(ComputePrIfsClock);
}

void HEkk::computeSimplexDualInfeasible() {
  analysis_.simplexTimerStart(ComputeDuIfsClock);
  // Computes num/max/sum of dual infeasibilities in phase 1 and phase
  // 2 according to nonbasicMove. The bounds are only used to identify
  // free variables. Fixed variables are assumed to have
  // nonbasicMove=0 so that no dual infeasibility is counted for them.
  const double scaled_dual_feasibility_tolerance =
      options_->dual_feasibility_tolerance;
  HighsInt& num_dual_infeasibility = info_.num_dual_infeasibilities;
  double& max_dual_infeasibility = info_.max_dual_infeasibility;
  double& sum_dual_infeasibility = info_.sum_dual_infeasibilities;
  num_dual_infeasibility = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibility = 0;

  for (HighsInt iCol = 0; iCol < lp_.num_col_ + lp_.num_row_; iCol++) {
    if (!basis_.nonbasicFlag_[iCol]) continue;
    // Nonbasic column
    const double dual = info_.workDual_[iCol];
    const double lower = info_.workLower_[iCol];
    const double upper = info_.workUpper_[iCol];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else {
      // Not free: any dual infeasibility is given by the dual value
      // signed by nonbasicMove
      dual_infeasibility = -basis_.nonbasicMove_[iCol] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance) {
        num_dual_infeasibility++;
      }
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  analysis_.simplexTimerStop(ComputeDuIfsClock);
}

void HEkk::computeSimplexLpDualInfeasible() {
  // Compute num/max/sum of dual infeasibilities according to the
  // bounds of the simplex LP. Assumes that boxed variables have
  // primal variable at the bound corresponding to the sign of the
  // dual so should only be used in dual phase 1 - where it's only
  // used for reporting after rebuilds.
  const double scaled_dual_feasibility_tolerance =
      options_->dual_feasibility_tolerance;
  HighsInt& num_dual_infeasibility =
      analysis_.num_dual_phase_1_lp_dual_infeasibility;
  double& max_dual_infeasibility =
      analysis_.max_dual_phase_1_lp_dual_infeasibility;
  double& sum_dual_infeasibility =
      analysis_.sum_dual_phase_1_lp_dual_infeasibility;
  num_dual_infeasibility = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibility = 0;

  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    HighsInt iVar = iCol;
    if (!basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double dual = info_.workDual_[iVar];
    const double lower = lp_.col_lower_[iCol];
    const double upper = lp_.col_upper_[iCol];
    double dual_infeasibility = 0;
    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free: any nonzero dual value is infeasible
        dual_infeasibility = fabs(dual);
      } else {
        // Only lower bounded: a negative dual is infeasible
        dual_infeasibility = -dual;
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded: a positive dual is infeasible
        dual_infeasibility = dual;
      } else {
        // Boxed or fixed: any dual value is feasible
        dual_infeasibility = 0;
      }
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    HighsInt iVar = lp_.num_col_ + iRow;
    if (!basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic row
    const double dual = -info_.workDual_[iVar];
    const double lower = lp_.row_lower_[iRow];
    const double upper = lp_.row_upper_[iRow];
    double dual_infeasibility = 0;
    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free: any nonzero dual value is infeasible
        dual_infeasibility = fabs(dual);
      } else {
        // Only lower bounded: a negative dual is infeasible
        dual_infeasibility = -dual;
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded: a positive dual is infeasible
        dual_infeasibility = dual;
      } else {
        // Boxed or fixed: any dual value is feasible
        dual_infeasibility = 0;
      }
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
}

void HEkk::invalidatePrimalMaxSumInfeasibilityRecord() {
  info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
}

void HEkk::invalidatePrimalInfeasibilityRecord() {
  info_.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
  invalidatePrimalMaxSumInfeasibilityRecord();
}

void HEkk::invalidateDualMaxSumInfeasibilityRecord() {
  info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
}

void HEkk::invalidateDualInfeasibilityRecord() {
  info_.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  invalidateDualMaxSumInfeasibilityRecord();
}

bool HEkk::bailout() {
  if (solve_bailout_) {
    // Bailout has already been decided: check that it's for one of these
    // reasons
    assert(model_status_ == HighsModelStatus::kTimeLimit ||
           model_status_ == HighsModelStatus::kIterationLimit ||
           model_status_ == HighsModelStatus::kObjectiveBound ||
           model_status_ == HighsModelStatus::kObjectiveTarget);
  } else if (options_->time_limit < kHighsInf &&
             timer_->read() > options_->time_limit) {
    solve_bailout_ = true;
    model_status_ = HighsModelStatus::kTimeLimit;
  } else if (iteration_count_ >= options_->simplex_iteration_limit) {
    solve_bailout_ = true;
    model_status_ = HighsModelStatus::kIterationLimit;
  } else if (callback_->user_callback &&
             callback_->active[kCallbackSimplexInterrupt]) {
    callback_->clearHighsCallbackOutput();
    callback_->data_out.simplex_iteration_count = iteration_count_;
    if (callback_->callbackAction(kCallbackSimplexInterrupt,
                                  "Simplex interrupt")) {
      highsLogDev(options_->log_options, HighsLogType::kInfo,
                  "User interrupt\n");
      solve_bailout_ = true;
      model_status_ = HighsModelStatus::kInterrupt;
    }
  }
  return solve_bailout_;
}

HighsStatus HEkk::returnFromEkkSolve(const HighsStatus return_status) {
  if (analysis_.analyse_simplex_time)
    analysis_.simplexTimerStop(SimplexTotalClock);
  // Restore any modified development or timing settings and analyse
  // solver timing
  if (debug_solve_report_) debugReporting(1);
  if (time_report_) timeReporting(1);
  // Note that in timeReporting(1), analysis_.analyse_simplex_time
  // reverts to its value given by options_
  if (analysis_.analyse_simplex_time) analysis_.reportSimplexTimer();
  simplex_stats_.valid = true;
  // Since HEkk::iteration_count_ includes iteration on presolved LP,
  // simplex_stats_.iteration_count is initialised to -
  // HEkk::iteration_count_
  simplex_stats_.iteration_count += iteration_count_;
  //  simplex_stats_.num_invert is incremented internally
  simplex_stats_.last_invert_num_el = simplex_nla_.factor_.invert_num_el;
  simplex_stats_.last_factored_basis_num_el =
      simplex_nla_.factor_.basis_matrix_num_el;
  simplex_stats_.col_aq_density = analysis_.col_aq_density;
  simplex_stats_.row_ep_density = analysis_.row_ep_density;
  simplex_stats_.row_ap_density = analysis_.row_ap_density;
  simplex_stats_.row_DSE_density = analysis_.row_DSE_density;
  return return_status;
}

HighsStatus HEkk::returnFromSolve(const HighsStatus return_status) {
  // Always called before returning from HEkkPrimal/Dual::solve()
  if (solve_bailout_) {
    // If bailout has already been decided: check that it's for one of
    // these reasons
    assert(model_status_ == HighsModelStatus::kTimeLimit ||
           model_status_ == HighsModelStatus::kIterationLimit ||
           model_status_ == HighsModelStatus::kObjectiveBound ||
           model_status_ == HighsModelStatus::kObjectiveTarget ||
           model_status_ == HighsModelStatus::kInterrupt);
  }
  // Check that returnFromSolve has not already been called: it should
  // be called exactly once per solve
  assert(!called_return_from_solve_);
  called_return_from_solve_ = true;
  info_.valid_backtracking_basis_ = false;

  // Initialise the status of the primal and dual solutions
  return_primal_solution_status_ = kSolutionStatusNone;
  return_dual_solution_status_ = kSolutionStatusNone;
  // Nothing more is known about the solve after an error return
  if (return_status == HighsStatus::kError) return return_status;

  // Check that an invert exists
  assert(status_.has_invert);

  // Determine a primal and dual solution, removing the effects of
  // perturbations and shifts
  //
  // Unless the solution is optimal, invalidate the infeasibility data
  if (model_status_ != HighsModelStatus::kOptimal) {
    invalidatePrimalInfeasibilityRecord();
    invalidateDualInfeasibilityRecord();
  }
  // The simplex algorithm used on exit should be set
  assert(exit_algorithm_ != SimplexAlgorithm::kNone);
  switch (model_status_) {
    case HighsModelStatus::kOptimal: {
      if (info_.num_primal_infeasibilities) {
        // Optimal - but not to desired primal feasibility tolerance
        return_primal_solution_status_ = kSolutionStatusInfeasible;
      } else {
        return_primal_solution_status_ = kSolutionStatusFeasible;
      }
      if (info_.num_dual_infeasibilities) {
        // Optimal - but not to desired dual feasibility tolerance
        return_dual_solution_status_ = kSolutionStatusInfeasible;
      } else {
        return_dual_solution_status_ = kSolutionStatusFeasible;
      }
      break;
    }
    case HighsModelStatus::kInfeasible: {
      // Primal infeasibility has been identified in primal phase 1,
      // or proved in dual phase 2. There should be no primal
      // perturbations or shifts
      assert(!info_.bounds_shifted);
      assert(!info_.bounds_perturbed);
      if (exit_algorithm_ == SimplexAlgorithm::kPrimal) {
        // Reset the simplex costs and recompute duals after primal
        // phase 1
        initialiseCost(SimplexAlgorithm::kDual, kSolvePhase2);
        computeDual();
      }
      computeSimplexInfeasible();
      // Primal solution shouldn't be feasible
      assert(info_.num_primal_infeasibilities > 0);
      break;
    }
    case HighsModelStatus::kUnboundedOrInfeasible: {
      // Dual simplex has identified dual infeasibility in phase
      // 1. There should be no dual perturbations
      assert(exit_algorithm_ == SimplexAlgorithm::kDual);
      assert(!info_.costs_shifted);
      assert(!info_.costs_perturbed);
      // Reset the simplex bounds and recompute primals
      initialiseBound(SimplexAlgorithm::kDual, kSolvePhase2);
      computePrimal();
      computeSimplexInfeasible();
      // Dual solution shouldn't be feasible
      assert(info_.num_dual_infeasibilities > 0);
      break;
    }
    case HighsModelStatus::kUnbounded: {
      // Primal simplex has identified unboundedness in phase 2. There
      // should be no primal or dual perturbations
      assert(exit_algorithm_ == SimplexAlgorithm::kPrimal);
      assert(!info_.costs_shifted);
      assert(!info_.costs_perturbed);
      assert(!info_.bounds_shifted);
      assert(!info_.bounds_perturbed);
      computeSimplexInfeasible();
      // Primal solution should be feasible
      assert(info_.num_primal_infeasibilities == 0);
      break;
    }
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
    case HighsModelStatus::kTimeLimit:
    case HighsModelStatus::kIterationLimit:
    case HighsModelStatus::kInterrupt:
    case HighsModelStatus::kUnknown: {
      // Simplex has failed to conclude a model property. Either it
      // has bailed out due to reaching the objective bound, target,
      // time, iteration limit or user interrupt, or it has not been
      // set (cycling is the only reason). Could happen anywhere.
      //
      // Reset the simplex bounds and recompute primals
      initialiseBound(SimplexAlgorithm::kDual, kSolvePhase2);
      initialiseNonbasicValueAndMove();
      computePrimal();
      // Reset the simplex costs and recompute duals
      initialiseCost(SimplexAlgorithm::kDual, kSolvePhase2);
      computeDual();
      computeSimplexInfeasible();
      break;
    }
    default: {
      highsLogDev(
          options_->log_options, HighsLogType::kError,
          "%s simplex solver returns status %s\n",
          exit_algorithm_ == SimplexAlgorithm::kPrimal ? "primal" : "dual",
          utilModelStatusToString(model_status_).c_str());
      return HighsStatus::kError;
      break;
    }
  }
  this->zeroBasicDuals();
  assert(info_.num_primal_infeasibilities >= 0);
  assert(info_.num_dual_infeasibilities >= 0);
  if (info_.num_primal_infeasibilities == 0) {
    return_primal_solution_status_ = kSolutionStatusFeasible;
  } else {
    return_primal_solution_status_ = kSolutionStatusInfeasible;
  }
  if (info_.num_dual_infeasibilities == 0) {
    return_dual_solution_status_ = kSolutionStatusFeasible;
  } else {
    return_dual_solution_status_ = kSolutionStatusInfeasible;
  }
  assert(debugNoShiftsOrPerturbations());
  computePrimalObjectiveValue();
  if (!options_->log_dev_level) {
    const bool force = true;
    analysis_.userInvertReport(force);
  }
  return return_status;
}

double HEkk::computeBasisCondition(const HighsLp& lp, const bool exact,
                                   const bool report) const {
  HighsInt solver_num_row = lp.num_row_;
  HighsInt solver_num_col = lp.num_col_;
  vector<double> bs_cond_x;
  vector<double> bs_cond_y;
  vector<double> bs_cond_z;
  vector<double> bs_cond_w;
  HVector row_ep;
  row_ep.setup(solver_num_row);

  const HighsInt* Astart = lp.a_matrix_.start_.data();
  const double* Avalue = lp.a_matrix_.value_.data();
  double exact_norm_Binv = 0;
  if (exact) {
    // Compute the exact norm of B^{-1}
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
      row_ep.clear();
      row_ep.index[row_ep.count] = r_n;
      row_ep.array[r_n] = 1.0;
      row_ep.count++;
      row_ep.packFlag = false;
      simplex_nla_.ftran(row_ep, 0.1);
      assert(row_ep.count <= solver_num_row);
      double c_norm = 0.0;
      for (HighsInt iX = 0; iX < row_ep.count; iX++)
        c_norm += std::fabs(row_ep.array[row_ep.index[iX]]);
      exact_norm_Binv = std::max(c_norm, exact_norm_Binv);
    }
  }
  // Compute the Hager condition number estimate for the basis matrix
  const double expected_density = 1;
  bs_cond_x.resize(solver_num_row);
  bs_cond_y.resize(solver_num_row);
  bs_cond_z.resize(solver_num_row);
  bs_cond_w.resize(solver_num_row);
  // x = ones(n,1)/n;
  // y = A\x;
  double mu = 1.0 / solver_num_row;
  double norm_Binv = 0;
  for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = mu;
  row_ep.clear();
  for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
    double value = bs_cond_x[r_n];
    if (value) {
      row_ep.index[row_ep.count] = r_n;
      row_ep.array[r_n] = value;
      row_ep.count++;
    }
  }
  for (HighsInt ps_n = 1; ps_n <= 5; ps_n++) {
    row_ep.packFlag = false;
    simplex_nla_.ftran(row_ep, expected_density);

    // zeta = sign(y);
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
      bs_cond_y[r_n] = row_ep.array[r_n];
      if (bs_cond_y[r_n] > 0)
        bs_cond_w[r_n] = 1.0;
      else if (bs_cond_y[r_n] < 0)
        bs_cond_w[r_n] = -1.0;
      else
        bs_cond_w[r_n] = 0.0;
    }
    // z=A'\zeta;
    row_ep.clear();
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
      double value = bs_cond_w[r_n];
      if (value) {
        row_ep.index[row_ep.count] = r_n;
        row_ep.array[r_n] = value;
        row_ep.count++;
      }
    }
    row_ep.packFlag = false;
    simplex_nla_.btran(row_ep, expected_density);
    double norm_z = 0.0;
    double ztx = 0.0;
    norm_Binv = 0.0;
    HighsInt argmax_z = -1;
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
      bs_cond_z[r_n] = row_ep.array[r_n];
      double abs_z_v = fabs(bs_cond_z[r_n]);
      if (abs_z_v > norm_z) {
        norm_z = abs_z_v;
        argmax_z = r_n;
      }
      ztx += bs_cond_z[r_n] * bs_cond_x[r_n];
      norm_Binv += fabs(bs_cond_y[r_n]);
    }
    if (norm_z <= ztx) break;
    // x = zeros(n,1);
    // x(fd_i) = 1;
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = 0.0;
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = argmax_z;
    row_ep.array[argmax_z] = 1.0;
    bs_cond_x[argmax_z] = 1.0;
  }
  double norm_B = 0.0;
  for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
    HighsInt vr_n = basis_.basicIndex_[r_n];
    double c_norm = 0.0;
    if (vr_n < solver_num_col)
      for (HighsInt el_n = Astart[vr_n]; el_n < Astart[vr_n + 1]; el_n++)
        c_norm += fabs(Avalue[el_n]);
    else
      c_norm += 1.0;
    norm_B = max(c_norm, norm_B);
  }
  double cond_B = norm_Binv * norm_B;
  double exact_cond_B = exact_norm_Binv * norm_B;
  if (exact) {
    assert(exact_norm_Binv > 0);
    if (report)
      highsLogUser(
          options_->log_options, HighsLogType::kInfo,
          "HEkk::computeBasisCondition: grep_kappa model,||B||_1,approx "
          "||B^{-1}||_1,approx_kappa,||B^{-1}||_1,kappa = ,%s,%g,%g,%g,%g,%g\n",
          lp.model_name_.c_str(), norm_B, norm_Binv, cond_B, exact_norm_Binv,
          exact_cond_B);
    return exact_cond_B;
  }
  return cond_B;
}

void HEkk::initialiseAnalysis() {
  analysis_.setup(lp_name_, lp_, *options_, iteration_count_);
}

std::string HEkk::rebuildReason(const HighsInt rebuild_reason) const {
  std::string rebuild_reason_string;
  if (rebuild_reason == kRebuildReasonCleanup) {
    rebuild_reason_string = "Perturbation cleanup";
  } else if (rebuild_reason == kRebuildReasonNo) {
    rebuild_reason_string = "No reason";
  } else if (rebuild_reason == kRebuildReasonUpdateLimitReached) {
    rebuild_reason_string = "Update limit reached";
  } else if (rebuild_reason == kRebuildReasonSyntheticClockSaysInvert) {
    rebuild_reason_string = "Synthetic clock";
  } else if (rebuild_reason == kRebuildReasonPossiblyOptimal) {
    rebuild_reason_string = "Possibly optimal";
  } else if (rebuild_reason == kRebuildReasonPossiblyPhase1Feasible) {
    rebuild_reason_string = "Possibly phase 1 feasible";
  } else if (rebuild_reason == kRebuildReasonPossiblyPrimalUnbounded) {
    rebuild_reason_string = "Possibly primal unbounded";
  } else if (rebuild_reason == kRebuildReasonPossiblyDualUnbounded) {
    rebuild_reason_string = "Possibly dual unbounded";
  } else if (rebuild_reason == kRebuildReasonPossiblySingularBasis) {
    rebuild_reason_string = "Possibly singular basis";
  } else if (rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex) {
    rebuild_reason_string = "Primal infeasible in primal simplex";
  } else if (rebuild_reason == kRebuildReasonChooseColumnFail) {
    rebuild_reason_string = "Choose column failure";
  } else if (rebuild_reason == kRebuildReasonExcessivePrimalValue) {
    rebuild_reason_string = "Excessive primal value";
  } else {
    rebuild_reason_string = "Unidentified";
    assert(1 == 0);
  }
  return rebuild_reason_string;
}

void HEkk::putIterate() {
  assert(this->status_.has_invert);
  SimplexIterate& iterate = this->simplex_nla_.simplex_iterate_;
  this->simplex_nla_.putInvert();
  iterate.basis_ = this->basis_;
  if (this->status_.has_dual_steepest_edge_weights) {
    // Copy the dual edge weights
    iterate.dual_edge_weight_ = this->dual_edge_weight_;
  } else {
    // Clear to indicate no weights
    iterate.dual_edge_weight_.clear();
  }
}

HighsStatus HEkk::getIterate() {
  SimplexIterate& iterate = this->simplex_nla_.simplex_iterate_;
  if (!iterate.valid_) return HighsStatus::kError;
  this->simplex_nla_.getInvert();
  this->basis_ = iterate.basis_;
  if (iterate.dual_edge_weight_.size()) {
    this->dual_edge_weight_ = iterate.dual_edge_weight_;
  } else {
    this->status_.has_dual_steepest_edge_weights = false;
  }
  this->status_.has_invert = true;
  return HighsStatus::kOk;
}

double HEkk::factorSolveError() {
  // Cheap assessment of factor accuracy.
  //
  // Forms a random solution with at most 50 nonzeros, solves for
  // the corresponding RHS, and then checks the 50 solution values.
  const HighsInt num_col = this->lp_.num_col_;
  const HighsInt num_row = this->lp_.num_row_;
  const HighsSparseMatrix& a_matrix = this->lp_.a_matrix_;
  const vector<HighsInt>& basic_index = this->basis_.basicIndex_;
  const HighsSparseMatrix& ar_matrix = this->ar_matrix_;
  HVector btran_rhs;
  HVector ftran_rhs;
  btran_rhs.setup(num_row);
  ftran_rhs.setup(num_row);

  // Solve for a random solution
  HighsRandom random(1);

  ftran_rhs.clear();
  const HighsInt ideal_solution_num_nz = 50;
  HighsInt solution_num_nz = min(ideal_solution_num_nz, (num_row + 1) / 2);
  assert(solution_num_nz > 0);
  vector<double> solution_value;
  vector<HighsInt> solution_index;
  vector<int8_t> solution_nonzero;
  solution_nonzero.assign(num_row, 0);
  for (;;) {
    HighsInt iRow = random.integer(num_row);
    assert(iRow < num_row);
    if (solution_nonzero[iRow]) continue;
    double value = random.fraction();
    solution_value.push_back(value);
    solution_index.push_back(iRow);
    solution_nonzero[iRow] = 1;
    HighsInt iCol = basic_index[iRow];
    a_matrix.collectAj(ftran_rhs, iCol, value);
    if ((int)solution_value.size() == solution_num_nz) break;
  }

  btran_rhs.clear();
  vector<double> btran_solution;
  btran_solution.assign(num_row, 0);
  for (size_t iX = 0; iX < solution_value.size(); iX++)
    btran_solution[solution_index[iX]] = solution_value[iX];
  vector<double> btran_scattered_rhs;
  btran_scattered_rhs.assign(num_col + num_row, 0);
  for (size_t iX = 0; iX < solution_value.size(); iX++) {
    HighsInt iRow = solution_index[iX];
    for (HighsInt iEl = ar_matrix.p_end_[iRow];
         iEl < ar_matrix.start_[iRow + 1]; iEl++) {
      HighsInt iCol = ar_matrix.index_[iEl];
      btran_scattered_rhs[iCol] += ar_matrix.value_[iEl] * solution_value[iX];
    }
    HighsInt iCol = num_col + iRow;
    if (this->basis_.nonbasicFlag_[iCol] == 0)
      btran_scattered_rhs[iCol] = solution_value[iX];
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iCol = basic_index[iRow];
    if (btran_scattered_rhs[iCol] == 0) continue;
    btran_rhs.array[iRow] = btran_scattered_rhs[iCol];
    btran_rhs.index[btran_rhs.count++] = iRow;
  }

  const double expected_density = solution_num_nz * info_.col_aq_density;
  ftran(ftran_rhs, expected_density);
  btran(btran_rhs, expected_density);

  double ftran_solution_error = 0;
  for (size_t iX = 0; iX < solution_value.size(); iX++)
    ftran_solution_error =
        max(fabs(ftran_rhs.array[solution_index[iX]] - solution_value[iX]),
            ftran_solution_error);
  double btran_solution_error = 0;
  for (size_t iX = 0; iX < solution_value.size(); iX++)
    btran_solution_error =
        max(fabs(btran_rhs.array[solution_index[iX]] - solution_value[iX]),
            btran_solution_error);
  double solution_error = max(ftran_solution_error, btran_solution_error);
  return solution_error;
}

void HEkk::clearBadBasisChange(const BadBasisChangeReason reason) {
  if (reason == BadBasisChangeReason::kAll) {
    bad_basis_change_.clear();
  } else {
    bad_basis_change_.erase(
        std::remove_if(
            bad_basis_change_.begin(), bad_basis_change_.end(),
            [reason](const HighsSimplexBadBasisChangeRecord& record) {
              return record.reason == reason;
            }),
        bad_basis_change_.end());
  }
}

void HEkk::updateBadBasisChange(const HVector& col_aq, double theta_primal) {
  if (!bad_basis_change_.empty()) {
    bad_basis_change_.erase(
        std::remove_if(bad_basis_change_.begin(), bad_basis_change_.end(),
                       [&](const HighsSimplexBadBasisChangeRecord& record) {
                         return std::fabs(col_aq.array[record.row_out] *
                                          theta_primal) >=
                                options_->primal_feasibility_tolerance;
                       }),
        bad_basis_change_.end());
  }
}

HighsInt HEkk::addBadBasisChange(const HighsInt row_out,
                                 const HighsInt variable_out,
                                 const HighsInt variable_in,
                                 const BadBasisChangeReason reason,
                                 const bool taboo) {
  assert(0 <= row_out && row_out <= lp_.num_row_);
  assert(0 <= variable_out && variable_out <= lp_.num_col_ + lp_.num_row_);
  if (variable_in >= 0) {
    assert(0 <= variable_in && variable_in <= lp_.num_col_ + lp_.num_row_);
  } else {
    assert(variable_in == -1);
  }
  // Check that this is not already on the list
  const HighsInt num_bad_basis_change = bad_basis_change_.size();
  HighsInt bad_basis_change_num = -1;
  for (HighsInt Ix = 0; Ix < num_bad_basis_change; Ix++) {
    HighsSimplexBadBasisChangeRecord& record = bad_basis_change_[Ix];
    if (record.row_out == row_out && record.variable_out == variable_out &&
        record.variable_in == variable_in && record.reason == reason) {
      bad_basis_change_num = Ix;
      break;
    }
  }
  if (bad_basis_change_num < 0) {
    // Not on the list so create record
    HighsSimplexBadBasisChangeRecord record;
    record.taboo = taboo;
    record.row_out = row_out;
    record.variable_out = variable_out;
    record.variable_in = variable_in;
    record.reason = reason;
    bad_basis_change_.push_back(record);
    bad_basis_change_num = bad_basis_change_.size() - 1;
  } else {
    // On the list so just update whether it is taboo
    bad_basis_change_[bad_basis_change_num].taboo = taboo;
  }
  return bad_basis_change_num;
}

void HEkk::clearBadBasisChangeTabooFlag() {
  for (auto& change : bad_basis_change_) change.taboo = false;
}

bool HEkk::tabooBadBasisChange() const {
  for (const auto& change : bad_basis_change_) {
    if (change.taboo) return true;
  }
  return false;
}

void HEkk::applyTabooRowOut(vector<double>& values,
                            const double overwrite_with) {
  assert(values.size() >= static_cast<size_t>(lp_.num_row_));
  for (auto& change : bad_basis_change_) {
    if (change.taboo) {
      HighsInt iRow = change.row_out;
      change.save_value = values[iRow];
      values[iRow] = overwrite_with;
    }
  }
}

void HEkk::unapplyTabooRowOut(vector<double>& values) {
  assert((HighsInt)values.size() >= lp_.num_row_);
  // Unapply taboo rows in opposite order in case the row appears
  // twice in the list. This way the first saved value for the row is
  // what remains, not overwrite_with
  for (HighsInt iX = (HighsInt)bad_basis_change_.size() - 1; iX >= 0; iX--) {
    if (bad_basis_change_[iX].taboo)
      values[bad_basis_change_[iX].row_out] = bad_basis_change_[iX].save_value;
  }
}

void HEkk::applyTabooVariableIn(vector<double>& values,
                                const double overwrite_with) {
  assert(values.size() >=
         static_cast<size_t>(lp_.num_col_) + static_cast<size_t>(lp_.num_row_));
  for (auto& change : bad_basis_change_) {
    if (change.taboo) {
      HighsInt iCol = change.variable_in;
      change.save_value = values[iCol];
      values[iCol] = overwrite_with;
    }
  }
}

void HEkk::unapplyTabooVariableIn(vector<double>& values) {
  assert((HighsInt)values.size() >= lp_.num_col_ + lp_.num_row_);
  // Unapply taboo variables in opposite order in case the row appears
  // twice in the list. This way the first saved value for the
  // variable is what remains, not overwrite_with
  for (HighsInt iX = (HighsInt)bad_basis_change_.size() - 1; iX >= 0; iX--) {
    if (bad_basis_change_[iX].taboo)
      values[bad_basis_change_[iX].variable_in] =
          bad_basis_change_[iX].save_value;
  }
}

bool HEkk::logicalBasis() const {
  for (HighsInt iRow = 0; iRow < this->lp_.num_row_; iRow++) {
    if (basis_.basicIndex_[iRow] < this->lp_.num_col_) return false;
  }
  return true;
}

bool HEkk::proofOfPrimalInfeasibility() {
  // To be called from outside HEkk when row_ep is not known
  assert(dual_ray_record_.index >= 0);
  HighsLp& lp = this->lp_;
  HighsInt move_out = dual_ray_record_.sign;
  HighsInt row_out = dual_ray_record_.index;
  // Compute the basis inverse row
  HVector row_ep;
  row_ep.setup(lp.num_row_);
  unitBtran(row_out, row_ep);
  return proofOfPrimalInfeasibility(row_ep, move_out, row_out);
}

bool HEkk::proofOfPrimalInfeasibility(HVector& row_ep, const HighsInt move_out,
                                      const HighsInt row_out) {
  // To be called from inside HEkkDual
  HighsLp& lp = this->lp_;

  HighsInt debug_product_report = kDebugReportOff;
  const bool debug_proof_report_on = true;
  bool debug_rows_report = false;
  bool debug_proof_report = false;
  if (debug_iteration_report_) {
    if (debug_proof_report_on)
      debug_product_report = kDebugReportOff;  // kDebugReportAll; //
    debug_rows_report = debug_proof_report_on;
    debug_proof_report = debug_proof_report_on;
  }

  const bool use_row_wise_matrix = status_.has_ar_matrix;
  const bool use_iterative_refinement = false;  // debug_iteration_report_;//
  if (use_iterative_refinement) {
    simplex_nla_.reportArray("Row e_p.0", lp.num_col_, &row_ep, true);
    unitBtranIterativeRefinement(row_out, row_ep);
    simplex_nla_.reportArray("Row e_p.1", lp.num_col_, &row_ep, true);
  }

  // Refine row_ep by removing relatively small values
  double row_ep_scale = 0;
  // if (use_refinement) refineArray(row_ep, row_ep_scale,
  // refinement_tolerance);
  // Determine the maximum absolute value in row_ep
  HighsCDouble proof_lower = 0.0;
  const HighsInt max_num_basic_proof_report = 25;
  const HighsInt max_num_zeroed_report = 25;
  HighsInt num_zeroed_for_small_report = 0;
  double max_zeroed_for_small_value = 0;
  HighsInt num_zeroed_for_lb_report = 0;
  double max_zeroed_for_lb_value = 0;
  HighsInt num_zeroed_for_ub_report = 0;
  double max_zeroed_for_ub_value = 0;
  for (HighsInt iX = 0; iX < row_ep.count; iX++) {
    HighsInt iRow = row_ep.index[iX];
    // Give row_ep the sign of the leaving row - as is done in
    // getDualRayInterface.
    const double row_ep_value = row_ep.array[iRow];
    assert(row_ep_value);
    if (std::abs(row_ep_value * getMaxAbsRowValue(iRow)) <=
        options_->small_matrix_value) {
      if (debug_proof_report) {
        const double abs_row_ep_value = std::abs(row_ep_value);
        if (num_zeroed_for_small_report < max_num_zeroed_report &&
            max_zeroed_for_small_value < abs_row_ep_value) {
          printf(
              "Zeroed row_ep.array[%6d] = %11.4g due to being small in "
              "contribution\n",
              (int)iRow, row_ep_value);
          num_zeroed_for_small_report++;
          max_zeroed_for_small_value = abs_row_ep_value;
        }
      }
      row_ep.array[iRow] = 0.0;
      continue;
    }

    row_ep.array[iRow] *= move_out;

    // make sure infinite sides are not used
    double rowBound;
    if (row_ep.array[iRow] > 0) {
      rowBound = lp.row_lower_[iRow];
      if (highs_isInfinity(-rowBound)) {
        // row lower bound is infinite
        if (debug_proof_report) {
          const double abs_row_ep_value = std::abs(row_ep_value);
          if (num_zeroed_for_lb_report < max_num_zeroed_report &&
              max_zeroed_for_lb_value < abs_row_ep_value) {
            printf(
                "Zeroed row_ep.array[%6d] = %11.4g due to infinite lower "
                "bound\n",
                (int)iRow, row_ep_value);
            num_zeroed_for_lb_report++;
            max_zeroed_for_lb_value = abs_row_ep_value;
          }
        }
        row_ep.array[iRow] = 0.0;
        continue;
      }

    } else {
      rowBound = lp.row_upper_[iRow];
      if (highs_isInfinity(rowBound)) {
        // row upper bound is infinite
        if (debug_proof_report) {
          const double abs_row_ep_value = std::abs(row_ep_value);
          if (num_zeroed_for_ub_report < max_num_zeroed_report &&
              max_zeroed_for_ub_value < abs_row_ep_value) {
            printf(
                "Zeroed row_ep.array[%6d] = %11.4g due to infinite upper "
                "bound\n",
                (int)iRow, row_ep_value);
            num_zeroed_for_ub_report++;
            max_zeroed_for_ub_value = abs_row_ep_value;
          }
        }
        row_ep.array[iRow] = 0.0;
        continue;
      }
    }
    // add up lower bound of proof constraint
    proof_lower += row_ep.array[iRow] * rowBound;
  }
  // Form the proof constraint coefficients
  proof_value_.clear();
  proof_index_.clear();
  vector<double>& proof_value = this->proof_value_;
  vector<HighsInt>& proof_index = this->proof_index_;
  if (use_row_wise_matrix) {
    this->ar_matrix_.productTransposeQuad(proof_value, proof_index, row_ep,
                                          debug_product_report);
  } else {
    lp.a_matrix_.productTransposeQuad(proof_value, proof_index, row_ep,
                                      debug_product_report);
  }

  HighsInt proof_num_nz = proof_index.size();
  if (debug_rows_report) {
    simplex_nla_.reportArray("Row e_p", lp.num_col_, &row_ep, true);
    simplex_nla_.reportVector("Proof", proof_num_nz, proof_value, proof_index,
                              true);
  }
  if (debug_proof_report) {
    printf(
        "HEkk::proofOfPrimalInfeasibility row_ep.count = %d; proof_num_nz = "
        "%d; row_ep_scale = %g\n",
        (int)row_ep.count, (int)proof_num_nz, row_ep_scale);
    HighsInt num_basic_proof_report = 0;
    double max_basic_proof_value = 0;
    for (HighsInt i = 0; i < proof_num_nz; ++i) {
      const HighsInt iCol = proof_index[i];
      const double value = proof_value[i];
      const double abs_value = std::abs(value);
      if (!basis_.nonbasicFlag_[iCol] && max_basic_proof_value < abs_value &&
          num_basic_proof_report < max_num_basic_proof_report) {
        printf("Proof entry %6d (Column %6d) is basic with value %11.4g\n",
               (int)i, (int)iCol, value);
        max_basic_proof_value = abs_value;
        num_basic_proof_report++;
      }
    }
  }
  HighsCDouble implied_upper = 0.0;
  HighsCDouble sumInf = 0.0;
  for (HighsInt i = 0; i < proof_num_nz; ++i) {
    const HighsInt iCol = proof_index[i];
    const double value = proof_value[i];
    if (value > 0) {
      if (highs_isInfinity(lp.col_upper_[iCol])) {
        sumInf += value;
        if (sumInf > options_->small_matrix_value) break;
        continue;
        // Commented out unreachable code
        //        if (value <= options_->small_matrix_value) continue;
      }
      implied_upper += value * lp.col_upper_[iCol];
    } else {
      if (highs_isInfinity(-lp.col_lower_[iCol])) {
        sumInf += -value;
        if (sumInf > options_->small_matrix_value) break;
        continue;
      }
      implied_upper += value * lp.col_lower_[iCol];
    }
  }
  bool infinite_implied_upper = sumInf > options_->small_matrix_value;
  const double gap = double(proof_lower - implied_upper);
  const bool gap_ok = gap > options_->primal_feasibility_tolerance;
  const bool proof_of_primal_infeasibility = !infinite_implied_upper && gap_ok;

  const double local_report = false;
  if (!proof_of_primal_infeasibility && local_report) {
    printf(
        "HEkk::proofOfPrimalInfeasibility: row %6d; gap = %11.4g (%s); "
        "sumInf = %11.4g (%s) so proof is %s\n",
        (int)row_out, gap, highsBoolToString(gap_ok).c_str(), (double)sumInf,
        highsBoolToString(infinite_implied_upper).c_str(),
        highsBoolToString(proof_of_primal_infeasibility).c_str());
  }
  if (debug_proof_report) {
    printf("HEkk::proofOfPrimalInfeasibility has %sfinite implied upper bound",
           infinite_implied_upper ? "in" : "");
    if (!infinite_implied_upper) printf(" and gap = %g", gap);
    printf(" so proof is %s\n",
           proof_of_primal_infeasibility ? "true" : "false");
  }
  return proof_of_primal_infeasibility;
}

double HEkk::getValueScale(const HighsInt count,
                           const vector<double>& value) const {
  if (count <= 0) return 1;
  double max_abs_value = 0;
  for (HighsInt iX = 0; iX < count; iX++)
    max_abs_value = std::max(fabs(value[iX]), max_abs_value);
  return nearestPowerOfTwoScale(max_abs_value);
}

double HEkk::getMaxAbsRowValue(HighsInt row) {
  if (!status_.has_ar_matrix) initialisePartitionedRowwiseMatrix();

  double val = -1.0;
  for (HighsInt i = ar_matrix_.start_[row]; i < ar_matrix_.start_[row + 1]; ++i)
    val = std::max(val, std::abs(ar_matrix_.value_[i]));

  return val;
}

void HEkk::unitBtranIterativeRefinement(const HighsInt row_out,
                                        HVector& row_ep) {
  // Perform an iteration of refinement
  HighsLp& lp = this->lp_;
  HVector residual;
  double residual_norm = 0;
  double correction_norm = 0;
  const double expected_density = 1;
  residual.setup(lp.num_row_);
  unitBtranResidual(row_out, row_ep, residual, residual_norm);
  const bool debug_iterative_refinement_report_on = false;
  bool debug_iterative_refinement_report = false;
  if (debug_iteration_report_) {
    debug_iterative_refinement_report = debug_iterative_refinement_report_on;
  }
  if (debug_iterative_refinement_report)
    printf(
        "HEkk::unitBtranIterativeRefinement: Residual   has %6d / %6d nonzeros "
        "and norm of %g\n",
        (int)residual.count, (int)lp.num_row_, residual_norm);
  if (!residual_norm) return;
  // Normalise using nearest power of 2 to ||correction_rhs|| so kHighsTiny
  // isn't used adversely
  const double residual_scale = nearestPowerOfTwoScale(residual_norm);
  for (HighsInt iEl = 0; iEl < residual.count; iEl++)
    residual.array[residual.index[iEl]] *= residual_scale;
  btran(residual, expected_density);
  row_ep.count = 0;
  correction_norm = 0;
  // Adding two (possibly sparse) vectors, so have to loop over all rows
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (residual.array[iRow]) {
      const double correction_value = residual.array[iRow] / residual_scale;
      correction_norm = max(fabs(correction_value), correction_norm);
      row_ep.array[iRow] -= correction_value;
    }
    if (fabs(row_ep.array[iRow]) < kHighsTiny) {
      row_ep.array[iRow] = 0;
    } else {
      row_ep.index[row_ep.count++] = iRow;
    }
  }
  if (debug_iterative_refinement_report)
    printf(
        "HEkk::unitBtranIterativeRefinement: Correction has %6d / %6d nonzeros "
        "and norm of %g\n",
        (int)residual.count, (int)lp.num_row_, correction_norm);
}

void HEkk::unitBtranResidual(const HighsInt row_out, const HVector& row_ep,
                             HVector& residual, double& residual_norm) {
  HighsLp& lp = this->lp_;
  vector<HighsCDouble> quad_residual;
  quad_residual.assign(lp.num_row_, 0);
  quad_residual[row_out] = -1.0;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    HighsInt iVar = basis_.basicIndex_[iRow];
    HighsCDouble value = quad_residual[iRow];
    if (iVar < lp.num_col_) {
      for (HighsInt iEl = lp.a_matrix_.start_[iVar];
           iEl < lp.a_matrix_.start_[iVar + 1]; iEl++)
        value +=
            lp.a_matrix_.value_[iEl] * row_ep.array[lp.a_matrix_.index_[iEl]];
    } else {
      value += row_ep.array[iVar - lp.num_col_];
    }
    quad_residual[iRow] = value;
  }
  residual.clear();
  residual.packFlag = false;
  residual_norm = 0;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    const double value = (double)quad_residual[iRow];
    if (value) {
      residual.array[iRow] = value;
      residual.index[residual.count++] = iRow;
    }
    residual_norm = max(fabs(residual.array[iRow]), residual_norm);
  }
}

void HighsSimplexStats::report(FILE* file, const std::string& message) const {
  fprintf(file, "\nSimplex stats: %s\n", message.c_str());
  fprintf(file, "   valid                      = %d\n", this->valid);
  fprintf(file, "   iteration_count            = %d\n",
          static_cast<int>(this->iteration_count));
  fprintf(file, "   num_invert                 = %d\n",
          static_cast<int>(this->num_invert));
  fprintf(file, "   last_invert_num_el         = %d\n",
          static_cast<int>(this->last_invert_num_el));
  fprintf(file, "   last_factored_basis_num_el = %d\n",
          static_cast<int>(this->last_factored_basis_num_el));
  fprintf(file, "   col_aq_density             = %g\n", this->col_aq_density);
  fprintf(file, "   row_ep_density             = %g\n", this->row_ep_density);
  fprintf(file, "   row_ap_density             = %g\n", this->row_ap_density);
  fprintf(file, "   row_DSE_density            = %g\n", this->row_DSE_density);
}

void HighsSimplexStats::initialise(const HighsInt iteration_count_) {
  valid = false;
  iteration_count = -iteration_count_;
  num_invert = 0;
  last_invert_num_el = 0;
  last_factored_basis_num_el = 0;
  col_aq_density = 0;
  row_ep_density = 0;
  row_ap_density = 0;
  row_DSE_density = 0;
}

HighsRayRecord HighsRayRecord::getRayRecord() const {
  HighsRayRecord record;
  record.index = this->index;
  record.sign = this->sign;
  record.value = this->value;
  return record;
}

void HighsRayRecord::setRayRecord(const HighsRayRecord& from_record) {
  this->index = from_record.index;
  this->sign = from_record.sign;
  this->value = from_record.value;
}

void HighsRayRecord::clear() {
  this->index = kNoRayIndex;
  this->sign = kNoRaySign;
  this->value.clear();
}

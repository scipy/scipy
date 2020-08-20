/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplex.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "simplex/HSimplex.h"

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsSolutionDebug.h"
#include "lp_data/HighsStatus.h"
#include "simplex/HCrash.h"
#include "simplex/HFactorDebug.h"
#include "simplex/HSimplexDebug.h"
#include "simplex/HVector.h"
#include "simplex/HighsSimplexInterface.h"
#include "simplex/SimplexConst.h"  // For simplex strategy constants
#include "simplex/SimplexTimer.h"
#include "util/HighsSort.h"
#include "util/HighsUtils.h"

using std::runtime_error;
#include <cassert>
#include <vector>

#ifdef OPENMP
#include "omp.h"
#endif

void setSimplexOptions(HighsModelObject& highs_model_object) {
  const HighsOptions& options = highs_model_object.options_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  //  HighsSolutionParams& scaled_solution_params =
  //  highs_model_object.scaled_solution_params_;
  //
  // Copy values of HighsOptions for the simplex solver
  //
  // Currently most of these options are straight copies, but they
  // will become valuable when "choose" becomes a HiGHS strategy value
  // that will need converting into a specific simplex strategy value.
  //
  simplex_info.simplex_strategy = options.simplex_strategy;
  simplex_info.dual_edge_weight_strategy =
      options.simplex_dual_edge_weight_strategy;
  simplex_info.price_strategy = options.simplex_price_strategy;
  simplex_info.dual_simplex_cost_perturbation_multiplier =
      options.dual_simplex_cost_perturbation_multiplier;
  simplex_info.update_limit = options.simplex_update_limit;

  // Set values of internal options
  simplex_info.store_squared_primal_infeasibility = true;
  // Option for analysing the LP solution
#ifdef HiGHSDEV
  bool useful_analysis = false;  // true;  //
  bool full_timing = false;
  // Options for reporting timing
  simplex_info.report_simplex_inner_clock = true;  // useful_analysis;
  simplex_info.report_simplex_outer_clock = full_timing;
  simplex_info.report_simplex_phases_clock = full_timing;
  simplex_info.report_HFactor_clock = useful_analysis;  // full_timing;//
  // Options for analysing the LP and simplex iterations
  simplex_info.analyse_lp = useful_analysis;  // false;  //
  simplex_info.analyse_iterations = useful_analysis;
  //  simplex_info.analyse_invert_form = useful_analysis;
  //  simplex_info.analyse_invert_condition = useful_analysis;
  simplex_info.analyse_invert_time = full_timing;
  simplex_info.analyse_rebuild_time = full_timing;
#endif
}

HighsStatus transition(HighsModelObject& highs_model_object) {
  // Perform the transition from whatever information is known about
  // the LP to a status where simplex data are set up for the initial
  // rebuild() of the chosen solver - primal, scalar dual or parallel
  // dual.
  //
  // First look at what basis and solution information is known. If a
  // simplex basis is known, then it's used, and there's no need for
  // the solution values. This will usually correspond to hot start
  // when solving MIP problems. Otherwise, generate a simplex basis,
  // thus:
  //
  // If there is a HiGHS basis: use it to determine what's basic and nonbasic
  // (nonbasicFlag).
  //
  // If there's no HiGHS basis: generate nonbasicFlag, possibly by dualising and
  // performing a crash.
  //
  // Use nonbasicFlag to generate basicIndex
  //
  // Use nonbasicFlag and any HiGHS solution to determine nonbasicMove
  //
  HighsStatus return_status = HighsStatus::OK;
  const HighsOptions& options = highs_model_object.options_;
  const HighsSolution& solution = highs_model_object.solution_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsScale& scale = highs_model_object.scale_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HFactor& factor = highs_model_object.factor_;
  HMatrix& matrix = highs_model_object.matrix_;
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;
  // First determine whether the HiGHS solution space has been
  // allocated, a necessary condition for its values to be used later
  bool have_highs_solution =
      isSolutionRightSize(highs_model_object.lp_, solution);
  if (!simplex_lp_status.valid) {
    // Simplex LP is not valid so initialise the simplex LP data
    analysis.simplexTimerStart(initialiseSimplexLpDefinitionClock);
    initialiseSimplexLpDefinition(highs_model_object);
    analysis.simplexTimerStop(initialiseSimplexLpDefinitionClock);
    // Initialise the real and integer random vectors
    analysis.simplexTimerStart(initialiseSimplexLpRandomVectorsClock);
    initialiseSimplexLpRandomVectors(highs_model_object);
    analysis.simplexTimerStop(initialiseSimplexLpRandomVectorsClock);
  }
  if (simplex_lp_status.has_basis) {
    // There is a simplex basis: it should be valid - since it's set
    // internally - but check
    if (debugSimplexBasisCorrect(highs_model_object) ==
        HighsDebugStatus::LOGICAL_ERROR) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Supposed to be a Simplex basis, but incorrect");
      highs_model_object.scaled_model_status_ = HighsModelStatus::SOLVE_ERROR;
      return HighsStatus::Error;
    }
  }

  // Now we know whether the simplex basis at least has the right number
  // of basic and nonbasic variables
  if (!simplex_lp_status.has_basis) {
    // There is no simplex basis so try to identify one
    if (basis.valid_) {
      // There is is HiGHS basis: use it to construct nonbasicFlag
      if (debugBasisConsistent(options, simplex_lp, basis) ==
          HighsDebugStatus::LOGICAL_ERROR) {
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Supposed to be a Highs basis, but not valid");
        highs_model_object.scaled_model_status_ = HighsModelStatus::SOLVE_ERROR;
        return HighsStatus::Error;
      }
      analysis.simplexTimerStart(setNonbasicFlagClock);
      // Allocate memory for nonbasicFlag and set it up from the HiGHS basis
      simplex_basis.nonbasicFlag_.resize(simplex_lp.numCol_ +
                                         simplex_lp.numRow_);
      setNonbasicFlag(simplex_lp, simplex_basis.nonbasicFlag_,
                      &basis.col_status[0], &basis.row_status[0]);
      analysis.simplexTimerStop(setNonbasicFlagClock);
    }
    // nonbasicFlag is valid if the HiGHS basis exists and has the correct
    // number of basic variables
    bool nonbasicFlag_valid = basis.valid_;
    if (!nonbasicFlag_valid) {
      // So, nonbasicFlag is not valid - either because there is no
      // simplex or HiGHS basis, or because what was claimed to be
      // valid has been found to have the wrong number of basic and
      // nonbasic variables
      //
      // This is taken to imply that this is a "new" LP to be solved, so
      //
      // 1. Set simplex options from HiGHS options. This is only done with a new
      // LP so that strategy and knowledge based on run-time experience with the
      // same LP should be preserved.
      //      setSimplexOptions(highs_model_object);
      //
      // 2. Initialise the simplex timing
      //      SimplexTimer simplex_timer;
      //      simplex_timer.initialiseSimplexClocks(highs_model_object);
      //
      // 3. Generate a simplex basis, possibly by performing a crash,
      // and possibly after dualising

      /*
      // Possibly dualise, making sure that no simplex or other data are used to
      initialise
      //
      if (options.simplex_dualise_strategy != OPTION_OFF) {
      dualiseSimplexLp(highs_model_object); have_highs_solution = false;
        // Initialise the real and integer random vectors
        initialiseSimplexLpRandomVectors(highs_model_object);
      }
      */
      // Possibly permute the columns of the LP to be used by the solver.
      if (options.simplex_permute_strategy != OPTION_OFF) {
        permuteSimplexLp(highs_model_object);
      }
      analysis.simplexTimerStart(setNonbasicFlagClock);
      // Allocate memory for nonbasicFlag and set it up for a logical basis
      simplex_basis.nonbasicFlag_.resize(simplex_lp.numCol_ +
                                         simplex_lp.numRow_);
      setNonbasicFlag(simplex_lp, simplex_basis.nonbasicFlag_);
      analysis.simplexTimerStop(setNonbasicFlagClock);

      // Possibly find a crash basis
      if (options.simplex_crash_strategy != SIMPLEX_CRASH_STRATEGY_OFF) {
        HCrash crash(highs_model_object);
        analysis.simplexTimerStart(CrashClock);
        crash.crash(options.simplex_crash_strategy);
        analysis.simplexTimerStop(CrashClock);
        int num_basic_structurals = 0;
        for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
          if (simplex_basis.nonbasicFlag_[iCol] == NONBASIC_FLAG_FALSE)
            num_basic_structurals++;
        }
        HighsLogMessage(options.logfile, HighsMessageType::INFO,
                        "Crash has created a basis with %d/%d structurals",
                        num_basic_structurals, simplex_lp.numRow_);
      }
    }
    // Now that the dimensions of the LP to be solved by the simplex
    // method are known, make sure that there is a postive number of
    // rows. ToDo: Ensure that LPs with no rows can still be solved
    assert(simplex_lp.numRow_ > 0);
    if (simplex_lp.numRow_ == 0) {
      printf("Solution of LPs with no rows shouldn't reach transition()\n");
      highs_model_object.scaled_model_status_ = HighsModelStatus::SOLVE_ERROR;
      return HighsStatus::Error;
    }
    // There is now a nonbasicFlag: it should be valid - since it's
    // just been set but check
    if (debugNonbasicFlagConsistent(options, simplex_lp, simplex_basis) ==
        HighsDebugStatus::LOGICAL_ERROR) {
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "Supposed to be a Simplex basis, but nonbasicFlag not valid");
      highs_model_object.scaled_model_status_ = HighsModelStatus::SOLVE_ERROR;
      return HighsStatus::Error;
    }
    // Use nonbasicFlag to form basicIndex
    // Allocate memory for basicIndex
    analysis.simplexTimerStart(setBasicIndexClock);
    simplex_basis.basicIndex_.resize(simplex_lp.numRow_);
    int num_basic_variables = 0;
    simplex_info.num_basic_logicals = 0;
    for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE) {
        simplex_basis.basicIndex_[num_basic_variables] = iVar;
        if (iVar >= simplex_lp.numCol_) simplex_info.num_basic_logicals++;
        num_basic_variables++;
      }
    }
    analysis.simplexTimerStop(setBasicIndexClock);
    // Double-check that we have the right number of basic variables
    nonbasicFlag_valid = num_basic_variables == simplex_lp.numRow_;
    assert(nonbasicFlag_valid);
    updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_BASIS);
  }
  // Execute from here for all calls
  // Note whether a HiGHS basis can be used to (try to) choose the better bound
  // for boxed variables
  bool have_highs_basis = basis.valid_;
  //
  // Possibly scale the LP to be solved
  //
  // If the LP to be solved isn't scaled then initialise unit scaling
  // factors, to simplify things if no scaling is performed. ToDo This
  // is inefficient if the LP isn't to be scales and is repeatedly
  // hot-started - but is this really going to happen?
  if (!simplex_lp_status.scaling_tried) scaleHighsModelInit(highs_model_object);
  //
  // Scale the LP to be used by the solver if scaling is to be used and the LP
  // is not already scaled
  bool scale_lp =
      options.simplex_scale_strategy != SIMPLEX_SCALE_STRATEGY_OFF &&
      !simplex_lp_status.scaling_tried;
  const bool force_no_scaling = false;  // true;//
  if (force_no_scaling) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Forcing no scaling");
    scale_lp = false;
  }
  if (scale_lp) {
    analysis.simplexTimerStart(ScaleClock);
    scaleSimplexLp(highs_model_object);
    analysis.simplexTimerStop(ScaleClock);
#ifdef HiGHSDEV
    // Analyse the scaled LP
    if (simplex_info.analyse_lp) {
      analyseLp(simplex_lp, "Unscaled");
      HighsScale& scale = highs_model_object.scale_;
      if (scale.is_scaled_) {
        analyseVectorValues("Column scaling factors", simplex_lp.numCol_,
                            scale.col_);
        analyseVectorValues("Row    scaling factors", simplex_lp.numRow_,
                            scale.row_);
        analyseLp(simplex_lp, "Scaled");
      }
    }
#endif
  }
  // Now there is a valid nonbasicFlag and basicIndex, possibly
  // reinvert to check for basis condition/singularity
  //
  // First setup the factor arrays if they don't exist
  if (!simplex_lp_status.has_factor_arrays) {
    analysis.simplexTimerStart(factorSetupClock);
    factor.setup(simplex_lp.numCol_, simplex_lp.numRow_, &simplex_lp.Astart_[0],
                 &simplex_lp.Aindex_[0], &simplex_lp.Avalue_[0],
                 &simplex_basis.basicIndex_[0], options.highs_debug_level,
                 options.logfile, options.output, options.message_level);
    simplex_lp_status.has_factor_arrays = true;
    analysis.simplexTimerStop(factorSetupClock);
  }
  // Reinvert if there isn't a fresh INVERT. ToDo Override this for MIP hot
  // start
  //  bool reinvert = !simplex_lp_status.has_fresh_invert;
  bool reinvert = !simplex_lp_status.has_invert;
  if (reinvert) {
    analysis.simplexTimerStart(InvertClock);
    computeFactor(highs_model_object);
    analysis.simplexTimerStop(InvertClock);
    simplex_lp_status.has_fresh_invert = true;
    assert(simplex_lp_status.has_invert);
  }
  // Possibly check for basis condition. ToDo Override this for MIP hot start
  bool basis_condition_ok = true;
  if (options.simplex_initial_condition_check)
    basis_condition_ok = basisConditionOk(highs_model_object);
  // ToDo Handle ill-conditioned basis with basis crash, in which case
  // ensure that HiGHS and simplex basis are invalidated and simplex
  // work and base arrays are re-populated
  //  assert(basis_condition_ok);
  if (!basis_condition_ok) {
    // Basis crash really doesn't work, so use logical basis
    simplex_basis.basicIndex_.resize(simplex_lp.numRow_);
    for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++)
      simplex_basis.nonbasicFlag_[iCol] = NONBASIC_FLAG_TRUE;
    for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
      int iVar = simplex_lp.numCol_ + iRow;
      simplex_basis.nonbasicFlag_[iVar] = NONBASIC_FLAG_FALSE;
      simplex_basis.basicIndex_[iRow] = iVar;
    }
    simplex_info.num_basic_logicals = simplex_lp.numRow_;
    analysis.simplexTimerStart(InvertClock);
    computeFactor(highs_model_object);
    analysis.simplexTimerStop(InvertClock);

    /*
    HCrash crash(highs_model_object);
    analysis.simplexTimerStart(CrashClock);
    crash.crash(SIMPLEX_CRASH_STRATEGY_BASIC);
    analysis.simplexTimerStop(CrashClock);
     HighsLogMessage(options.logfile,
    HighsMessageType::INFO, "Performed crash to prioritise previously basic
    variables " "in well-conditioned basis");
    // Use nonbasicFlag to form basicIndex
    // Allocate memory for basicIndex
    simplex_basis.basicIndex_.resize(simplex_lp.numRow_);
    int num_basic_variables = 0;
    simplex_info.num_basic_logicals = 0;
    for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE) {
        simplex_basis.basicIndex_[num_basic_variables] = iVar;
        if (iVar >= simplex_lp.numCol_) simplex_info.num_basic_logicals++;
        num_basic_variables++;
      }
    }
    // Double-check that we have the right number of basic variables
    assert(num_basic_variables == simplex_lp.numRow_);
    updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_BASIS);
    // Report on the outcome of crash
    int num_basic_structurals =
        simplex_lp.numRow_ - simplex_info.num_basic_logicals;
    HighsLogMessage(options.logfile, HighsMessageType::INFO,
                    "Crash has created a basis with %d/%d structurals",
                    num_basic_structurals, simplex_lp.numRow_);
    // Now reinvert
    int rank_deficiency = computeFactor(highs_model_object);
    if (rank_deficiency) {
      // ToDo Handle rank deficiency by replacing singular columns with logicals
      throw runtime_error("Transition has singular basis matrix");
    }
    // Check the condition after the basis crash
    basis_condition_ok = basisConditionOk(highs_model_object);
    */
    updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_BASIS);
    simplex_lp_status.has_fresh_invert = true;
  }

  // Now there are nonbasicFlag and basicIndex corresponding to a
  // basis with well-conditioned invertible representation - and
  // corresponding nonbasicMove if simplex_lp_status.has_basis is true
  //
  // Possibly set up the HMatrix column-wise and row-wise copies of the matrix
  if (!simplex_lp_status.has_matrix_col_wise ||
      !simplex_lp_status.has_matrix_row_wise) {
    analysis.simplexTimerStart(matrixSetupClock);
    matrix.setup(simplex_lp.numCol_, simplex_lp.numRow_, &simplex_lp.Astart_[0],
                 &simplex_lp.Aindex_[0], &simplex_lp.Avalue_[0],
                 &simplex_basis.nonbasicFlag_[0]);
    simplex_lp_status.has_matrix_col_wise = true;
    simplex_lp_status.has_matrix_row_wise = true;
    analysis.simplexTimerStop(matrixSetupClock);
  }
  // If simplex basis isn't complete, set up nonbasicMove
  if (!simplex_lp_status.has_basis) {
    analysis.simplexTimerStart(setNonbasicMoveClock);
    simplex_basis.nonbasicMove_.resize(simplex_lp.numCol_ + simplex_lp.numRow_);
    setNonbasicMove(simplex_lp, scale, have_highs_basis, basis,
                    have_highs_solution, solution, simplex_basis);
    simplex_lp_status.has_basis = true;
    analysis.simplexTimerStop(setNonbasicMoveClock);
  }

  // Set up the simplex work and base arrays
  analysis.simplexTimerStart(allocateSimplexArraysClock);
  allocateWorkAndBaseArrays(highs_model_object);
  analysis.simplexTimerStop(allocateSimplexArraysClock);

  analysis.simplexTimerStart(initialiseSimplexCostBoundsClock);
  initialiseCost(highs_model_object);
  initialiseBound(highs_model_object);
  analysis.simplexTimerStop(initialiseSimplexCostBoundsClock);

  // Possibly solve for the basic primal and nonbasic dual values to determine
  // which simplex solver to use, unless it's forced
  //  if (simplex_lp_status.has_basic_primal_values) {
  initialiseNonbasicWorkValue(simplex_lp, simplex_basis, simplex_info);
  analysis.simplexTimerStart(ComputePrimalClock);
  computePrimal(highs_model_object);
  analysis.simplexTimerStop(ComputePrimalClock);
  simplex_lp_status.has_basic_primal_values = true;
  //}
  //  if (simplex_lp_status.has_basic_dual_values) {
  analysis.simplexTimerStart(ComputeDualClock);
  computeDual(highs_model_object);
  analysis.simplexTimerStop(ComputeDualClock);
  simplex_lp_status.has_nonbasic_dual_values = true;
  //}

  if (have_highs_solution) {
    // There is a HiGHS solution so possibly determine the changes in
    // basic and nonbasic values and duals for columns and rows
    if (debugSimplexHighsSolutionDifferences(highs_model_object) ==
        HighsDebugStatus::LOGICAL_ERROR)
      return HighsStatus::Error;
  }

  // Store, analyse and possibly report the number of primal and dual
  // infeasiblities and the simplex status
  computeSimplexInfeasible(highs_model_object);
  copySimplexInfeasible(highs_model_object);

  HighsSolutionParams& scaled_solution_params =
      highs_model_object.scaled_solution_params_;
  analysis.simplexTimerStart(ComputeDuObjClock);
  computeDualObjectiveValue(highs_model_object);
  analysis.simplexTimerStop(ComputeDuObjClock);
  analysis.simplexTimerStart(ComputePrObjClock);
  computePrimalObjectiveValue(highs_model_object);
  analysis.simplexTimerStop(ComputePrObjClock);
  simplex_lp_status.valid = true;

  bool primal_feasible = scaled_solution_params.num_primal_infeasibilities == 0;
  bool dual_feasible = scaled_solution_params.num_dual_infeasibilities == 0;
  if (primal_feasible && dual_feasible) {
    highs_model_object.scaled_model_status_ = HighsModelStatus::OPTIMAL;
    scaled_solution_params.primal_status =
        PrimalDualStatus::STATUS_FEASIBLE_POINT;
    scaled_solution_params.dual_status =
        PrimalDualStatus::STATUS_FEASIBLE_POINT;
  }

  scaled_solution_params.objective_function_value =
      simplex_info.primal_objective_value;

  if (debugSimplexBasicSolution("After transition", highs_model_object) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return HighsStatus::Error;

  return return_status;
}

void setNonbasicFlag(const HighsLp& simplex_lp, vector<int>& nonbasicFlag,
                     const HighsBasisStatus* col_status,
                     const HighsBasisStatus* row_status) {
  if (col_status == NULL || row_status == NULL) {
    // Initialise a logical basis
    for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
      int iVar = iCol;
      nonbasicFlag[iVar] = NONBASIC_FLAG_TRUE;
    }
    for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
      int iVar = simplex_lp.numCol_ + iRow;
      nonbasicFlag[iVar] = NONBASIC_FLAG_FALSE;
    }
  } else {
    // Initialise from HiGHS basis
    for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
      int iVar = iCol;
      if (col_status[iCol] == HighsBasisStatus::BASIC) {
        nonbasicFlag[iVar] = NONBASIC_FLAG_FALSE;
      } else {
        nonbasicFlag[iVar] = NONBASIC_FLAG_TRUE;
      }
    }
    for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
      int iVar = simplex_lp.numCol_ + iRow;
      if (row_status[iRow] == HighsBasisStatus::BASIC) {
        nonbasicFlag[iVar] = NONBASIC_FLAG_FALSE;
      } else {
        nonbasicFlag[iVar] = NONBASIC_FLAG_TRUE;
      }
    }
  }
}

void setNonbasicMove(const HighsLp& simplex_lp, const HighsScale& scale,
                     const bool have_highs_basis, const HighsBasis& basis,
                     const bool have_highs_solution,
                     const HighsSolution& solution,
                     SimplexBasis& simplex_basis) {
  // Don't have a simplex basis since nonbasicMove is not set up.
  const int illegal_move_value = -99;

  // Assign nonbasicMove using as much information as is available
  double lower;
  double upper;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;

  for (int iVar = 0; iVar < numTot; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) {
      // Basic variable
      simplex_basis.nonbasicMove_[iVar] = 0;
      continue;
    }
    // Nonbasic variable
    if (iVar < simplex_lp.numCol_) {
      lower = simplex_lp.colLower_[iVar];
      upper = simplex_lp.colUpper_[iVar];
    } else {
      int iRow = iVar - simplex_lp.numCol_;
      lower = -simplex_lp.rowUpper_[iRow];
      upper = -simplex_lp.rowLower_[iRow];
    }
    int move = illegal_move_value;
    if (lower == upper) {
      // Fixed
      move = NONBASIC_MOVE_ZE;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        //
        // Determine the bound to set the value to according to, in order of
        // priority
        //
        // 1. Any valid HiGHS basis status
        if (have_highs_basis) {
          if (iVar < simplex_lp.numCol_) {
            if (basis.col_status[iVar] == HighsBasisStatus::LOWER) {
              move = NONBASIC_MOVE_UP;
            } else if (basis.col_status[iVar] == HighsBasisStatus::UPPER) {
              move = NONBASIC_MOVE_DN;
            }
          } else {
            int iRow = iVar - simplex_lp.numCol_;
            if (basis.row_status[iRow] == HighsBasisStatus::LOWER) {
              move = NONBASIC_MOVE_DN;
            } else if (basis.row_status[iRow] == HighsBasisStatus::UPPER) {
              move = NONBASIC_MOVE_UP;
            }
          }
        }
        // 2. Any HiGHS solution value
        if (move == illegal_move_value && have_highs_solution) {
          // Reach here if there is no HiGHS basis or the HiGHS
          // nonbasic status is just NONBASIC.
          double midpoint = 0.5 * (lower + upper);
          double value_from_highs_solution;
          if (iVar < simplex_lp.numCol_) {
            assert(!have_highs_basis ||
                   basis.col_status[iVar] == HighsBasisStatus::NONBASIC);
            value_from_highs_solution =
                solution.col_value[iVar] / scale.col_[iVar];
          } else {
            int iRow = iVar - simplex_lp.numCol_;
            assert(!have_highs_basis ||
                   basis.row_status[iRow] == HighsBasisStatus::NONBASIC);
            value_from_highs_solution =
                -solution.row_value[iRow] * scale.row_[iRow];
          }
          if (value_from_highs_solution < midpoint) {
            move = NONBASIC_MOVE_UP;
          } else {
            move = NONBASIC_MOVE_DN;
          }
        }
        // 3. Bound of original LP that is closer to zero
        if (move == illegal_move_value) {
          if (fabs(lower) < fabs(upper)) {
            move = NONBASIC_MOVE_UP;
          } else {
            move = NONBASIC_MOVE_DN;
          }
        }
      } else {
        // Lower (since upper bound is infinite)
        move = NONBASIC_MOVE_UP;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = NONBASIC_MOVE_DN;
    } else {
      // FREE
      move = NONBASIC_MOVE_ZE;
    }
    assert(move != illegal_move_value);
    simplex_basis.nonbasicMove_[iVar] = move;
  }
}

void initialiseNonbasicWorkValue(const HighsLp& simplex_lp,
                                 const SimplexBasis& simplex_basis,
                                 HighsSimplexInfo& simplex_info) {
  // Assign nonbasic values from bounds and (if necessary) nonbasicMove
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int iVar = 0; iVar < numTot; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic variable
    const double lower = simplex_info.workLower_[iVar];
    const double upper = simplex_info.workUpper_[iVar];
    double value;
    if (lower == upper) {
      value = lower;
    } else if (simplex_basis.nonbasicMove_[iVar] == NONBASIC_MOVE_UP) {
      value = lower;
    } else if (simplex_basis.nonbasicMove_[iVar] == NONBASIC_MOVE_DN) {
      value = upper;
    } else {
      assert(simplex_basis.nonbasicMove_[iVar] == NONBASIC_MOVE_ZE);
      value = 0;
    }
    simplex_info.workValue_[iVar] = value;
  }
}

bool basisConditionOk(HighsModelObject& highs_model_object) {
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;
  bool basis_condition_ok;
  analysis.simplexTimerStart(BasisConditionClock);
  double basis_condition = computeBasisCondition(highs_model_object);
  analysis.simplexTimerStop(BasisConditionClock);
  double basis_condition_tolerance =
      highs_model_object.options_.simplex_initial_condition_tolerance;
  basis_condition_ok = basis_condition < basis_condition_tolerance;
  HighsMessageType message_type = HighsMessageType::INFO;
  std::string condition_comment;
  if (basis_condition_ok) {
    condition_comment = "is within";
  } else {
    message_type = HighsMessageType::WARNING;
    condition_comment = "exceeds";
  }
  HighsLogMessage(
      highs_model_object.options_.logfile, message_type,
      "Initial basis condition estimate of %11.4g %s the tolerance of %g",
      basis_condition, condition_comment.c_str(), basis_condition_tolerance);
  return basis_condition_ok;
}
bool dual_infeasible(const double value, const double lower, const double upper,
                     const double dual, const double value_tolerance,
                     const double dual_tolerance) {
  double midpoint = (lower + upper) * 0.5;
  double residual = max(lower - value, value - upper);
  bool infeasible = false;
  if (highs_isInfinity(-lower)) {
    // Infinite lower bound
    if (highs_isInfinity(upper)) {
      // Infinite upper bound
      // Free
      infeasible = fabs(dual) >= dual_tolerance;
    } else {
      // Finite upper bound
      // Upper bounded - and assumed to be nonbasic at that bound
      if (fabs(residual) >= value_tolerance) {
        printf("dual_infeasible: %12g %12g %12g %12g %12g\n", value, lower,
               upper, residual, value_tolerance);
      }
      assert(fabs(residual) < value_tolerance);
      infeasible = dual >= dual_tolerance;
    }
  } else {
    // Finite lower bound
    if (highs_isInfinity(upper)) {
      // Infinite upper bound
      // Lower bounded - and assumed to be nonbasic at that bound
      assert(fabs(residual) < value_tolerance);
      infeasible = dual <= -dual_tolerance;
    } else {
      // Finite upper bound
      // Assumed to be nonbasic at that bound
      assert(fabs(residual) < value_tolerance);
      if (lower < upper) {
        // Boxed
        if (value < midpoint) {
          // At lower bound
          infeasible = dual <= -dual_tolerance;
        } else {
          // At upper bound
          infeasible = dual >= dual_tolerance;
        }
      } else {
        // Fixed
        infeasible = false;
      }
    }
  }
  return infeasible;
}

void appendNonbasicColsToBasis(HighsLp& lp, HighsBasis& basis, int XnumNewCol) {
  assert(basis.valid_);
  if (!basis.valid_) {
    printf("\n!!Appending columns to invalid basis!!\n\n");
  }
  // Add nonbasic structurals
  if (XnumNewCol == 0) return;
  int newNumCol = lp.numCol_ + XnumNewCol;
  basis.col_status.resize(newNumCol);
  // Make any new columns nonbasic
  for (int iCol = lp.numCol_; iCol < newNumCol; iCol++) {
    if (!highs_isInfinity(-lp.colLower_[iCol])) {
      // Has finite lower bound so set it there
      basis.col_status[iCol] = HighsBasisStatus::LOWER;
    } else if (!highs_isInfinity(lp.colUpper_[iCol])) {
      // Has finite upper bound so set it there
      basis.col_status[iCol] = HighsBasisStatus::UPPER;
    } else {
      // Free variable so set to zero
      basis.col_status[iCol] = HighsBasisStatus::ZERO;
    }
  }
}

void appendNonbasicColsToBasis(HighsLp& lp, SimplexBasis& basis,
                               int XnumNewCol) {
  // Add nonbasic structurals
  if (XnumNewCol == 0) return;
  int newNumCol = lp.numCol_ + XnumNewCol;
  int newNumTot = newNumCol + lp.numRow_;
  basis.nonbasicFlag_.resize(newNumTot);
  basis.nonbasicMove_.resize(newNumTot);
  // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
  // necessary
  for (int iRow = lp.numRow_ - 1; iRow >= 0; iRow--) {
    int iCol = basis.basicIndex_[iRow];
    if (iCol >= lp.numCol_) {
      // This basic variable is a row, so shift its index
      basis.basicIndex_[iRow] += XnumNewCol;
    }
    basis.nonbasicFlag_[newNumCol + iRow] =
        basis.nonbasicFlag_[lp.numCol_ + iRow];
    basis.nonbasicMove_[newNumCol + iRow] =
        basis.nonbasicMove_[lp.numCol_ + iRow];
  }
  // Make any new columns nonbasic
  const int illegal_move_value = -99;
  for (int iCol = lp.numCol_; iCol < newNumCol; iCol++) {
    basis.nonbasicFlag_[iCol] = NONBASIC_FLAG_TRUE;
    double lower = lp.colLower_[iCol];
    double upper = lp.colUpper_[iCol];
    int move = illegal_move_value;
    if (lower == upper) {
      // Fixed
      move = NONBASIC_MOVE_ZE;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (fabs(lower) < fabs(upper)) {
          move = NONBASIC_MOVE_UP;
        } else {
          move = NONBASIC_MOVE_DN;
        }
      } else {
        // Lower (since upper bound is infinite)
        move = NONBASIC_MOVE_UP;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = NONBASIC_MOVE_DN;
    } else {
      // FREE
      move = NONBASIC_MOVE_ZE;
    }
    assert(move != illegal_move_value);
    basis.nonbasicMove_[iCol] = move;
  }
}

void appendBasicRowsToBasis(HighsLp& lp, HighsBasis& basis, int XnumNewRow) {
  assert(basis.valid_);
  if (!basis.valid_) {
    printf("\n!!Appending columns to invalid basis!!\n\n");
  }
  // Add basic logicals
  if (XnumNewRow == 0) return;
  int newNumRow = lp.numRow_ + XnumNewRow;
  basis.row_status.resize(newNumRow);
  // Make the new rows basic
  for (int iRow = lp.numRow_; iRow < newNumRow; iRow++) {
    basis.row_status[iRow] = HighsBasisStatus::BASIC;
  }
}

void appendBasicRowsToBasis(HighsLp& lp, SimplexBasis& basis, int XnumNewRow) {
  // Add basic logicals
  if (XnumNewRow == 0) return;

  int newNumRow = lp.numRow_ + XnumNewRow;
  int newNumTot = lp.numCol_ + newNumRow;
  basis.nonbasicFlag_.resize(newNumTot);
  basis.nonbasicMove_.resize(newNumTot);
  basis.basicIndex_.resize(newNumRow);
  // Make the new rows basic
  for (int iRow = lp.numRow_; iRow < newNumRow; iRow++) {
    basis.nonbasicFlag_[lp.numCol_ + iRow] = NONBASIC_FLAG_FALSE;
    basis.nonbasicMove_[lp.numCol_ + iRow] = 0;
    basis.basicIndex_[iRow] = lp.numCol_ + iRow;
  }
}

void reportBasis(const HighsOptions options, const HighsLp& lp,
                 const HighsBasis& basis) {
  if (lp.numCol_ > 0)
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "HighsBasis\n   Col Status\n");
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "%6d %6d\n", iCol, (int)basis.col_status[iCol]);
  }
  if (lp.numRow_ > 0)
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "   Row Status\n");
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "%6d %6d\n", iRow, (int)basis.row_status[iRow]);
  }
}

void reportBasis(const HighsOptions options, const HighsLp& lp,
                 const SimplexBasis& simplex_basis) {
  if (lp.numCol_ > 0)
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "SimplexBasis\n   Var    Col   Flag\n");
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    int iVar = iCol;
    if (simplex_basis.nonbasicFlag_[iVar])
      HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                        "%6d %6d %6d\n", iVar, iCol,
                        simplex_basis.nonbasicFlag_[iVar]);
    else
      HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                        "%6d %6d %6d\n", iVar, iCol,
                        simplex_basis.nonbasicFlag_[iVar]);
  }
  if (lp.numRow_ > 0)
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "   Var    Row   Flag  Basic\n");
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    int iVar = lp.numCol_ + iRow;
    if (simplex_basis.nonbasicFlag_[iVar])
      HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                        "%6d %6d %6d %6d\n", iVar, iRow,
                        simplex_basis.nonbasicFlag_[iVar],
                        simplex_basis.basicIndex_[iRow]);
    else
      HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                        "%6d %6d %6d %6d\n", iVar, iRow,
                        simplex_basis.nonbasicFlag_[iVar],
                        simplex_basis.basicIndex_[iRow]);
  }
}

/**
 * @brief Simplex utilities
 */

void computeDualObjectiveValue(HighsModelObject& highs_model_object,
                               int phase) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;

  simplex_info.dual_objective_value = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++) {
    if (highs_model_object.simplex_basis_.nonbasicFlag_[i]) {
      const double term =
          simplex_info.workValue_[i] * simplex_info.workDual_[i];
      if (term) {
        simplex_info.dual_objective_value +=
            simplex_info.workValue_[i] * simplex_info.workDual_[i];
      }
    }
  }
  simplex_info.dual_objective_value *= highs_model_object.scale_.cost_;
  if (phase != 1) {
    // In phase 1 the dual objective has no objective
    // shift. Otherwise, if minimizing the shift is added. If
    // maximizing, workCost (and hence workDual) are negated, so the
    // shift is subtracted. Hence the shift is added according to the
    // sign implied by sense_
    simplex_info.dual_objective_value +=
        ((int)simplex_lp.sense_) * simplex_lp.offset_;
  }
  // Now have dual objective value
  simplex_lp_status.has_dual_objective_value = true;
}

int setSourceOutFmBd(const HighsModelObject& highs_model_object,
                     const int columnOut) {
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  int sourceOut = 0;
  if (simplex_info.workLower_[columnOut] !=
      simplex_info.workUpper_[columnOut]) {
    if (!highs_isInfinity(-simplex_info.workLower_[columnOut])) {
      // Finite LB so sourceOut = -1 ensures value set to LB if LB < UB
      sourceOut = -1;
      //      printf("STRANGE: variable %d leaving the basis is [%11.4g, %11.4g]
      //      so setting sourceOut = -1\n", columnOut,
      //      simplex_info.workLower_[columnOut],
      //      simplex_info.workUpper_[columnOut]);
    } else {
      // Infinite LB so sourceOut = 1 ensures value set to UB
      sourceOut = 1;
      if (!highs_isInfinity(simplex_info.workUpper_[columnOut])) {
        // Free variable => trouble!
        printf("TROUBLE: variable %d leaving the basis is free!\n", columnOut);
      }
    }
  }
  return sourceOut;
}
void computePrimalObjectiveValue(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  simplex_info.primal_objective_value = 0;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    int iVar = simplex_basis.basicIndex_[iRow];
    if (iVar < simplex_lp.numCol_) {
      simplex_info.primal_objective_value +=
          simplex_info.baseValue_[iRow] * simplex_lp.colCost_[iVar];
    }
  }
  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    if (simplex_basis.nonbasicFlag_[iCol])
      simplex_info.primal_objective_value +=
          simplex_info.workValue_[iCol] * simplex_lp.colCost_[iCol];
  }
  simplex_info.primal_objective_value *= highs_model_object.scale_.cost_;
  // Objective value calculation is done using primal values and
  // original costs so offset is vanilla
  simplex_info.primal_objective_value += simplex_lp.offset_;
  // Now have primal objective value
  simplex_lp_status.has_primal_objective_value = true;
}

#ifdef HiGHSDEV
void getPrimalValue(const HighsModelObject& highs_model_object,
                    vector<double>& primal_value) {
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  // Copy all of workValue to get all the nonbasic values
  primal_value.resize(simplex_lp.numCol_ + simplex_lp.numRow_);
  for (int iCol = 0; iCol < simplex_lp.numCol_ + simplex_lp.numRow_; iCol++)
    primal_value[iCol] = simplex_info.workValue_[iCol];
  // Over-write the value of the nonbasic variables
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
    primal_value[simplex_basis.basicIndex_[iRow]] =
        simplex_info.baseValue_[iRow];
}

void analysePrimalObjectiveValue(const HighsModelObject& highs_model_object) {
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  HighsValueDistribution objective_value_term_distribution;
  HighsValueDistribution basic_value_distribution;
  HighsValueDistribution basic_cost_distribution;
  initialiseValueDistribution("Nonzero objective terms", "", 1e-16, 1e16, 10.0,
                              objective_value_term_distribution);
  initialiseValueDistribution("Basic values", "", 1e-16, 1e16, 10.0,
                              basic_value_distribution);
  initialiseValueDistribution("Nonzero basic costs", "", 1e-16, 1e16, 10.0,
                              basic_cost_distribution);

  double primal_objective_value = 0;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    int iVar = simplex_basis.basicIndex_[iRow];
    const double value = simplex_info.baseValue_[iRow];
    updateValueDistribution(value, basic_value_distribution);
    if (iVar < simplex_lp.numCol_) {
      const double cost = simplex_lp.colCost_[iVar];
      if (cost) {
        updateValueDistribution(cost, basic_cost_distribution);
        const double term = value * cost;
        primal_objective_value += term;
        const double abs_term = fabs(term);
        updateValueDistribution(abs_term, objective_value_term_distribution);
      }
    }
  }
  HighsValueDistribution nonbasic_value_distribution;
  HighsValueDistribution nonbasic_cost_distribution;
  initialiseValueDistribution("Nonbasic values", "", 1e-16, 1e16, 10.0,
                              nonbasic_value_distribution);
  initialiseValueDistribution("Nonzero nonbasic costs", "", 1e-16, 1e16, 10.0,
                              nonbasic_cost_distribution);
  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    if (simplex_basis.nonbasicFlag_[iCol]) {
      const double value = simplex_info.workValue_[iCol];
      updateValueDistribution(value, nonbasic_value_distribution);
      const double cost = simplex_lp.colCost_[iCol];
      if (cost) {
        updateValueDistribution(cost, nonbasic_cost_distribution);
        const double term = value * cost;
        primal_objective_value += term;
        const double abs_term = fabs(term);
        updateValueDistribution(abs_term, objective_value_term_distribution);
      }
    }
  }
  for (int iCol = simplex_lp.numCol_;
       iCol < simplex_lp.numCol_ + simplex_lp.numRow_; iCol++) {
    if (simplex_basis.nonbasicFlag_[iCol]) {
      const double value = simplex_info.workValue_[iCol];
      updateValueDistribution(value, nonbasic_value_distribution);
    }
  }
  printf("\nAnalysis of values, costs and objective terms:\n");
  printValueDistribution(nonbasic_value_distribution);
  printValueDistribution(basic_value_distribution);
  printValueDistribution(nonbasic_cost_distribution);
  printValueDistribution(basic_cost_distribution);
  printValueDistribution(objective_value_term_distribution);
  printf("Linear objective value: %g\n", primal_objective_value);
  primal_objective_value *= highs_model_object.scale_.cost_;
  printf("Scaled objective value: %g\n", primal_objective_value);
  // Objective value calculation is done using primal values and
  // original costs so offset is vanilla
  primal_objective_value -= simplex_lp.offset_;
  printf("Offset objective value: %g\n", primal_objective_value);
}
#endif

void initialiseSimplexLpDefinition(HighsModelObject& highs_model_object) {
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  // Ensure that the simplex LP is fully invalidated
  invalidateSimplexLp(simplex_lp_status);
  // Copy the LP to the structure to be used by the solver
  highs_model_object.simplex_lp_ = highs_model_object.lp_;
}

void initialiseSimplexLpRandomVectors(HighsModelObject& highs_model_object) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const int numCol = highs_model_object.simplex_lp_.numCol_;
  const int numTot = highs_model_object.simplex_lp_.numCol_ +
                     highs_model_object.simplex_lp_.numRow_;
  // Instantiate and (re-)initialise the random number generator
  HighsRandom& random = highs_model_object.random_;
  random.initialise();
  //
  // Generate a random permutation of the column indices
  simplex_info.numColPermutation_.resize(numCol);
  int* numColPermutation = &simplex_info.numColPermutation_[0];
  for (int i = 0; i < numCol; i++) numColPermutation[i] = i;
  for (int i = numCol - 1; i >= 1; i--) {
    int j = random.integer() % (i + 1);
    std::swap(numColPermutation[i], numColPermutation[j]);
  }

  // Re-initialise the random number generator and generate the
  // random vectors in the same order as hsol to maintain repeatable
  // performance
  random.initialise();
  //
  // Generate a random permutation of all the indices
  simplex_info.numTotPermutation_.resize(numTot);
  int* numTotPermutation = &simplex_info.numTotPermutation_[0];
  for (int i = 0; i < numTot; i++) numTotPermutation[i] = i;
  for (int i = numTot - 1; i >= 1; i--) {
    int j = random.integer() % (i + 1);
    std::swap(numTotPermutation[i], numTotPermutation[j]);
  }

  // Generate a vector of random reals
  simplex_info.numTotRandomValue_.resize(numTot);
  double* numTotRandomValue = &simplex_info.numTotRandomValue_[0];
  for (int i = 0; i < numTot; i++) {
    numTotRandomValue[i] = random.fraction();
  }
}

// SCALING:
#ifdef HiGHSDEV
// Information on large costs
const double tlLargeCo = 1e5;
int numLargeCo;
vector<int> largeCostFlag;
double largeCostScale;
#endif

void scaleHighsModelInit(HighsModelObject& highs_model_object) {
  HighsScale& scale = highs_model_object.scale_;
  scale.is_scaled_ = false;
  scale.col_.assign(highs_model_object.simplex_lp_.numCol_, 1);
  scale.row_.assign(highs_model_object.simplex_lp_.numRow_, 1);
  scale.cost_ = 1;
#ifdef HiGHSDEV
  //  largeCostScale = 1;
#endif
}

void scaleCosts(HighsModelObject& highs_model_object) {
  // Scale the costs by no less than minAlwCostScale
  double max_allowed_cost_scale =
      pow(2.0, highs_model_object.options_.allowed_simplex_cost_scale_factor);
  double cost_scale;
  double max_nonzero_cost = 0;
  for (int iCol = 0; iCol < highs_model_object.simplex_lp_.numCol_; iCol++) {
    if (highs_model_object.simplex_lp_.colCost_[iCol]) {
      max_nonzero_cost =
          max(fabs(highs_model_object.simplex_lp_.colCost_[iCol]),
              max_nonzero_cost);
    }
  }
  // Scaling the costs up effectively increases the dual tolerance to
  // which the problem is solved - so, if the max cost is small the
  // scaling factor pushes it up by a power of 2 so it's close to 1
  // Scaling the costs down effectively decreases the dual tolerance
  // to which the problem is solved - so this can't be done too much
  cost_scale = 1;
  const double ln2 = log(2.0);
  // Scale if the max cost is positive and outside the range [1/16, 16]
  if ((max_nonzero_cost > 0) &&
      ((max_nonzero_cost < (1.0 / 16)) || (max_nonzero_cost > 16))) {
    cost_scale = max_nonzero_cost;
    cost_scale = pow(2.0, floor(log(cost_scale) / ln2 + 0.5));
    cost_scale = min(cost_scale, max_allowed_cost_scale);
  }
  highs_model_object.scale_.cost_ = cost_scale;
  if (cost_scale == 1) return;
  // Scale the costs (and record of max_nonzero_cost) by cost_scale, being at
  // most max_allowed_cost_scale
  for (int iCol = 0; iCol < highs_model_object.simplex_lp_.numCol_; iCol++) {
    highs_model_object.simplex_lp_.colCost_[iCol] /= cost_scale;
  }
  max_nonzero_cost /= cost_scale;

#ifdef HiGHSDEV
  /*
    bool alwLargeCostScaling = false;
    if (alwLargeCostScaling && (numLargeCo > 0)) {
    // Scale any large costs by largeCostScale, being at most (a further)
    // max_allowed_cost_scale
    largeCostScale = max_nonzero_cost;
    largeCostScale = pow(2.0, floor(log(largeCostScale) / ln2 + 0.5));
    largeCostScale = min(largeCostScale, max_allowed_cost_scale);
    printf(
    "   Scaling all |cost| > %11.4g by %11.4g\ngrep_LargeCostScale,%g,%g\n",
    tlLargeCo, largeCostScale, tlLargeCo, largeCostScale);
    for (int iCol = 0; iCol < highs_model_object.simplex_lp_.numCol_; iCol++) {
    if (largeCostFlag[iCol]) {
    highs_model_object.simplex_lp_.colCost_[iCol] /= largeCostScale;
    }
    }
    }
  */
  //  utils.analyseVectorValues("Column costs",
  //  highs_model_object.simplex_lp_.numCol_,
  //  highs_model_object.simplex_lp_.colCost_);
#endif
}

void scaleFactorRanges(HighsModelObject& highs_model_object,
                       double& min_col_scale, double& max_col_scale,
                       double& min_row_scale, double& max_row_scale) {
  int numCol = highs_model_object.simplex_lp_.numCol_;
  int numRow = highs_model_object.simplex_lp_.numRow_;
  double* colScale = &highs_model_object.scale_.col_[0];
  double* rowScale = &highs_model_object.scale_.row_[0];
  // Determine the max and min row and column scaling factors
  min_col_scale = HIGHS_CONST_INF;
  max_col_scale = 1 / HIGHS_CONST_INF;
  min_row_scale = HIGHS_CONST_INF;
  max_row_scale = 1 / HIGHS_CONST_INF;
  for (int iCol = 0; iCol < numCol; iCol++) {
    min_col_scale = min(colScale[iCol], min_col_scale);
    max_col_scale = max(colScale[iCol], max_col_scale);
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
    min_row_scale = min(rowScale[iRow], min_row_scale);
    max_row_scale = max(rowScale[iRow], max_row_scale);
  }
}

void scaleSimplexLp(HighsModelObject& highs_model_object) {
  //  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.scaling_tried) return;
  // Scale the LP highs_model_object.simplex_lp_, assuming all data are in place
  HighsScale& scale = highs_model_object.scale_;
  // Reset all scaling to 1
  scaleHighsModelInit(highs_model_object);
  int numCol = highs_model_object.simplex_lp_.numCol_;
  int numRow = highs_model_object.simplex_lp_.numRow_;
  double* colScale = &highs_model_object.scale_.col_[0];
  double* rowScale = &highs_model_object.scale_.row_[0];
  int* Astart = &highs_model_object.simplex_lp_.Astart_[0];
  double* Avalue = &highs_model_object.simplex_lp_.Avalue_[0];
  double* colCost = &highs_model_object.simplex_lp_.colCost_[0];
  double* colLower = &highs_model_object.simplex_lp_.colLower_[0];
  double* colUpper = &highs_model_object.simplex_lp_.colUpper_[0];
  double* rowLower = &highs_model_object.simplex_lp_.rowLower_[0];
  double* rowUpper = &highs_model_object.simplex_lp_.rowUpper_[0];

  // Allow a switch to/from the original scaling rules
  int simplex_scale_strategy =
      highs_model_object.options_.simplex_scale_strategy;
  bool allow_cost_scaling =
      highs_model_object.options_.allowed_simplex_cost_scale_factor > 0;
  // Find out range of matrix values and skip matrix scaling if all
  // |values| are in [0.2, 5]
  const double no_scaling_original_matrix_min_value = 0.2;
  const double no_scaling_original_matrix_max_value = 5.0;
  double original_matrix_min_value = HIGHS_CONST_INF;
  double original_matrix_max_value = 0;
  for (int k = 0, AnX = Astart[numCol]; k < AnX; k++) {
    double value = fabs(Avalue[k]);
    original_matrix_min_value = min(original_matrix_min_value, value);
    original_matrix_max_value = max(original_matrix_max_value, value);
  }
  bool no_scaling =
      (original_matrix_min_value >= no_scaling_original_matrix_min_value) &&
      (original_matrix_max_value <= no_scaling_original_matrix_max_value);
  const bool force_scaling = false;
  if (force_scaling) {
    no_scaling = false;
    printf("!!!! FORCE SCALING !!!!\n");
  }
  bool scaled_matrix = false;
  if (no_scaling) {
    // No matrix scaling, but possible cost scaling
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Scaling: Matrix has [min, max] values of [%g, %g] within "
                    "[%g, %g] so no scaling performed",
                    original_matrix_min_value, original_matrix_max_value,
                    no_scaling_original_matrix_min_value,
                    no_scaling_original_matrix_max_value);
  } else {
    const bool equilibration_scaling =
        simplex_scale_strategy == SIMPLEX_SCALE_STRATEGY_HIGHS ||
        simplex_scale_strategy == SIMPLEX_SCALE_STRATEGY_HIGHS_FORCED;
    if (equilibration_scaling) {
      scaled_matrix = equilibrationScaleMatrix(highs_model_object);
    } else {
      scaled_matrix = maxValueScaleMatrix(highs_model_object);
    }
    scale.is_scaled_ = scaled_matrix;
    if (scaled_matrix) {
      // Matrix is scaled, so scale the bounds and costs
      for (int iCol = 0; iCol < numCol; iCol++) {
        colLower[iCol] /= colScale[iCol];
        colUpper[iCol] /= colScale[iCol];
        colCost[iCol] *= colScale[iCol];
      }
      for (int iRow = 0; iRow < numRow; iRow++) {
        rowLower[iRow] *= rowScale[iRow];
        rowUpper[iRow] *= rowScale[iRow];
      }
    }
  }
  // Possibly scale the costs
  if (allow_cost_scaling) scaleCosts(highs_model_object);

  // If matrix is unscaled, then LP is only scaled if there is a cost scaling
  // factor
  if (!scaled_matrix) scale.is_scaled_ = scale.cost_ != 1;

  // Deduce the consequences of scaling the LP
  if (scale.is_scaled_)
    updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                          LpAction::SCALE);
}

bool equilibrationScaleMatrix(HighsModelObject& highs_model_object) {
  int numCol = highs_model_object.simplex_lp_.numCol_;
  int numRow = highs_model_object.simplex_lp_.numRow_;
  double* colScale = &highs_model_object.scale_.col_[0];
  double* rowScale = &highs_model_object.scale_.row_[0];
  int* Astart = &highs_model_object.simplex_lp_.Astart_[0];
  int* Aindex = &highs_model_object.simplex_lp_.Aindex_[0];
  double* Avalue = &highs_model_object.simplex_lp_.Avalue_[0];
  double* colCost = &highs_model_object.simplex_lp_.colCost_[0];

  int simplex_scale_strategy =
      highs_model_object.options_.simplex_scale_strategy;

  double original_matrix_min_value = HIGHS_CONST_INF;
  double original_matrix_max_value = 0;
  for (int k = 0, AnX = Astart[numCol]; k < AnX; k++) {
    double value = fabs(Avalue[k]);
    original_matrix_min_value = min(original_matrix_min_value, value);
    original_matrix_max_value = max(original_matrix_max_value, value);
  }

  // Include cost in scaling if minimum nonzero cost is less than 0.1
  double min_nonzero_cost = HIGHS_CONST_INF;
  for (int i = 0; i < numCol; i++) {
    if (colCost[i]) min_nonzero_cost = min(fabs(colCost[i]), min_nonzero_cost);
  }
  bool include_cost_in_scaling = false;
  include_cost_in_scaling = min_nonzero_cost < 0.1;

  // Limits on scaling factors
  double max_allow_scale;
  double min_allow_scale;
  // Now that HIGHS_CONST_INF =
  // std::numeric_limits<double>::infinity(), this Qi-trick doesn't
  // work so, in recognition, use the old value of HIGHS_CONST_INF
  const double finite_infinity = 1e200;
  max_allow_scale =
      pow(2.0, highs_model_object.options_.allowed_simplex_matrix_scale_factor);
  min_allow_scale = 1 / max_allow_scale;

  double min_allow_col_scale = min_allow_scale;
  double max_allow_col_scale = max_allow_scale;
  double min_allow_row_scale = min_allow_scale;
  double max_allow_row_scale = max_allow_scale;

  // Search up to 6 times
  vector<double> row_min_value(numRow, finite_infinity);
  vector<double> row_max_value(numRow, 1 / finite_infinity);
  for (int search_count = 0; search_count < 6; search_count++) {
    // Find column scale, prepare row data
    for (int iCol = 0; iCol < numCol; iCol++) {
      // For column scale (find)
      double col_min_value = finite_infinity;
      double col_max_value = 1 / finite_infinity;
      double abs_col_cost = fabs(colCost[iCol]);
      if (include_cost_in_scaling && abs_col_cost != 0) {
        col_min_value = min(col_min_value, abs_col_cost);
        col_max_value = max(col_max_value, abs_col_cost);
      }
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        double value = fabs(Avalue[k]) * rowScale[Aindex[k]];
        col_min_value = min(col_min_value, value);
        col_max_value = max(col_max_value, value);
      }
      double col_equilibration = 1 / sqrt(col_min_value * col_max_value);
      // Ensure that column scale factor is not excessively large or small
      colScale[iCol] =
          min(max(min_allow_col_scale, col_equilibration), max_allow_col_scale);
      // For row scale (only collect)
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        int iRow = Aindex[k];
        double value = fabs(Avalue[k]) * colScale[iCol];
        row_min_value[iRow] = min(row_min_value[iRow], value);
        row_max_value[iRow] = max(row_max_value[iRow], value);
      }
    }
    // For row scale (find)
    for (int iRow = 0; iRow < numRow; iRow++) {
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
  for (int iCol = 0; iCol < numCol; iCol++) {
    colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / log2 + 0.5));
    min_col_scale = min(colScale[iCol], min_col_scale);
    max_col_scale = max(colScale[iCol], max_col_scale);
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
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
  for (int iCol = 0; iCol < numCol; iCol++) {
    double original_col_min_value = finite_infinity;
    double original_col_max_value = 1 / finite_infinity;
    double col_min_value = finite_infinity;
    double col_max_value = 1 / finite_infinity;
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      int iRow = Aindex[k];
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

  for (int iRow = 0; iRow < numRow; iRow++) {
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
#ifdef HiGHSDEV
  HighsLogMessage(
      highs_model_object.options_.logfile, HighsMessageType::INFO,
      "Scaling: Original equilibration: min/mean/max %11.4g/%11.4g/%11.4g "
      "(cols); min/mean/max %11.4g/%11.4g/%11.4g (rows)",
      min_original_col_equilibration, geomean_original_col_equilibration,
      max_original_col_equilibration, min_original_row_equilibration,
      geomean_original_row_equilibration, max_original_row_equilibration);
  HighsLogMessage(
      highs_model_object.options_.logfile, HighsMessageType::INFO,
      "Scaling: Final    equilibration: min/mean/max %11.4g/%11.4g/%11.4g "
      "(cols); min/mean/max %11.4g/%11.4g/%11.4g (rows)",
      min_col_equilibration, geomean_col_equilibration, max_col_equilibration,
      min_row_equilibration, geomean_row_equilibration, max_row_equilibration);
#endif

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
      (geomean_original_col * geomean_original_row) /
      (geomean_col * geomean_row);
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
#ifdef HiGHSDEV
  HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                  "Scaling: Extreme equilibration improvement = ( %11.4g + "
                  "%11.4g) / ( %11.4g + %11.4g) = %11.4g / %11.4g = %11.4g",
                  original_col_ratio, original_row_ratio, col_ratio, row_ratio,
                  (original_col_ratio + original_row_ratio),
                  (col_ratio + row_ratio), extreme_equilibration_improvement);
  HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                  "Scaling:    Mean equilibration improvement = ( %11.4g * "
                  "%11.4g) / ( %11.4g * %11.4g) = %11.4g / %11.4g = %11.4g",
                  geomean_original_col, geomean_original_row, geomean_col,
                  geomean_row, (geomean_original_col * geomean_original_row),
                  (geomean_col * geomean_row), mean_equilibration_improvement);
  HighsLogMessage(
      highs_model_object.options_.logfile, HighsMessageType::INFO,
      "Scaling: Yields [min, max, ratio] matrix values of [%0.4g, %0.4g, "
      "%0.4g]; Originally [%0.4g, %0.4g, %0.4g]: Improvement of %0.4g",
      matrix_min_value, matrix_max_value, matrix_value_ratio,
      original_matrix_min_value, original_matrix_max_value,
      original_matrix_value_ratio, matrix_value_ratio_improvement);
  HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                  "Scaling: Improves    mean equilibration by a factor %0.4g",
                  mean_equilibration_improvement);
  HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                  "Scaling: Improves extreme equilibration by a factor %0.4g",
                  extreme_equilibration_improvement);
  HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                  "Scaling: Improves max/min matrix values by a factor %0.4g",
                  matrix_value_ratio_improvement);
#endif
  const bool possibly_abandon_scaling =
      simplex_scale_strategy != SIMPLEX_SCALE_STRATEGY_HIGHS_FORCED;
  const double improvement_factor = extreme_equilibration_improvement *
                                    mean_equilibration_improvement *
                                    matrix_value_ratio_improvement;

  const double improvement_factor_required = 1.0;
  const bool poor_improvement =
      improvement_factor < improvement_factor_required;

  // Possibly abandon scaling if it's not improved equlibration significantly
  if (possibly_abandon_scaling && poor_improvement) {
    // Unscale the matrix
    for (int iCol = 0; iCol < numCol; iCol++) {
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        int iRow = Aindex[k];
        Avalue[k] /= (colScale[iCol] * rowScale[iRow]);
      }
    }
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Scaling: Improvement factor %0.4g < %0.4g required, so no "
                    "scaling applied",
                    improvement_factor, improvement_factor_required);
    scaleHighsModelInit(highs_model_object);
    return false;
  } else {
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Scaling: Improvement factor is %0.4g >= %0.4g so scale LP",
                    improvement_factor, improvement_factor_required);
#ifdef HiGHSDEV
    if (extreme_equilibration_improvement < 1.0) {
      HighsLogMessage(
          highs_model_object.options_.logfile, HighsMessageType::WARNING,
          "Scaling: Applying scaling with extreme improvement of %0.4g",
          extreme_equilibration_improvement);
    }
    if (mean_equilibration_improvement < 1.0) {
      HighsLogMessage(
          highs_model_object.options_.logfile, HighsMessageType::WARNING,
          "Scaling: Applying scaling with mean improvement of %0.4g",
          mean_equilibration_improvement);
    }
    if (matrix_value_ratio_improvement < 1.0) {
      HighsLogMessage(highs_model_object.options_.logfile,
                      HighsMessageType::WARNING,
                      "Scaling: Applying scaling with matrix value ratio "
                      "improvement of %0.4g",
                      matrix_value_ratio_improvement);
    }
    if (improvement_factor < 10 * improvement_factor_required) {
      HighsLogMessage(highs_model_object.options_.logfile,
                      HighsMessageType::WARNING,
                      "Scaling: Applying scaling with improvement factor %0.4g "
                      "< 10*(%0.4g) improvement",
                      improvement_factor, improvement_factor_required);
    }
#endif
  }
  return true;
}

bool maxValueScaleMatrix(HighsModelObject& highs_model_object) {
  int numCol = highs_model_object.simplex_lp_.numCol_;
  int numRow = highs_model_object.simplex_lp_.numRow_;
  vector<double>& colScale = highs_model_object.scale_.col_;
  vector<double>& rowScale = highs_model_object.scale_.row_;
  vector<int>& Astart = highs_model_object.simplex_lp_.Astart_;
  vector<int>& Aindex = highs_model_object.simplex_lp_.Aindex_;
  vector<double>& Avalue = highs_model_object.simplex_lp_.Avalue_;

  assert(highs_model_object.options_.simplex_scale_strategy ==
             SIMPLEX_SCALE_STRATEGY_015 ||
         highs_model_object.options_.simplex_scale_strategy ==
             SIMPLEX_SCALE_STRATEGY_0157);
  const double log2 = log(2.0);
  const double max_allow_scale =
      pow(2.0, highs_model_object.options_.allowed_simplex_matrix_scale_factor);
  const double min_allow_scale = 1 / max_allow_scale;

  const double min_allow_col_scale = min_allow_scale;
  const double max_allow_col_scale = max_allow_scale;
  const double min_allow_row_scale = min_allow_scale;
  const double max_allow_row_scale = max_allow_scale;

  double min_row_scale = HIGHS_CONST_INF;
  double max_row_scale = 0;
  double original_matrix_min_value = HIGHS_CONST_INF;
  double original_matrix_max_value = 0;
  // Determine the row scaling. Also determine the max/min row scaling
  // factors, and max/min original matrix values
  vector<double> row_max_value(numRow, 0);
  for (int iCol = 0; iCol < numCol; iCol++) {
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      const int iRow = Aindex[k];
      const double value = fabs(Avalue[k]);
      row_max_value[iRow] = max(row_max_value[iRow], value);
      original_matrix_min_value = min(original_matrix_min_value, value);
      original_matrix_max_value = max(original_matrix_max_value, value);
    }
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
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
  double min_col_scale = HIGHS_CONST_INF;
  double max_col_scale = 0;
  double matrix_min_value = HIGHS_CONST_INF;
  double matrix_max_value = 0;
  for (int iCol = 0; iCol < numCol; iCol++) {
    double col_max_value = 0;
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      const int iRow = Aindex[k];
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
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
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
  HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                  "Scaling: Factors are in [%0.4g, %0.4g] for columns and in "
                  "[%0.4g, %0.4g] for rows",
                  min_col_scale, max_col_scale, min_row_scale, max_row_scale);
  HighsLogMessage(
      highs_model_object.options_.logfile, HighsMessageType::INFO,
      "Scaling: Yields [min, max, ratio] matrix values of [%0.4g, %0.4g, "
      "%0.4g]; Originally [%0.4g, %0.4g, %0.4g]: Improvement of %0.4g",
      matrix_min_value, matrix_max_value, matrix_value_ratio,
      original_matrix_min_value, original_matrix_max_value,
      original_matrix_value_ratio, matrix_value_ratio_improvement);
  return true;
}

HighsStatus deleteScale(const HighsOptions& options, vector<double>& scale,
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
                         index_collection.set_num_entries_, 0,
                         index_collection.dimension_ - 1, true))
      return HighsStatus::Error;
  }
  if (from_k > to_k) return HighsStatus::OK;

  int delete_from_col;
  int delete_to_col;
  int keep_from_col;
  int keep_to_col = -1;
  int current_set_entry = 0;

  int col_dim = index_collection.dimension_;
  int new_num_col = 0;
  for (int k = from_k; k <= to_k; k++) {
    updateIndexCollectionOutInIndex(index_collection, delete_from_col,
                                    delete_to_col, keep_from_col, keep_to_col,
                                    current_set_entry);
    // Account for the initial columns being kept
    if (k == from_k) new_num_col = delete_from_col;
    if (delete_to_col >= col_dim - 1) break;
    assert(delete_to_col < col_dim);
    for (int col = keep_from_col; col <= keep_to_col; col++) {
      scale[new_num_col] = scale[col];
      new_num_col++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
  return HighsStatus::OK;
}

// PERMUTE:

void permuteSimplexLp(HighsModelObject& highs_model_object) {
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
#ifdef HiGHSDEV
  printf("Called permuteSimplexLp: simplex_lp_status.is_permuted = %d\n",
         simplex_lp_status.is_permuted);
#endif
  if (simplex_lp_status.is_permuted) return;

  int numCol = highs_model_object.simplex_lp_.numCol_;
  vector<int>& numColPermutation =
      highs_model_object.simplex_info_.numColPermutation_;
  vector<int>& Astart = highs_model_object.simplex_lp_.Astart_;
  vector<int>& Aindex = highs_model_object.simplex_lp_.Aindex_;
  vector<double>& Avalue = highs_model_object.simplex_lp_.Avalue_;
  vector<double>& colCost = highs_model_object.simplex_lp_.colCost_;
  vector<double>& colLower = highs_model_object.simplex_lp_.colLower_;
  vector<double>& colUpper = highs_model_object.simplex_lp_.colUpper_;

  // 2. Duplicate the original data to copy from
  vector<int> saveAstart = highs_model_object.simplex_lp_.Astart_;
  vector<int> saveAindex = highs_model_object.simplex_lp_.Aindex_;
  vector<double> saveAvalue = highs_model_object.simplex_lp_.Avalue_;
  vector<double> saveColCost = highs_model_object.simplex_lp_.colCost_;
  vector<double> saveColLower = highs_model_object.simplex_lp_.colLower_;
  vector<double> saveColUpper = highs_model_object.simplex_lp_.colUpper_;

  // 3. Generate the permuted matrix and corresponding vectors of column data
  int countX = 0;
  for (int i = 0; i < numCol; i++) {
    int fromCol = numColPermutation[i];
    Astart[i] = countX;
    for (int k = saveAstart[fromCol]; k < saveAstart[fromCol + 1]; k++) {
      Aindex[countX] = saveAindex[k];
      Avalue[countX] = saveAvalue[k];
      countX++;
    }
    colCost[i] = saveColCost[fromCol];
    colLower[i] = saveColLower[fromCol];
    colUpper[i] = saveColUpper[fromCol];
  }
  if (highs_model_object.scale_.is_scaled_) {
    // Permute any columns scaling factors
    vector<double>& colScale = highs_model_object.scale_.col_;
    vector<double> saveColScale = highs_model_object.scale_.col_;
    for (int i = 0; i < numCol; i++) {
      int fromCol = numColPermutation[i];
      colScale[i] = saveColScale[fromCol];
    }
  }
  assert(Astart[numCol] == countX);
  // Deduce the consequences of permuting the LP
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                        LpAction::PERMUTE);
}

#ifdef HiGHSDEV
// Only used to analyse the row and column status after Crash
void initialise_basic_index(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  int num_basic_variables = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int iVar = 0; iVar < numTot; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) {
      assert(num_basic_variables < simplex_lp.numRow_);
      simplex_basis.basicIndex_[num_basic_variables] = iVar;
      num_basic_variables++;
    }
  }
  assert(num_basic_variables == simplex_lp.numRow_);
}
#endif

void allocateWorkAndBaseArrays(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  // Allocate bounds and solution spaces
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  simplex_info.workCost_.resize(numTot);
  simplex_info.workDual_.resize(numTot);
  simplex_info.workShift_.resize(numTot);

  simplex_info.workLower_.resize(numTot);
  simplex_info.workUpper_.resize(numTot);
  simplex_info.workRange_.resize(numTot);
  simplex_info.workValue_.resize(numTot);

  // Feel that it should be possible to resize this with in dual
  // solver, and only if Devex is being used, but a pointer to it
  // needs to be set up when constructing HDual
  simplex_info.devex_index_.resize(numTot);

  simplex_info.baseLower_.resize(simplex_lp.numRow_);
  simplex_info.baseUpper_.resize(simplex_lp.numRow_);
  simplex_info.baseValue_.resize(simplex_lp.numRow_);
}

void initialiseValueAndNonbasicMove(HighsModelObject& highs_model_object) {
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const int numTot = highs_model_object.simplex_lp_.numCol_ +
                     highs_model_object.simplex_lp_.numRow_;
  for (int iVar = 0; iVar < numTot; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) {
      // Basic variable
      simplex_basis.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
    } else {
      // Nonbasic variable
      if (simplex_info.workLower_[iVar] == simplex_info.workUpper_[iVar]) {
        // Fixed
        simplex_info.workValue_[iVar] = simplex_info.workLower_[iVar];
        simplex_basis.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
      } else if (!highs_isInfinity(-simplex_info.workLower_[iVar])) {
        // Finite lower bound so boxed or lower
        if (!highs_isInfinity(simplex_info.workUpper_[iVar])) {
          // Finite upper bound so boxed
          if (simplex_basis.nonbasicMove_[iVar] == NONBASIC_MOVE_UP) {
            // Set at lower
            simplex_info.workValue_[iVar] = simplex_info.workLower_[iVar];
          } else if (simplex_basis.nonbasicMove_[iVar] == NONBASIC_MOVE_DN) {
            // Set at upper
            simplex_info.workValue_[iVar] = simplex_info.workUpper_[iVar];
          } else {
            // Invalid nonbasicMove: correct and set value at lower
            simplex_basis.nonbasicMove_[iVar] = NONBASIC_MOVE_UP;
            simplex_info.workValue_[iVar] = simplex_info.workLower_[iVar];
          }
        } else {
          // Lower
          simplex_info.workValue_[iVar] = simplex_info.workLower_[iVar];
          simplex_basis.nonbasicMove_[iVar] = NONBASIC_MOVE_UP;
        }
      } else if (!highs_isInfinity(simplex_info.workUpper_[iVar])) {
        // Upper
        simplex_info.workValue_[iVar] = simplex_info.workUpper_[iVar];
        simplex_basis.nonbasicMove_[iVar] = NONBASIC_MOVE_DN;
      } else {
        // FREE
        simplex_info.workValue_[iVar] = 0;
        simplex_basis.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
      }
    }
  }
}

void initialisePhase2ColBound(HighsModelObject& highs_model_object) {
  // Copy bounds and compute ranges
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    simplex_info.workLower_[iCol] = simplex_lp.colLower_[iCol];
    simplex_info.workUpper_[iCol] = simplex_lp.colUpper_[iCol];
    simplex_info.workRange_[iCol] =
        simplex_info.workUpper_[iCol] - simplex_info.workLower_[iCol];
  }
}

void initialisePhase2RowBound(HighsModelObject& highs_model_object) {
  // Copy bounds and compute ranges
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = simplex_lp.numCol_ + row;
    simplex_info.workLower_[var] = -simplex_lp.rowUpper_[row];
    simplex_info.workUpper_[var] = -simplex_lp.rowLower_[row];
    simplex_info.workRange_[var] =
        simplex_info.workUpper_[var] - simplex_info.workLower_[var];
  }
}

void initialiseBound(HighsModelObject& highs_model_object, int phase) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  // Initialise the Phase 2 bounds (and ranges). NB Phase 2 bounds
  // necessary to compute Phase 1 bounds
  initialisePhase2ColBound(highs_model_object);
  initialisePhase2RowBound(highs_model_object);
  if (phase == 2) return;

  // The dual objective is the sum of products of primal and dual
  // values for nonbasic variables. For dual simplex phase 1, the
  // primal bounds are set so that when the dual value is feasible, the
  // primal value is set to zero. Otherwise the value is +1/-1
  // according to the required sign of the dual, except for free
  // variables, where the bounds are [-1000, 1000]. Hence the dual
  // objective is the negation of the sum of infeasibilities, unless there are
  // free In Phase 1: change to dual phase 1 bound.
  const double inf = HIGHS_CONST_INF;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++) {
    if (simplex_info.workLower_[i] == -inf &&
        simplex_info.workUpper_[i] == inf) {
      // Don't change for row variables: they should never become
      // nonbasic when starting from a logical basis, and no crash
      // should make a free row nonbasic, but could an advanced basis
      // make a free row nonbasic.
      // But what it it happened?
      if (i >= simplex_lp.numCol_) continue;
      simplex_info.workLower_[i] = -1000,
      simplex_info.workUpper_[i] = 1000;  // FREE
    } else if (simplex_info.workLower_[i] == -inf) {
      simplex_info.workLower_[i] = -1, simplex_info.workUpper_[i] = 0;  // UPPER
    } else if (simplex_info.workUpper_[i] == inf) {
      simplex_info.workLower_[i] = 0, simplex_info.workUpper_[i] = 1;  // LOWER
    } else {
      simplex_info.workLower_[i] = 0,
      simplex_info.workUpper_[i] = 0;  // BOXED or FIXED
    }
    simplex_info.workRange_[i] =
        simplex_info.workUpper_[i] - simplex_info.workLower_[i];
  }
}

void initialisePhase2ColCost(HighsModelObject& highs_model_object) {
  // Copy the Phase 2 cost and zero the shift
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  for (int col = 0; col < simplex_lp.numCol_; col++) {
    int var = col;
    simplex_info.workCost_[var] =
        (int)simplex_lp.sense_ * simplex_lp.colCost_[col];
    simplex_info.workShift_[var] = 0;
  }
}

void initialisePhase2RowCost(HighsModelObject& highs_model_object) {
  // Zero the cost and shift
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  for (int iVar = simplex_lp.numCol_;
       iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    simplex_info.workCost_[iVar] = 0;
    simplex_info.workShift_[iVar] = 0;
  }
}

void initialiseCost(HighsModelObject& highs_model_object, int perturb) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
#ifdef HiGHSDEV
  HighsSimplexAnalysis* analysis = &highs_model_object.simplex_analysis_;
#endif
  // Copy the cost
  initialisePhase2ColCost(highs_model_object);
  initialisePhase2RowCost(highs_model_object);
  // See if we want to skip perturbation
  simplex_info.costs_perturbed = 0;
  if (perturb == 0 ||
      simplex_info.dual_simplex_cost_perturbation_multiplier == 0)
    return;
  simplex_info.costs_perturbed = 1;

  // Perturb the original costs, scale down if is too big
#ifdef HiGHSDEV
  printf("grep_DuPtrb: Cost perturbation for %s\n",
         highs_model_object.simplex_lp_.model_name_.c_str());
  int num_original_nonzero_cost = 0;
#endif
  double bigc = 0;
  for (int i = 0; i < simplex_lp.numCol_; i++) {
    const double abs_cost = fabs(simplex_info.workCost_[i]);
    bigc = max(bigc, abs_cost);
#ifdef HiGHSDEV
    if (abs_cost) num_original_nonzero_cost++;
#endif
  }
#ifdef HiGHSDEV
  const int pct0 = (100 * num_original_nonzero_cost) / simplex_lp.numCol_;
  double average_cost = 0;
  if (num_original_nonzero_cost) {
    average_cost = bigc / num_original_nonzero_cost;
  } else {
    printf("grep_DuPtrb:    STRANGE initial workCost has non nonzeros\n");
  }
  printf(
      "grep_DuPtrb:    Initially have %d nonzero costs (%3d%%) with bigc = %g "
      "and average = %g\n",
      num_original_nonzero_cost, pct0, bigc, average_cost);
#endif
  if (bigc > 100) {
    bigc = sqrt(sqrt(bigc));
#ifdef HiGHSDEV
    printf("grep_DuPtrb:    Large so set bigc = sqrt(bigc) = %g\n", bigc);
#endif
  }

  // If there's few boxed variables, we will just use simple perturbation
  double boxedRate = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++)
    boxedRate += (simplex_info.workRange_[i] < 1e30);
  boxedRate /= numTot;
  if (boxedRate < 0.01) {
    bigc = min(bigc, 1.0);
#ifdef HiGHSDEV
    printf(
        "grep_DuPtrb:    small boxedRate (%g) so set bigc = min(bigc, 1.0) = "
        "%g\n",
        boxedRate, bigc);
#endif
  }
  // Determine the perturbation base
  double base = 5e-7 * bigc;
#ifdef HiGHSDEV
  printf("grep_DuPtrb:    Perturbation base = %g\n", base);
#endif

  // Now do the perturbation
  for (int i = 0; i < simplex_lp.numCol_; i++) {
    double lower = simplex_lp.colLower_[i];
    double upper = simplex_lp.colUpper_[i];
    double xpert = (fabs(simplex_info.workCost_[i]) + 1) * base *
                   simplex_info.dual_simplex_cost_perturbation_multiplier *
                   (1 + simplex_info.numTotRandomValue_[i]);
#ifdef HiGHSDEV
    const double previous_cost = simplex_info.workCost_[i];
#endif
    if (lower <= -HIGHS_CONST_INF && upper >= HIGHS_CONST_INF) {
      // Free - no perturb
    } else if (upper >= HIGHS_CONST_INF) {  // Lower
      simplex_info.workCost_[i] += xpert;
    } else if (lower <= -HIGHS_CONST_INF) {  // Upper
      simplex_info.workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      simplex_info.workCost_[i] +=
          (simplex_info.workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
#ifdef HiGHSDEV
    const double perturbation1 =
        fabs(simplex_info.workCost_[i] - previous_cost);
    if (perturbation1)
      updateValueDistribution(perturbation1,
                              analysis->cost_perturbation1_distribution);
#endif
  }
  for (int i = simplex_lp.numCol_; i < numTot; i++) {
    double perturbation2 =
        (0.5 - simplex_info.numTotRandomValue_[i]) *
        simplex_info.dual_simplex_cost_perturbation_multiplier * 1e-12;
    simplex_info.workCost_[i] += perturbation2;
#ifdef HiGHSDEV
    perturbation2 = fabs(perturbation2);
    updateValueDistribution(perturbation2,
                            analysis->cost_perturbation2_distribution);
#endif
  }
}

void populateWorkArrays(HighsModelObject& highs_model_object) {
  // Initialize the values
  initialiseCost(highs_model_object);
  initialiseBound(highs_model_object);
  initialiseValueAndNonbasicMove(highs_model_object);
}

#ifdef HiGHSDEV
void reportSimplexProfiling(HighsModelObject& highs_model_object) {
  HighsTimer& timer = highs_model_object.timer_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;
  SimplexTimer simplex_timer;
  int simplex_strategy_for_profiling = simplex_info.simplex_strategy;
  if (simplex_strategy_for_profiling == SIMPLEX_STRATEGY_CHOOSE) {
    // Simplex strategy not chosen - probably due to solution after
    // postsolve being optimal - so profile as if
    // SIMPLEX_STRATEGY_DUAL_PLAIN has been used
    simplex_strategy_for_profiling = SIMPLEX_STRATEGY_DUAL_PLAIN;
  }

  if (simplex_strategy_for_profiling == SIMPLEX_STRATEGY_PRIMAL) {
    if (simplex_info.report_simplex_inner_clock) {
      simplex_timer.reportSimplexInnerClock(analysis.thread_simplex_clocks[0]);
    }
  } else if (simplex_strategy_for_profiling == SIMPLEX_STRATEGY_DUAL_PLAIN) {
    if (simplex_info.report_simplex_inner_clock) {
      simplex_timer.reportSimplexInnerClock(analysis.thread_simplex_clocks[0]);
      simplex_timer.reportSimplexChuzc3Clock(analysis.thread_simplex_clocks[0]);
    }
    if (simplex_info.report_simplex_outer_clock) {
      simplex_timer.reportDualSimplexIterateClock(
          analysis.thread_simplex_clocks[0]);
      simplex_timer.reportDualSimplexOuterClock(
          analysis.thread_simplex_clocks[0]);
    }
  }

  if (simplex_strategy_for_profiling == SIMPLEX_STRATEGY_DUAL_MULTI) {
    if (simplex_info.report_simplex_inner_clock) {
      simplex_timer.reportSimplexMultiInnerClock(
          analysis.thread_simplex_clocks[0]);
    }
    printf("PAMI   %-20s    CUTOFF  %6g    PERSISTENSE  %6g\n",
           highs_model_object.lp_.model_name_.c_str(), simplex_info.pami_cutoff,
           highs_model_object.iteration_counts_.simplex /
               (1.0 + simplex_info.multi_iteration));
  }

  if (simplex_info.report_simplex_phases_clock) {
    simplex_timer.reportSimplexTotalClock(analysis.thread_simplex_clocks[0]);
    simplex_timer.reportSimplexPhasesClock(analysis.thread_simplex_clocks[0]);
  }

  if (simplex_info.analyse_invert_time) {
    double current_run_highs_time = timer.readRunHighsClock();
    simplex_info.total_inverts = analysis.simplexTimerNumCall(InvertClock);
    simplex_info.total_invert_time = analysis.simplexTimerRead(InvertClock);
    printf(
        "Time: Total inverts =  %4d; Total invert  time = %11.4g of Total time "
        "= %11.4g",
        simplex_info.total_inverts, simplex_info.total_invert_time,
        current_run_highs_time);
    if (current_run_highs_time > 0.001) {
      printf(" (%6.2f%%)\n",
             (100 * simplex_info.total_invert_time) / current_run_highs_time);
    } else {
      printf("\n");
    }
  }
  /*
  if (simplex_info.analyse_rebuild_time) {
    double current_run_highs_time = timer.readRunHighsClock();
    HighsClockRecord totalRebuildClock;
    timer.clockInit(totalRebuildClock);
    timer.clockAdd(totalRebuildClock,
                   simplex_info.clock_[IterateDualRebuildClock]);
    timer.clockAdd(totalRebuildClock,
                   simplex_info.clock_[IteratePrimalRebuildClock]);
    int totalRebuilds = 0;
    double totalRebuildTime = 0;
    printf("Time: Total rebuild time = %11.4g (%4d) of Total time = %11.4g",
           totalRebuildTime, totalRebuilds, current_run_highs_time);
    if (current_run_highs_time > 0.001) {
      printf(" (%6.2f%%)\n", (100 * totalRebuildTime) / current_run_highs_time);
    } else {
      printf("\n");
    }
  }
  */
}
#endif

void setRunQuiet(HighsModelObject& highs_model_object) {
  highs_model_object.simplex_info_.run_quiet =
      highs_model_object.options_.output == NULL &&
      highs_model_object.options_.logfile == NULL;
}

double computeBasisCondition(const HighsModelObject& highs_model_object) {
  int solver_num_row = highs_model_object.simplex_lp_.numRow_;
  int solver_num_col = highs_model_object.simplex_lp_.numCol_;
  vector<double> bs_cond_x;
  vector<double> bs_cond_y;
  vector<double> bs_cond_z;
  vector<double> bs_cond_w;
  HVector row_ep;
  row_ep.setup(solver_num_row);

  const HFactor& factor = highs_model_object.factor_;
  const int* Astart = &highs_model_object.simplex_lp_.Astart_[0];
  const double* Avalue = &highs_model_object.simplex_lp_.Avalue_[0];
  // Compute the Hager condition number estimate for the basis matrix
  const double NoDensity = 1;
  bs_cond_x.resize(solver_num_row);
  bs_cond_y.resize(solver_num_row);
  bs_cond_z.resize(solver_num_row);
  bs_cond_w.resize(solver_num_row);
  // x = ones(n,1)/n;
  // y = A\x;
  double mu = 1.0 / solver_num_row;
  double norm_Binv;
  for (int r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = mu;
  row_ep.clear();
  for (int r_n = 0; r_n < solver_num_row; r_n++) {
    double value = bs_cond_x[r_n];
    if (value) {
      row_ep.index[row_ep.count] = r_n;
      row_ep.array[r_n] = value;
      row_ep.count++;
    }
  }
  for (int ps_n = 1; ps_n <= 5; ps_n++) {
    row_ep.packFlag = false;
    factor.ftran(row_ep, NoDensity);
    // zeta = sign(y);
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
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
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
      double value = bs_cond_w[r_n];
      if (value) {
        row_ep.index[row_ep.count] = r_n;
        row_ep.array[r_n] = value;
        row_ep.count++;
      }
    }
    row_ep.packFlag = false;
    factor.btran(row_ep, NoDensity);
    double norm_z = 0.0;
    double ztx = 0.0;
    norm_Binv = 0.0;
    int argmax_z = -1;
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
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
    for (int r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = 0.0;
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = argmax_z;
    row_ep.array[argmax_z] = 1.0;
    bs_cond_x[argmax_z] = 1.0;
  }
  double norm_B = 0.0;
  for (int r_n = 0; r_n < solver_num_row; r_n++) {
    int vr_n = highs_model_object.simplex_basis_.basicIndex_[r_n];
    double c_norm = 0.0;
    if (vr_n < solver_num_col)
      for (int el_n = Astart[vr_n]; el_n < Astart[vr_n + 1]; el_n++)
        c_norm += fabs(Avalue[el_n]);
    else
      c_norm += 1.0;
    norm_B = max(c_norm, norm_B);
  }
  double cond_B = norm_Binv * norm_B;
  return cond_B;
}

void flip_bound(HighsModelObject& highs_model_object, int iCol) {
  int* nonbasicMove = &highs_model_object.simplex_basis_.nonbasicMove_[0];
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const int move = nonbasicMove[iCol] = -nonbasicMove[iCol];
  simplex_info.workValue_[iCol] =
      move == 1 ? simplex_info.workLower_[iCol] : simplex_info.workUpper_[iCol];
}

int simplexHandleRankDeficiency(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HFactor& factor = highs_model_object.factor_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  int rank_deficiency = factor.rank_deficiency;
  vector<int>& noPvC = factor.noPvC;
  printf("Returned %d = factor.build();\n", rank_deficiency);
  fflush(stdout);
  vector<int> basicRows;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  basicRows.resize(numTot);
  //    printf("Before - simplex_basis.basicIndex_:");
  // for (int iRow=0; iRow<simplex_lp.numRow_; iRow++)
  //    printf(" %2d", simplex_basis.basicIndex_[iRow]); printf("\n");
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
    basicRows[simplex_basis.basicIndex_[iRow]] = iRow;
  for (int k = 0; k < rank_deficiency; k++) {
    //      printf("noPvR[%2d] = %d; noPvC[%2d] = %d; \n", k, factor.noPvR[k],
    //      k, noPvC[k]);fflush(stdout);
    int columnIn = simplex_lp.numCol_ + factor.noPvR[k];
    int columnOut = noPvC[k];
    int rowOut = basicRows[columnOut];
    //      printf("columnIn = %6d; columnOut = %6d; rowOut = %6d [%11.4g,
    //      %11.4g]\n", columnIn, columnOut, rowOut,
    //      simplex_info.workLower_[columnOut],
    //      simplex_info.workUpper_[columnOut]);
    assert(simplex_basis.basicIndex_[rowOut] == columnOut);
    if (highs_model_object.simplex_info_.initialised) {
      int sourceOut = setSourceOutFmBd(highs_model_object, columnOut);
      update_pivots(highs_model_object, columnIn, rowOut, sourceOut);
      update_matrix(highs_model_object, columnIn, columnOut);
    } else {
      simplex_basis.basicIndex_[rowOut] = columnIn;
      simplex_basis.nonbasicFlag_[columnIn] = NONBASIC_FLAG_FALSE;
      simplex_basis.nonbasicFlag_[columnOut] = NONBASIC_FLAG_TRUE;
    }
  }
  //    printf("After  - simplex_basis.basicIndex_:");
  // for (int iRow=0; iRow<simplex_lp.numRow_; iRow++)
  //    printf(" %2d", simplex_basis.basicIndex_[iRow]); printf("\n");
  debugCheckInvert(highs_model_object.options_, highs_model_object.factor_);
  return 0;
}

int computeFactor(HighsModelObject& highs_model_object) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  HFactor& factor = highs_model_object.factor_;
#ifdef HiGHSDEV
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;
  double tt0 = 0;
  if (simplex_info.analyse_invert_time)
    tt0 = analysis.simplexTimerRead(InvertClock);
#endif
  HighsTimerClock* factor_timer_clock_pointer = NULL;
  // TODO Understand why handling noPvC and noPvR in what seem to be
  // different ways ends up equivalent.
#ifdef HiGHSDEV
  int thread_id = 0;
#ifdef OPENMP
  thread_id = omp_get_thread_num();
  //  printf("Hello world from computeFactor: thread %d\n", thread_id);
#endif
  factor_timer_clock_pointer =
      highs_model_object.simplex_analysis_.getThreadFactorTimerClockPtr(
          thread_id);
#endif
  int rank_deficiency = factor.build(factor_timer_clock_pointer);
  if (rank_deficiency) {
    simplexHandleRankDeficiency(highs_model_object);
  }
#ifdef HiGHSDEV
  if (simplex_info.analyse_invert_form) {
    const bool report_kernel = false;
    simplex_info.num_invert++;
    assert(factor.basis_matrix_num_el);
    double invert_fill_factor =
        ((1.0 * factor.invert_num_el) / factor.basis_matrix_num_el);
    if (report_kernel) printf("INVERT fill = %6.2f", invert_fill_factor);
    simplex_info.sum_invert_fill_factor += invert_fill_factor;
    simplex_info.running_average_invert_fill_factor =
        0.95 * simplex_info.running_average_invert_fill_factor +
        0.05 * invert_fill_factor;

    double kernel_relative_dim =
        (1.0 * factor.kernel_dim) / highs_model_object.simplex_lp_.numRow_;
    if (report_kernel) printf("; kernel dim = %11.4g", kernel_relative_dim);
    if (factor.kernel_dim) {
      simplex_info.num_kernel++;
      simplex_info.max_kernel_dim =
          max(kernel_relative_dim, simplex_info.max_kernel_dim);
      simplex_info.sum_kernel_dim += kernel_relative_dim;
      simplex_info.running_average_kernel_dim =
          0.95 * simplex_info.running_average_kernel_dim +
          0.05 * kernel_relative_dim;

      int kernel_invert_num_el =
          factor.invert_num_el -
          (factor.basis_matrix_num_el - factor.kernel_num_el);
      assert(factor.kernel_num_el);
      double kernel_fill_factor =
          (1.0 * kernel_invert_num_el) / factor.kernel_num_el;
      simplex_info.sum_kernel_fill_factor += kernel_fill_factor;
      simplex_info.running_average_kernel_fill_factor =
          0.95 * simplex_info.running_average_kernel_fill_factor +
          0.05 * kernel_fill_factor;
      if (report_kernel) printf("; fill = %6.2f", kernel_fill_factor);
      if (kernel_relative_dim >
          simplex_info.major_kernel_relative_dim_threshhold) {
        simplex_info.num_major_kernel++;
        simplex_info.sum_major_kernel_fill_factor += kernel_fill_factor;
        simplex_info.running_average_major_kernel_fill_factor =
            0.95 * simplex_info.running_average_major_kernel_fill_factor +
            0.05 * kernel_fill_factor;
      }
    }
    if (report_kernel) printf("\n");
  }
#endif
  simplex_info.update_count = 0;
  debugCheckInvert(highs_model_object.options_, highs_model_object.factor_);

#ifdef HiGHSDEV
  if (simplex_info.analyse_invert_time) {
    simplex_info.total_inverts = analysis.simplexTimerNumCall(InvertClock);
    simplex_info.total_invert_time = analysis.simplexTimerRead(InvertClock);
    const double invert_time = simplex_info.total_invert_time - tt0;
    printf(
        "           INVERT  %4d     on iteration %9d: INVERT  time = %11.4g; "
        "Total INVERT  time = %11.4g\n",
        simplex_info.total_inverts,
        highs_model_object.iteration_counts_.simplex, invert_time,
        simplex_info.total_invert_time);
  }
#endif

  // Now have a representation of B^{-1}, and it is fresh!
  simplex_lp_status.has_invert = true;
  simplex_lp_status.has_fresh_invert = true;

#ifdef HiGHSDEV
  if (simplex_info.analyse_invert_condition) {
    analysis.simplexTimerStart(BasisConditionClock);
    simplex_info.invert_condition = computeBasisCondition(highs_model_object);
    analysis.simplexTimerStop(BasisConditionClock);
  }
#endif

  return rank_deficiency;
}

// Compute the primal values (in baseValue) and set the lower and upper bounds
// of basic variables
void computePrimal(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HMatrix& matrix = highs_model_object.matrix_;
  HFactor& factor = highs_model_object.factor_;
  HighsSimplexAnalysis* analysis = &highs_model_object.simplex_analysis_;
  // Setup a local buffer for the values of basic variables
  HVector primal_col;
  primal_col.setup(simplex_lp.numRow_);
  primal_col.clear();
  for (int i = 0; i < simplex_lp.numCol_ + simplex_lp.numRow_; i++) {
    if (simplex_basis.nonbasicFlag_[i] && simplex_info.workValue_[i] != 0) {
      matrix.collect_aj(primal_col, i, simplex_info.workValue_[i]);
    }
  }
  // If debugging, take a copy of the RHS
  vector<double> debug_primal_rhs;
  if (highs_model_object.options_.highs_debug_level >= HIGHS_DEBUG_LEVEL_COSTLY)
    debug_primal_rhs = primal_col.array;

  // It's possible that the buffer has no nonzeros, so performing
  // FTRAN is unnecessary. Not much of a saving, but the zero density
  // looks odd in the analysis!
  if (primal_col.count) {
    factor.ftran(primal_col, analysis->primal_col_density,
                 analysis->pointer_serial_factor_clocks);
    const double local_primal_col_density =
        (double)primal_col.count / simplex_lp.numRow_;
    analysis->updateOperationResultDensity(local_primal_col_density,
                                           analysis->primal_col_density);
  }
  for (int i = 0; i < simplex_lp.numRow_; i++) {
    int iCol = simplex_basis.basicIndex_[i];
    simplex_info.baseValue_[i] = -primal_col.array[i];
    simplex_info.baseLower_[i] = simplex_info.workLower_[iCol];
    simplex_info.baseUpper_[i] = simplex_info.workUpper_[iCol];
  }
  debugComputePrimal(highs_model_object, debug_primal_rhs);
  // Now have basic primals
  simplex_lp_status.has_basic_primal_values = true;
}

void computeSimplexInfeasible(HighsModelObject& highs_model_object) {
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;
  analysis.simplexTimerStart(ComputePrIfsClock);
  computeSimplexPrimalInfeasible(highs_model_object);
  analysis.simplexTimerStop(ComputePrIfsClock);

  analysis.simplexTimerStart(ComputeDuIfsClock);
  computeSimplexDualInfeasible(highs_model_object);
  analysis.simplexTimerStop(ComputeDuIfsClock);
}

void computeSimplexPrimalInfeasible(HighsModelObject& highs_model_object) {
  // Computes num/max/sum of primal infeasibliities according to the
  // simplex bounds. This is used to determine optimality in dual
  // phase 1 and dual phase 2, albeit using different bounds in
  // workLower/Upper.
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  const double scaled_primal_feasibility_tolerance =
      highs_model_object.scaled_solution_params_.primal_feasibility_tolerance;
  int& num_primal_infeasibilities =
      highs_model_object.simplex_info_.num_primal_infeasibilities;
  double& max_primal_infeasibility =
      highs_model_object.simplex_info_.max_primal_infeasibility;
  double& sum_primal_infeasibilities =
      highs_model_object.simplex_info_.sum_primal_infeasibilities;
  num_primal_infeasibilities = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibilities = 0;

  for (int i = 0; i < simplex_lp.numCol_ + simplex_lp.numRow_; i++) {
    if (simplex_basis.nonbasicFlag_[i]) {
      // Nonbasic column
      double value = simplex_info.workValue_[i];
      double lower = simplex_info.workLower_[i];
      double upper = simplex_info.workUpper_[i];
      double primal_infeasibility = max(lower - value, value - upper);
      if (primal_infeasibility > 0) {
        if (primal_infeasibility > scaled_primal_feasibility_tolerance)
          num_primal_infeasibilities++;
        max_primal_infeasibility =
            std::max(primal_infeasibility, max_primal_infeasibility);
        sum_primal_infeasibilities += primal_infeasibility;
      }
    }
  }
  for (int i = 0; i < simplex_lp.numRow_; i++) {
    // Basic variable
    double value = simplex_info.baseValue_[i];
    double lower = simplex_info.baseLower_[i];
    double upper = simplex_info.baseUpper_[i];
    double primal_infeasibility = max(lower - value, value - upper);
    if (primal_infeasibility > 0) {
      if (primal_infeasibility > scaled_primal_feasibility_tolerance)
        num_primal_infeasibilities++;
      max_primal_infeasibility =
          std::max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibilities += primal_infeasibility;
    }
  }
}

void computeSimplexDualInfeasible(HighsModelObject& highs_model_object) {
  // Computes num/max/sum of dual infeasibilities in phase 1 and phase
  // 2 according to nonbasicMove. The bounds are only used to identify
  // free variables. Fixed variables are assumed to have
  // nonbasicMove=0 so that no dual infeasibility is counted for them.
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  const double scaled_dual_feasibility_tolerance =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;
  // Possibly verify that nonbasicMove is correct for fixed variables
  debugFixedNonbasicMove(highs_model_object);

  int& num_dual_infeasibilities =
      highs_model_object.simplex_info_.num_dual_infeasibilities;
  double& max_dual_infeasibility =
      highs_model_object.simplex_info_.max_dual_infeasibility;
  double& sum_dual_infeasibilities =
      highs_model_object.simplex_info_.sum_dual_infeasibilities;
  num_dual_infeasibilities = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibilities = 0;

  for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double dual = simplex_info.workDual_[iVar];
    const double lower = simplex_info.workLower_[iVar];
    const double upper = simplex_info.workUpper_[iVar];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else {
      // Not free: any dual infeasibility is given by the dual value
      // signed by nonbasicMove
      dual_infeasibility = -simplex_basis.nonbasicMove_[iVar] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
}

void computeSimplexLpDualInfeasible(HighsModelObject& highs_model_object) {
  // Compute num/max/sum of dual infeasibliities according to the
  // bounds of the simplex LP. Assumes that boxed variables have
  // primal variable at the bound corresponding to the sign of the
  // dual so should only be used in dual phase 1 - where it's only
  // used for reporting after rebuilds and to determine whether the LP
  // is dual infeasible and, hence, primal unbounded.
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  // Possibly verify that nonbasicMove is correct for fixed variables
  debugFixedNonbasicMove(highs_model_object);
  const double scaled_dual_feasibility_tolerance =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;
  int& num_dual_infeasibilities =
      highs_model_object.scaled_solution_params_.num_dual_infeasibilities;
  double& max_dual_infeasibility =
      highs_model_object.scaled_solution_params_.max_dual_infeasibility;
  double& sum_dual_infeasibilities =
      highs_model_object.scaled_solution_params_.sum_dual_infeasibilities;
  num_dual_infeasibilities = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibilities = 0;

  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    int iVar = iCol;
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double dual = simplex_info.workDual_[iVar];
    const double lower = simplex_lp.colLower_[iCol];
    const double upper = simplex_lp.colUpper_[iCol];
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
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    int iVar = simplex_lp.numCol_ + iRow;
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic row
    const double dual = -simplex_info.workDual_[iVar];
    const double lower = simplex_lp.rowLower_[iRow];
    const double upper = simplex_lp.rowUpper_[iRow];
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
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
}

void copySimplexInfeasible(HighsModelObject& highs_model_object) {
  copySimplexPrimalInfeasible(highs_model_object);
  copySimplexDualInfeasible(highs_model_object);
}

void copySimplexPrimalInfeasible(HighsModelObject& highs_model_object) {
  highs_model_object.scaled_solution_params_.num_primal_infeasibilities =
      highs_model_object.simplex_info_.num_primal_infeasibilities;
  highs_model_object.scaled_solution_params_.max_primal_infeasibility =
      highs_model_object.simplex_info_.max_primal_infeasibility;
  highs_model_object.scaled_solution_params_.sum_primal_infeasibilities =
      highs_model_object.simplex_info_.sum_primal_infeasibilities;
}

void copySimplexDualInfeasible(HighsModelObject& highs_model_object) {
  highs_model_object.scaled_solution_params_.num_dual_infeasibilities =
      highs_model_object.simplex_info_.num_dual_infeasibilities;
  highs_model_object.scaled_solution_params_.max_dual_infeasibility =
      highs_model_object.simplex_info_.max_dual_infeasibility;
  highs_model_object.scaled_solution_params_.sum_dual_infeasibilities =
      highs_model_object.simplex_info_.sum_dual_infeasibilities;
}

void computeDualInfeasibleWithFlips(HighsModelObject& highs_model_object) {
  // Computes num/max/sum of dual infeasibliities according to
  // nonbasicMove, using the bounds only to identify free variables
  // and non-boxed. Fixed variables are assumed to have nonbasicMove=0
  // so that no dual infeasibility is counted for them. Indeed, when
  // called from cleanup() at the end of dual phase 1, nonbasicMove
  // relates to the phase 1 bounds, but workLower and workUpper will
  // have been set to phase 2 values!
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSolutionParams& scaled_solution_params =
      highs_model_object.scaled_solution_params_;
  const double scaled_dual_feasibility_tolerance =
      scaled_solution_params.dual_feasibility_tolerance;
  // Possibly verify that nonbasicMove is correct for fixed variables
  debugFixedNonbasicMove(highs_model_object);

  int num_dual_infeasibilities = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibilities = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;

  for (int iVar = 0; iVar < numTot; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double lower = simplex_info.workLower_[iVar];
    const double upper = simplex_info.workUpper_[iVar];
    const double dual = simplex_info.workDual_[iVar];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else if (highs_isInfinity(-lower) || highs_isInfinity(upper)) {
      // Not free or boxed: any dual infeasibility is given by value
      // signed by nonbasicMove.
      //
      // For boxed variables, nonbasicMove may have the wrong sign for
      // dual, but nonbasicMove and the primal value can be flipped to
      // achieve dual feasiblility.
      dual_infeasibility = -simplex_basis.nonbasicMove_[iVar] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
  scaled_solution_params.num_dual_infeasibilities = num_dual_infeasibilities;
  scaled_solution_params.max_dual_infeasibility = max_dual_infeasibility;
  scaled_solution_params.sum_dual_infeasibilities = sum_dual_infeasibilities;
}

void choosePriceTechnique(const int price_strategy, const double row_ep_density,
                          bool& use_col_price, bool& use_row_price_w_switch) {
  // By default switch to column PRICE when pi_p has at least this
  // density
  const double density_for_column_price_switch = 0.75;
  use_col_price =
      (price_strategy == SIMPLEX_PRICE_STRATEGY_COL) ||
      (price_strategy == SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH &&
       row_ep_density > density_for_column_price_switch);
  use_row_price_w_switch =
      price_strategy == SIMPLEX_PRICE_STRATEGY_ROW_SWITCH ||
      price_strategy == SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH;
}

void computeTableauRowFromPiP(HighsModelObject& highs_model_object,
                              const HVector& row_ep, HVector& row_ap) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const HMatrix* matrix = &highs_model_object.matrix_;
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;

  const int solver_num_row = highs_model_object.simplex_lp_.numRow_;
  const double local_density = 1.0 * row_ep.count / solver_num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  choosePriceTechnique(simplex_info.price_strategy, local_density,
                       use_col_price, use_row_price_w_switch);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations) {
    if (use_col_price) {
      analysis.operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep,
                                     0.0);
      analysis.num_col_price++;
    } else if (use_row_price_w_switch) {
      analysis.operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep,
                                     analysis.row_ep_density);
      analysis.num_row_price_with_switch++;
    } else {
      analysis.operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep,
                                     analysis.row_ep_density);
      analysis.num_row_price++;
    }
  }
#endif
  analysis.simplexTimerStart(PriceClock);
  row_ap.clear();
  if (use_col_price) {
    // Perform column-wise PRICE
    matrix->priceByColumn(row_ap, row_ep);
  } else if (use_row_price_w_switch) {
    // Perform hyper-sparse row-wise PRICE, but switch if the density of row_ap
    // becomes extreme
    const double switch_density = matrix->hyperPRICE;
    matrix->priceByRowSparseResultWithSwitch(
        row_ap, row_ep, analysis.row_ap_density, 0, switch_density);
  } else {
    // Perform hyper-sparse row-wise PRICE
    matrix->priceByRowSparseResult(row_ap, row_ep);
  }

  const int solver_num_col = highs_model_object.simplex_lp_.numCol_;
  if (use_col_price) {
    // Column-wise PRICE computes components of row_ap corresponding
    // to basic variables, so zero these by exploiting the fact that,
    // for basic variables, nonbasicFlag[*]=0
    const int* nonbasicFlag =
        &highs_model_object.simplex_basis_.nonbasicFlag_[0];
    for (int col = 0; col < solver_num_col; col++)
      row_ap.array[col] = nonbasicFlag[col] * row_ap.array[col];
  }
#ifdef HiGHSDEV
  // Possibly analyse the error in the result of PRICE
  const bool analyse_price_error = false;
  if (analyse_price_error) matrix->price_er_ck(row_ap, row_ep);
#endif
  // Update the record of average row_ap density
  const double local_row_ap_density = (double)row_ap.count / solver_num_col;
  analysis.updateOperationResultDensity(local_row_ap_density,
                                        analysis.row_ap_density);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis.operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ap);
#endif
  analysis.simplexTimerStop(PriceClock);
}

void computeDual(HighsModelObject& highs_model_object) {
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HMatrix& matrix = highs_model_object.matrix_;
  HFactor& factor = highs_model_object.factor_;

  // Create a local buffer for the pi vector
  HVector dual_col;
  dual_col.setup(simplex_lp.numRow_);
  dual_col.clear();
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    const double value =
        simplex_info.workCost_[simplex_basis.basicIndex_[iRow]] +
        simplex_info.workShift_[simplex_basis.basicIndex_[iRow]];
    if (value) {
      dual_col.count++;
      dual_col.index[iRow] = iRow;
      dual_col.array[iRow] = value;
    }
  }
  // If debugging, take a copy of the basic costs and any previous duals
  vector<double> debug_previous_workDual;
  vector<double> debug_basic_costs;
  if (highs_model_object.options_.highs_debug_level >=
      HIGHS_DEBUG_LEVEL_COSTLY) {
    debug_basic_costs = dual_col.array;
    if (simplex_lp_status.has_nonbasic_dual_values)
      debug_previous_workDual = simplex_info.workDual_;
  }
  // Copy the costs in case the basic costs are all zero
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++)
    simplex_info.workDual_[i] = simplex_info.workCost_[i];
  if (dual_col.count) {
    // RHS of row dual calculation is nonzero
#ifdef HiGHSDEV
    if (simplex_info.analyse_iterations)
      analysis.operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_FULL,
                                     dual_col, analysis.dual_col_density);
#endif
    factor.btran(dual_col, analysis.dual_col_density,
                 analysis.pointer_serial_factor_clocks);
#ifdef HiGHSDEV
    if (simplex_info.analyse_iterations)
      analysis.operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_FULL,
                                    dual_col);
#endif
    const double local_dual_col_density =
        (double)dual_col.count / simplex_lp.numRow_;
    analysis.updateOperationResultDensity(local_dual_col_density,
                                          analysis.dual_col_density);
    // Create a local buffer for the values of reduced costs
    HVector dual_row;
    dual_row.setup(simplex_lp.numCol_);
    dual_row.clear();
#ifdef HiGHSDEV
    double price_full_historical_density = 1;
    if (simplex_info.analyse_iterations)
      analysis.operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_FULL,
                                     dual_row, price_full_historical_density);
#endif
    matrix.priceByColumn(dual_row, dual_col);
#ifdef HiGHSDEV
    if (simplex_info.analyse_iterations)
      analysis.operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_FULL,
                                    dual_row);
#endif
    for (int i = 0; i < simplex_lp.numCol_; i++)
      simplex_info.workDual_[i] -= dual_row.array[i];
    for (int i = simplex_lp.numCol_; i < numTot; i++)
      simplex_info.workDual_[i] -= dual_col.array[i - simplex_lp.numCol_];
    // Possibly analyse the computed dual values
    debugComputeDual(highs_model_object, debug_previous_workDual,
                     debug_basic_costs, dual_col.array);
  }
  // Now have nonbasic duals
  simplex_lp_status.has_nonbasic_dual_values = true;
}

void correctDual(HighsModelObject& highs_model_object,
                 int* free_infeasibility_count) {
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsRandom& random = highs_model_object.random_;
  const double tau_d =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;
  const double inf = HIGHS_CONST_INF;
  int workCount = 0;
  double flip_dual_objective_value_change = 0;
  double shift_dual_objective_value_change = 0;
  int num_flip = 0;
  int num_shift = 0;
  double sum_flip = 0;
  double sum_shift = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++) {
    if (simplex_basis.nonbasicFlag_[i]) {
      if (simplex_info.workLower_[i] == -inf &&
          simplex_info.workUpper_[i] == inf) {
        // FREE variable
        workCount += (fabs(simplex_info.workDual_[i]) >= tau_d);
      } else if (simplex_basis.nonbasicMove_[i] * simplex_info.workDual_[i] <=
                 -tau_d) {
        if (simplex_info.workLower_[i] != -inf &&
            simplex_info.workUpper_[i] != inf) {
          // Boxed variable = flip
          const int move = simplex_basis.nonbasicMove_[i];
          flip_bound(highs_model_object, i);
          double flip = simplex_info.workUpper_[i] - simplex_info.workLower_[i];
          // Negative dual at lower bound (move=1): flip to upper
          // bound so objective contribution is change in value (flip)
          // times dual, being move*flip*dual
          //
          // Positive dual at upper bound (move=-1): flip to lower
          // bound so objective contribution is change in value
          // (-flip) times dual, being move*flip*dual
          double local_dual_objective_change =
              move * flip * simplex_info.workDual_[i];
          local_dual_objective_change *= highs_model_object.scale_.cost_;
          flip_dual_objective_value_change += local_dual_objective_change;
          num_flip++;
          sum_flip += fabs(flip);
        } else if (simplex_info.allow_cost_perturbation) {
          // Other variable = shift
          //
          // Before 07/01/20, these shifts were always done, but doing
          // it after cost perturbation has been removed can lead to
          // cycling when primal infeasibility has been detecteed in
          // Phase 2, since the shift below removes dual
          // infeasibilities, which are then reinstated after the dual
          // values are recomputed.
          //
          // ToDo: Not shifting leads to dual infeasibilities when an
          // LP is declared to be (primal) infeasible. Should go to
          // phase 1 primal simplex to "prove" infeasibility.
          simplex_info.costs_perturbed = 1;
          std::string direction;
          double shift;
          if (simplex_basis.nonbasicMove_[i] == 1) {
            direction = "  up";
            double dual = (1 + random.fraction()) * tau_d;
            shift = dual - simplex_info.workDual_[i];
            simplex_info.workDual_[i] = dual;
            simplex_info.workCost_[i] = simplex_info.workCost_[i] + shift;
          } else {
            direction = "down";
            double dual = -(1 + random.fraction()) * tau_d;
            shift = dual - simplex_info.workDual_[i];
            simplex_info.workDual_[i] = dual;
            simplex_info.workCost_[i] = simplex_info.workCost_[i] + shift;
          }
          double local_dual_objective_change =
              shift * simplex_info.workValue_[i];
          local_dual_objective_change *= highs_model_object.scale_.cost_;
          shift_dual_objective_value_change += local_dual_objective_change;
          num_shift++;
          sum_shift += fabs(shift);
          HighsPrintMessage(
              highs_model_object.options_.output,
              highs_model_object.options_.message_level, ML_VERBOSE,
              "Move %s: cost shift = %g; objective change = %g\n",
              direction.c_str(), shift, local_dual_objective_change);
        }
      }
    }
  }
  if (num_flip)
    HighsPrintMessage(
        highs_model_object.options_.output,
        highs_model_object.options_.message_level, ML_VERBOSE,
        "Performed %d flip(s): total = %g; objective change = %g\n", num_flip,
        sum_flip, flip_dual_objective_value_change);
  if (num_shift)
    HighsPrintMessage(
        highs_model_object.options_.output,
        highs_model_object.options_.message_level, ML_DETAILED,
        "Performed %d cost shift(s): total = %g; objective change = %g\n",
        num_shift, sum_shift, shift_dual_objective_value_change);
  *free_infeasibility_count = workCount;
}

// Record the shift in the cost of a particular column
void shift_cost(HighsModelObject& highs_model_object, int iCol, double amount) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  simplex_info.costs_perturbed = 1;
  assert(simplex_info.workShift_[iCol] == 0);
  simplex_info.workShift_[iCol] = amount;
}

// Undo the shift in the cost of a particular column
void shift_back(HighsModelObject& highs_model_object, int iCol) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  simplex_info.workDual_[iCol] -= simplex_info.workShift_[iCol];
  /*
  if (simplex_info.workShift_[iCol]) {
    printf("shift_back: column %d: shift = %g; value = %g\n",
           iCol, simplex_info.workShift_[iCol],
           simplex_info.workValue_[iCol]);
    simplex_info.updated_dual_objective_value -=
      simplex_info.workShift_[iCol] * simplex_info.workValue_[iCol];
  }
  */
  simplex_info.workShift_[iCol] = 0;
}

// The major model updates. Factor calls factor.update; Matrix
// calls matrix.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void update_factor(HighsModelObject& highs_model_object, HVector* column,
                   HVector* row_ep, int* iRow, int* hint) {
  //    HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  HFactor& factor = highs_model_object.factor_;
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;

  analysis.simplexTimerStart(UpdateFactorClock);
  factor.update(column, row_ep, iRow, hint);
  // Now have a representation of B^{-1}, but it is not fresh
  simplex_lp_status.has_invert = true;
  if (simplex_info.update_count >= simplex_info.update_limit)
    *hint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  analysis.simplexTimerStop(UpdateFactorClock);
}

void update_pivots(HighsModelObject& highs_model_object, int columnIn,
                   int rowOut, int sourceOut) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;

  analysis.simplexTimerStart(UpdatePivotsClock);
  int columnOut = simplex_basis.basicIndex_[rowOut];

  // Incoming variable
  simplex_basis.basicIndex_[rowOut] = columnIn;
  simplex_basis.nonbasicFlag_[columnIn] = 0;
  simplex_basis.nonbasicMove_[columnIn] = 0;
  simplex_info.baseLower_[rowOut] = simplex_info.workLower_[columnIn];
  simplex_info.baseUpper_[rowOut] = simplex_info.workUpper_[columnIn];

  // Outgoing variable
  simplex_basis.nonbasicFlag_[columnOut] = 1;
  if (simplex_info.workLower_[columnOut] ==
      simplex_info.workUpper_[columnOut]) {
    simplex_info.workValue_[columnOut] = simplex_info.workLower_[columnOut];
    simplex_basis.nonbasicMove_[columnOut] = 0;
  } else if (sourceOut == -1) {
    simplex_info.workValue_[columnOut] = simplex_info.workLower_[columnOut];
    simplex_basis.nonbasicMove_[columnOut] = 1;
  } else {
    simplex_info.workValue_[columnOut] = simplex_info.workUpper_[columnOut];
    simplex_basis.nonbasicMove_[columnOut] = -1;
  }
  // Update the dual objective value
  double nwValue = simplex_info.workValue_[columnOut];
  double vrDual = simplex_info.workDual_[columnOut];
  double dl_dual_objective_value = nwValue * vrDual;
  simplex_info.updated_dual_objective_value += dl_dual_objective_value;
  simplex_info.update_count++;
  // Update the number of basic logicals
  if (columnOut < simplex_lp.numCol_) simplex_info.num_basic_logicals -= 1;
  if (columnIn < simplex_lp.numCol_) simplex_info.num_basic_logicals += 1;
  // No longer have a representation of B^{-1}, and certainly not
  // fresh!
  simplex_lp_status.has_invert = false;
  simplex_lp_status.has_fresh_invert = false;
  // Data are no longer fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = false;
  analysis.simplexTimerStop(UpdatePivotsClock);
}

void update_matrix(HighsModelObject& highs_model_object, int columnIn,
                   int columnOut) {
  HMatrix& matrix = highs_model_object.matrix_;
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;

  analysis.simplexTimerStart(UpdateMatrixClock);
  matrix.update(columnIn, columnOut);
  analysis.simplexTimerStop(UpdateMatrixClock);
}

bool reinvertOnNumericalTrouble(const std::string method_name,
                                const HighsModelObject& highs_model_object,
                                double& numerical_trouble_measure,
                                const double alpha_from_col,
                                const double alpha_from_row,
                                const double numerical_trouble_tolerance) {
  double abs_alpha_from_col = fabs(alpha_from_col);
  double abs_alpha_from_row = fabs(alpha_from_row);
  double min_abs_alpha = min(abs_alpha_from_col, abs_alpha_from_row);
  double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  numerical_trouble_measure = abs_alpha_diff / min_abs_alpha;
  const int update_count = highs_model_object.simplex_info_.update_count;
  // Reinvert if the relative difference is large enough, and updates have been
  // performed
  const bool reinvert =
      numerical_trouble_measure > numerical_trouble_tolerance &&
      update_count > 0;
#ifdef HiGHSDEV
  const int iteration_count = highs_model_object.iteration_counts_.simplex;
  string model_name = highs_model_object.simplex_lp_.model_name_;
  const bool rp_numerical_trouble = false;  // true;//
  if (rp_numerical_trouble)
    printf("%s Measure %11.4g from [Col: %11.4g; Row: %11.4g; Diff = %11.4g]\n",
           method_name.c_str(), numerical_trouble_measure, abs_alpha_from_col,
           abs_alpha_from_row, abs_alpha_diff);
#endif
  if (reinvert) {
#ifdef HiGHSDEV
    printf(
        "%s has identified numerical trouble solving LP %s in iteration %d so "
        "reinvert\n",
        method_name.c_str(), model_name.c_str(), iteration_count);
#else
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "HiGHS has identified numerical trouble so reinvert");
#endif
  }
  return reinvert;
}

// If the scaled LP's model status is optimal, gets suggested
// feasibility tolerances for resolving the scaled LP. Assumes that
// the unscaled primal and dual infeasibilities are not known in
// highs_model_object, so also passes get_unscaled_solution_params as
// values to check from
HighsStatus getNewInfeasibilityTolerancesFromSimplexBasicSolution(
    const HighsModelObject& highs_model_object,
    HighsSolutionParams& get_unscaled_solution_params,
    double& new_scaled_primal_feasibility_tolerance,
    double& new_scaled_dual_feasibility_tolerance) {
  HighsSolutionParams get_scaled_solution_params =
      highs_model_object.scaled_solution_params_;
  return getInfeasibilitiesAndNewTolerances(
      highs_model_object.options_, highs_model_object.lp_,
      highs_model_object.scale_, highs_model_object.simplex_basis_,
      highs_model_object.simplex_info_, highs_model_object.scaled_model_status_,
      get_unscaled_solution_params, highs_model_object.scaled_solution_params_,
      get_unscaled_solution_params, get_scaled_solution_params,
      new_scaled_primal_feasibility_tolerance,
      new_scaled_dual_feasibility_tolerance);
}

// Gets the unscaled and scaled primal and dual infeasibilities from a
// simplex basic solution. The values in unscaled_solution_params and
// scaled_solution_params are checked against them. If the scaled LP's
// model status is optimal, gets suggested feasibility tolerances for
// resolving the scaled LP
HighsStatus getInfeasibilitiesAndNewTolerances(
    const HighsOptions& options, const HighsLp& lp, const HighsScale& scale,
    const SimplexBasis& basis, const HighsSimplexInfo& simplex_info,
    const HighsModelStatus scaled_model_status,
    const HighsSolutionParams& unscaled_solution_params,
    const HighsSolutionParams& scaled_solution_params,
    HighsSolutionParams& get_unscaled_solution_params,
    HighsSolutionParams& get_scaled_solution_params,
    double& new_scaled_primal_feasibility_tolerance,
    double& new_scaled_dual_feasibility_tolerance) {
  const double unscaled_primal_feasibility_tolerance =
      unscaled_solution_params.primal_feasibility_tolerance;
  const double unscaled_dual_feasibility_tolerance =
      unscaled_solution_params.dual_feasibility_tolerance;

  get_unscaled_solution_params = unscaled_solution_params;
  get_scaled_solution_params = scaled_solution_params;

  int& num_unscaled_primal_infeasibilities =
      get_unscaled_solution_params.num_primal_infeasibilities;
  double& max_unscaled_primal_infeasibility =
      get_unscaled_solution_params.max_primal_infeasibility;
  double& sum_unscaled_primal_infeasibilities =
      get_unscaled_solution_params.sum_primal_infeasibilities;
  int& num_unscaled_dual_infeasibilities =
      get_unscaled_solution_params.num_dual_infeasibilities;
  double& max_unscaled_dual_infeasibility =
      get_unscaled_solution_params.max_dual_infeasibility;
  double& sum_unscaled_dual_infeasibilities =
      get_unscaled_solution_params.sum_dual_infeasibilities;

  int& num_scaled_primal_infeasibilities =
      get_scaled_solution_params.num_primal_infeasibilities;
  double& max_scaled_primal_infeasibility =
      get_scaled_solution_params.max_primal_infeasibility;
  double& sum_scaled_primal_infeasibilities =
      get_scaled_solution_params.sum_primal_infeasibilities;
  int& num_scaled_dual_infeasibilities =
      get_scaled_solution_params.num_dual_infeasibilities;
  double& max_scaled_dual_infeasibility =
      get_scaled_solution_params.max_dual_infeasibility;
  double& sum_scaled_dual_infeasibilities =
      get_scaled_solution_params.sum_dual_infeasibilities;

  // Invalidate the unscaled and scaled infeasibility params
  invalidateSolutionInfeasibilityParams(get_unscaled_solution_params);
  invalidateSolutionInfeasibilityParams(get_scaled_solution_params);
  // Zero the counts of unscaled and scaled primal and dual infeasibilities
  num_unscaled_primal_infeasibilities = 0;
  num_unscaled_dual_infeasibilities = 0;
  num_scaled_primal_infeasibilities = 0;
  num_scaled_dual_infeasibilities = 0;

  // If the scaled LP has beeen solved to optimality, look at the
  // scaled solution and, if there are infeasibilities, identify new
  // feasibility tolerances for the scaled LP
  const bool get_new_scaled_feasibility_tolerances =
      scaled_model_status == HighsModelStatus::OPTIMAL;
  // The scaled infeasibility parameters are not known if the dual
  // objective upper bound has been reached, the time limit has been
  // reached, or the iteration limit has been reached,
  const bool check_scaled_solution_params =
      scaled_model_status !=
          HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND &&
      scaled_model_status != HighsModelStatus::REACHED_TIME_LIMIT &&
      scaled_model_status != HighsModelStatus::REACHED_ITERATION_LIMIT;

  const double scaled_primal_feasibility_tolerance =
      scaled_solution_params.primal_feasibility_tolerance;
  const double scaled_dual_feasibility_tolerance =
      scaled_solution_params.dual_feasibility_tolerance;

  if (get_new_scaled_feasibility_tolerances) {
    new_scaled_primal_feasibility_tolerance =
        scaled_primal_feasibility_tolerance;
    new_scaled_dual_feasibility_tolerance = scaled_dual_feasibility_tolerance;
  }
  for (int iVar = 0; iVar < lp.numCol_ + lp.numRow_; iVar++) {
    // Look at the dual infeasibilities of nonbasic variables
    if (basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE) continue;
    // No dual infeasiblity for fixed rows and columns
    if (simplex_info.workLower_[iVar] == simplex_info.workUpper_[iVar])
      continue;
    bool col = iVar < lp.numCol_;
    double scale_mu;
    int iCol = 0;
    int iRow = 0;
    if (col) {
      iCol = iVar;
      scale_mu = 1 / (scale.col_[iCol] / scale.cost_);
    } else {
      iRow = iVar - lp.numCol_;
      scale_mu = scale.row_[iRow] * scale.cost_;
    }
    const double scaled_dual = simplex_info.workDual_[iVar];
    const double unscaled_dual = scaled_dual * scale_mu;
    const double lower = simplex_info.workLower_[iVar];
    const double upper = simplex_info.workUpper_[iVar];

    double scaled_dual_infeasibility;
    double unscaled_dual_infeasibility;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      scaled_dual_infeasibility = fabs(scaled_dual);
      unscaled_dual_infeasibility = fabs(unscaled_dual);
    } else {
      // Not fixed: any dual infeasibility is given by value signed by
      // nonbasicMove. This assumes that nonbasicMove=0 for fixed
      // variables
      scaled_dual_infeasibility = -basis.nonbasicMove_[iVar] * scaled_dual;
      unscaled_dual_infeasibility = -basis.nonbasicMove_[iVar] * unscaled_dual;
    }
    if (scaled_dual_infeasibility > 0) {
      if (scaled_dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_scaled_dual_infeasibilities++;
      max_scaled_dual_infeasibility =
          max(scaled_dual_infeasibility, max_scaled_dual_infeasibility);
      sum_scaled_dual_infeasibilities += scaled_dual_infeasibility;
    }
    if (unscaled_dual_infeasibility > 0) {
      if (unscaled_dual_infeasibility >= unscaled_dual_feasibility_tolerance) {
        num_unscaled_dual_infeasibilities++;
        if (get_new_scaled_feasibility_tolerances) {
          double multiplier = unscaled_dual_feasibility_tolerance / scale_mu;
#ifdef HiGHSDEV
          /*
          double value = simplex_info.workValue_[iVar];
          HighsLogMessage(logfile, HighsMessageType::INFO,
                          "Var %6d (%6d, %6d): [%11.4g, %11.4g, %11.4g] %11.4g
          s=%11.4g %11.4g: Mu = %g", iVar, iCol, iRow, lower, value, upper,
                          scaled_dual_infeasibility, scale_mu,
          unscaled_dual_infeasibility, multiplier);
          */
#endif
          new_scaled_dual_feasibility_tolerance =
              min(multiplier, new_scaled_dual_feasibility_tolerance);
        }
      }
      max_unscaled_dual_infeasibility =
          max(unscaled_dual_infeasibility, max_unscaled_dual_infeasibility);
      sum_unscaled_dual_infeasibilities += unscaled_dual_infeasibility;
    }
  }
  // Look at the primal infeasibilities of basic variables
  for (int ix = 0; ix < lp.numRow_; ix++) {
    int iVar = basis.basicIndex_[ix];
    bool col = iVar < lp.numCol_;
    double scale_mu;
    int iCol = 0;
    int iRow = 0;
    if (col) {
      iCol = iVar;
      scale_mu = scale.col_[iCol];
    } else {
      iRow = iVar - lp.numCol_;
      scale_mu = 1 / scale.row_[iRow];
    }
    // Look at the basic primal infeasibilities

    double lower = simplex_info.baseLower_[ix];
    double upper = simplex_info.baseUpper_[ix];
    double value = simplex_info.baseValue_[ix];
    double scaled_primal_infeasibility =
        max(max(lower - value, value - upper), 0.);
    double unscaled_primal_infeasibility =
        scaled_primal_infeasibility * scale_mu;
    if (scaled_primal_infeasibility > scaled_primal_feasibility_tolerance) {
      num_scaled_primal_infeasibilities++;
    }
    max_scaled_primal_infeasibility =
        max(scaled_primal_infeasibility, max_scaled_primal_infeasibility);
    sum_scaled_primal_infeasibilities += scaled_primal_infeasibility;
    if (unscaled_primal_infeasibility > unscaled_primal_feasibility_tolerance) {
      num_unscaled_primal_infeasibilities++;
      if (get_new_scaled_feasibility_tolerances) {
        double multiplier = unscaled_primal_feasibility_tolerance / scale_mu;
#ifdef HiGHSDEV
        /*
         HighsLogMessage(logfile, HighsMessageType::INFO,
                        "Var %6d (%6d, %6d): [%11.4g, %11.4g, %11.4g] %11.4g
         s=%11.4g %11.4g: Mu = %g", iVar, iCol, iRow, lower, value, upper,
                        scaled_primal_infeasibility, scale_mu,
         unscaled_primal_infeasibility, multiplier);
        */
#endif
        new_scaled_primal_feasibility_tolerance =
            min(multiplier, new_scaled_primal_feasibility_tolerance);
      }
    }
    max_unscaled_primal_infeasibility =
        max(unscaled_primal_infeasibility, max_unscaled_primal_infeasibility);
    sum_unscaled_primal_infeasibilities += unscaled_primal_infeasibility;
  }

  HighsDebugStatus debug_status;
  debug_status = debugCompareSolutionInfeasibilityParams(
      options, get_unscaled_solution_params, unscaled_solution_params);
  if (debug_status != HighsDebugStatus::OK) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Unequal unscaled solution infeasibility params in "
                    "getPrimalDualInfeasibilitiesFromSimplexBasicSolution");
    assert(debug_status == HighsDebugStatus::OK);
    return HighsStatus::Error;
  }
  if (check_scaled_solution_params) {
    debug_status = debugCompareSolutionInfeasibilityParams(
        options, get_scaled_solution_params, scaled_solution_params);
    if (debug_status != HighsDebugStatus::OK) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Unequal scaled solution infeasibility params in "
                      "getPrimalDualInfeasibilitiesFromSimplexBasicSolution");
      assert(debug_status == HighsDebugStatus::OK);
      return HighsStatus::Error;
    }
  }
  return HighsStatus::OK;
}

void logRebuild(HighsModelObject& highs_model_object, const bool primal,
                const int solve_phase) {
  HighsSolutionParams& scaled_solution_params =
      highs_model_object.scaled_solution_params_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  double objective_value;
  string simplex_variant;
  if (primal) {
    simplex_variant = "Pr";
    objective_value = simplex_info.primal_objective_value;
  } else {
    simplex_variant = "Du";
    objective_value = simplex_info.dual_objective_value;
  }
  if (solve_phase < 2) {
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Iter %10d: %20.10e %sPh%1d",
                    highs_model_object.iteration_counts_.simplex,
                    objective_value, simplex_variant.c_str(), solve_phase);
  } else if (!primal && scaled_solution_params.sum_dual_infeasibilities == 0) {
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Iter %10d: %20.10e %sPh%1d Pr: %d(%g)",
                    highs_model_object.iteration_counts_.simplex,
                    objective_value, simplex_variant.c_str(), solve_phase,
                    scaled_solution_params.num_primal_infeasibilities,
                    scaled_solution_params.sum_primal_infeasibilities);
  } else if (primal && scaled_solution_params.num_primal_infeasibilities) {
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Iter %10d: %20.10e %sPh%1d Pr: %d(%g); Du: %d(%g)",
                    highs_model_object.iteration_counts_.simplex,
                    objective_value, simplex_variant.c_str(), 1,
                    scaled_solution_params.num_primal_infeasibilities,
                    scaled_solution_params.sum_primal_infeasibilities,
                    scaled_solution_params.num_dual_infeasibilities,
                    scaled_solution_params.sum_dual_infeasibilities);
  } else {
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Iter %10d: %20.10e %sPh%1d Pr: %d(%g); Du: %d(%g)",
                    highs_model_object.iteration_counts_.simplex,
                    objective_value, simplex_variant.c_str(), solve_phase,
                    scaled_solution_params.num_primal_infeasibilities,
                    scaled_solution_params.sum_primal_infeasibilities,
                    scaled_solution_params.num_dual_infeasibilities,
                    scaled_solution_params.sum_dual_infeasibilities);
  }
}

void reportSimplexLpStatus(HighsSimplexLpStatus& simplex_lp_status,
                           const char* message) {
  printf("\nReporting solver status and flags: %s\n\n", message);
  printf("  valid =                          %d\n", simplex_lp_status.valid);
  printf("  is_dualised =                    %d\n",
         simplex_lp_status.is_dualised);
  printf("  is_permuted =                    %d\n",
         simplex_lp_status.is_permuted);
  printf("  is_scaled =                      %d\n",
         simplex_lp_status.scaling_tried);
  printf("  has_basis =                      %d\n",
         simplex_lp_status.has_basis);
  printf("  has_matrix_col_wise =            %d\n",
         simplex_lp_status.has_matrix_col_wise);
  printf("  has_matrix_row_wise =            %d\n",
         simplex_lp_status.has_matrix_row_wise);
  printf("  has_factor_arrays =              %d\n",
         simplex_lp_status.has_factor_arrays);
  printf("  has_dual_steepest_edge_weights = %d\n",
         simplex_lp_status.has_dual_steepest_edge_weights);
  printf("  has_nonbasic_dual_values =       %d\n",
         simplex_lp_status.has_nonbasic_dual_values);
  printf("  has_basic_primal_values =        %d\n",
         simplex_lp_status.has_basic_primal_values);
  printf("  has_invert =                     %d\n",
         simplex_lp_status.has_invert);
  printf("  has_fresh_invert =               %d\n",
         simplex_lp_status.has_fresh_invert);
  printf("  has_fresh_rebuild =              %d\n",
         simplex_lp_status.has_fresh_rebuild);
  printf("  has_dual_objective_value =       %d\n",
         simplex_lp_status.has_dual_objective_value);
  printf("  has_primal_objective_value =     %d\n",
         simplex_lp_status.has_primal_objective_value);
}

void invalidateSimplexLpBasisArtifacts(
    HighsSimplexLpStatus& simplex_lp_status) {
  // Invalidate the artifacts of the basis of the simplex LP
  simplex_lp_status.has_matrix_col_wise = false;
  simplex_lp_status.has_matrix_row_wise = false;
  simplex_lp_status.has_factor_arrays = false;
  simplex_lp_status.has_dual_steepest_edge_weights = false;
  simplex_lp_status.has_nonbasic_dual_values = false;
  simplex_lp_status.has_basic_primal_values = false;
  simplex_lp_status.has_invert = false;
  simplex_lp_status.has_fresh_invert = false;
  simplex_lp_status.has_fresh_rebuild = false;
  simplex_lp_status.has_dual_objective_value = false;
  simplex_lp_status.has_primal_objective_value = false;
}

void invalidateSimplexLpBasis(HighsSimplexLpStatus& simplex_lp_status) {
  // Invalidate the basis of the simplex LP, and all its other
  // properties - since they are basis-related
  simplex_lp_status.has_basis = false;
  invalidateSimplexLpBasisArtifacts(simplex_lp_status);
}

void invalidateSimplexLp(HighsSimplexLpStatus& simplex_lp_status) {
  simplex_lp_status.valid = false;
  simplex_lp_status.is_dualised = false;
  simplex_lp_status.is_permuted = false;
  simplex_lp_status.scaling_tried = false;
  invalidateSimplexLpBasis(simplex_lp_status);
}

void updateSimplexLpStatus(HighsSimplexLpStatus& simplex_lp_status,
                           LpAction action) {
  switch (action) {
    case LpAction::DUALISE:
#ifdef HIGHSDEV
      printf(" LpAction::DUALISE\n");
#endif
      simplex_lp_status.is_dualised = true;
      invalidateSimplexLpBasis(simplex_lp_status);
      break;
    case LpAction::PERMUTE:
#ifdef HIGHSDEV
      printf(" LpAction::PERMUTE\n");
#endif
      simplex_lp_status.is_permuted = true;
      invalidateSimplexLpBasis(simplex_lp_status);
      break;
    case LpAction::SCALE:
#ifdef HIGHSDEV
      printf(" LpAction::SCALE\n");
#endif
      simplex_lp_status.scaling_tried = true;
      invalidateSimplexLpBasis(simplex_lp_status);
      break;
    case LpAction::NEW_COSTS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_COSTS\n");
#endif
      simplex_lp_status.has_nonbasic_dual_values = false;
      simplex_lp_status.has_fresh_rebuild = false;
      simplex_lp_status.has_dual_objective_value = false;
      simplex_lp_status.has_primal_objective_value = false;
      break;
    case LpAction::NEW_BOUNDS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_BOUNDS\n");
#endif
      simplex_lp_status.has_basic_primal_values = false;
      simplex_lp_status.has_fresh_rebuild = false;
      simplex_lp_status.has_dual_objective_value = false;
      simplex_lp_status.has_primal_objective_value = false;
      break;
    case LpAction::NEW_BASIS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_BASIS\n");
#endif
      invalidateSimplexLpBasis(simplex_lp_status);
      break;
    case LpAction::NEW_COLS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_COLS\n");
#endif
      invalidateSimplexLpBasisArtifacts(simplex_lp_status);
      break;
    case LpAction::NEW_ROWS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_ROWS\n");
#endif
      invalidateSimplexLpBasisArtifacts(simplex_lp_status);
      break;
    case LpAction::DEL_COLS:
#ifdef HIGHSDEV
      printf(" LpAction::DEL_COLS\n");
#endif
      invalidateSimplexLpBasis(simplex_lp_status);
      break;
    case LpAction::DEL_ROWS:
#ifdef HIGHSDEV
      printf(" LpAction::DEL_ROWS\n");
#endif
      invalidateSimplexLpBasis(simplex_lp_status);
      break;
    case LpAction::DEL_ROWS_BASIS_OK:
#ifdef HIGHSDEV
      printf(" LpAction::DEL_ROWS_BASIS_OK\n");
#endif
      //      simplex_info.simplex_lp_ = true;
      break;
    default:
#ifdef HIGHSDEV
      printf(" Unrecognised LpAction::%d\n", (int)action);
#endif
      break;
  }
}

bool isBasisRightSize(const HighsLp& lp, const SimplexBasis& basis) {
  bool right_size = true;
  right_size =
      (int)basis.nonbasicFlag_.size() == lp.numCol_ + lp.numRow_ && right_size;
  right_size =
      (int)basis.nonbasicMove_.size() == lp.numCol_ + lp.numRow_ && right_size;
  right_size = (int)basis.basicIndex_.size() == lp.numRow_ && right_size;
  return right_size;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolverData.h"

#include <random>
#include <sstream>

#include "../extern/pdqsort/pdqsort.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsRedcostFixing.h"
#include "mip/MipTimer.h"
#include "parallel/HighsParallel.h"
#include "presolve/HPresolve.h"
#include "util/HighsIntegers.h"

HighsMipSolverData::HighsMipSolverData(HighsMipSolver& mipsolver)
    : mipsolver(mipsolver),
      pseudocosts(1),
      parallel_lock(false),
      heuristics(mipsolver),
      cliquetable(mipsolver.numCol()),
      implications(mipsolver),
      objectiveFunction(mipsolver),
      presolve_status(HighsPresolveStatus::kNotSet),
      cliquesExtracted(false),
      rowMatrixSet(false),
      analyticCenterComputed(false),
      analyticCenterStatus(HighsModelStatus::kNotset),
      detectSymmetries(false),
      numRestarts(0),
      numRestartsRoot(0),
      numCliqueEntriesAfterPresolve(0),
      numCliqueEntriesAfterFirstPresolve(0),
      feastol(0.0),
      epsilon(0.0),
      heuristic_effort(0.0),
      dispfreq(0),
      firstlpsolobj(-kHighsInf),
      rootlpsolobj(-kHighsInf),
      numintegercols(0),
      maxTreeSizeLog2(0),
      pruned_treeweight(0),
      avgrootlpiters(0.0),
      disptime(0.0),
      last_disptime(0.0),
      firstrootlpiters(0),
      num_nodes(0),
      num_leaves(0),
      num_leaves_before_run(0),
      num_nodes_before_run(0),
      total_repair_lp(0),
      total_repair_lp_feasible(0),
      total_repair_lp_iterations(0),
      total_lp_iterations(0),
      heuristic_lp_iterations(0),
      sepa_lp_iterations(0),
      sb_lp_iterations(0),
      total_lp_iterations_before_run(0),
      heuristic_lp_iterations_before_run(0),
      sepa_lp_iterations_before_run(0),
      sb_lp_iterations_before_run(0),
      num_disp_lines(0),
      numImprovingSols(0),
      lower_bound(-kHighsInf),
      upper_bound(kHighsInf),
      upper_limit(kHighsInf),
      optimality_limit(kHighsInf),
      debugSolution(mipsolver) {
  lps.emplace_back(mipsolver);
  domains.emplace_back(mipsolver);
  conflictpools.emplace_back(5 * mipsolver.options_mip_->mip_pool_age_limit,
                             mipsolver.options_mip_->mip_pool_soft_limit);
  cutpools.emplace_back(mipsolver.numCol(),
                        mipsolver.options_mip_->mip_pool_age_limit,
                        mipsolver.options_mip_->mip_pool_soft_limit, 0);
  getDomain().addCutpool(getCutPool());
  getDomain().addConflictPool(getConflictPool());
  cliquetable.setAllowParallel(!mipsolver.submip);
  worker_lp_iterations_stop.store(std::numeric_limits<int64_t>::max(),
                                  std::memory_order_relaxed);
}

std::string HighsMipSolverData::solutionSourceToString(
    const int solution_source, const bool code) const {
  if (solution_source == kSolutionSourceNone) {
    if (code) return " ";
    return "None";
    //  } else if (solution_source == kSolutionSourceInitial) {
    //    if (code) return "0";
    //    return "Initial";
  } else if (solution_source == kSolutionSourceBranching) {
    if (code) return "B";
    return "Branching";
  } else if (solution_source == kSolutionSourceCentralRounding) {
    if (code) return "C";
    return "Central rounding";
  } else if (solution_source == kSolutionSourceFeasibilityPump) {
    if (code) return "F";
    return "Feasibility pump";
  } else if (solution_source == kSolutionSourceHeuristic) {
    if (code) return "H";
    return "Heuristic";
  } else if (solution_source == kSolutionSourceShifting) {
    if (code) return "I";
    return "Shifting";
  } else if (solution_source == kSolutionSourceFeasibilityJump) {
    if (code) return "J";
    return "Feasibility jump";
  } else if (solution_source == kSolutionSourceSubMip) {
    if (code) return "L";
    return "Sub-MIP";
  } else if (solution_source == kSolutionSourceEmptyMip) {
    if (code) return "P";
    return "Empty MIP";
  } else if (solution_source == kSolutionSourceRandomizedRounding) {
    if (code) return "R";
    return "Randomized rounding";
  } else if (solution_source == kSolutionSourceSolveLp) {
    if (code) return "S";
    return "Solve LP";
  } else if (solution_source == kSolutionSourceEvaluateNode) {
    if (code) return "T";
    return "Evaluate node";
  } else if (solution_source == kSolutionSourceUnbounded) {
    if (code) return "U";
    return "Unbounded";
  } else if (solution_source == kSolutionSourceUserSolution) {
    if (code) return "X";
    return "User solution";
  } else if (solution_source == kSolutionSourceHighsSolution) {
    if (code) return "Y";
    return "HiGHS solution";
  } else if (solution_source == kSolutionSourceZiRound) {
    if (code) return "Z";
    return "ZI Round";
  } else if (solution_source == kSolutionSourceTrivialZ) {
    if (code) return "z";
    return "Trivial zero";
  } else if (solution_source == kSolutionSourceTrivialL) {
    if (code) return "l";
    return "Trivial lower";
  } else if (solution_source == kSolutionSourceTrivialU) {
    if (code) return "u";
    return "Trivial upper";
  } else if (solution_source == kSolutionSourceTrivialP) {
    if (code) return "p";
    return "Trivial point";
  } else if (solution_source == kSolutionSourceCleanup) {
    if (code) return " ";
    return "";
  } else {
    printf("HighsMipSolverData::solutionSourceToString: Unknown source = %d\n",
           solution_source);
    assert(0 == 111);
    if (code) return "*";
    return "None";
  }
}

bool HighsMipSolverData::checkSolution(
    const std::vector<double>& solution) const {
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    if (solution[i] < mipsolver.model_->col_lower_[i] - feastol) return false;
    if (solution[i] > mipsolver.model_->col_upper_[i] + feastol) return false;
    if (mipsolver.isColInteger(i) && fractionality(solution[i]) > feastol)
      return false;
  }

  for (HighsInt i = 0; i != mipsolver.numRow(); ++i) {
    double rowactivity = 0.0;

    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];

    for (HighsInt j = start; j != end; ++j)
      rowactivity += solution[ARindex_[j]] * ARvalue_[j];

    if (rowactivity > mipsolver.rowUpper(i) + feastol) return false;
    if (rowactivity < mipsolver.rowLower(i) - feastol) return false;
  }

  return true;
}

std::vector<std::tuple<HighsInt, HighsInt, double>>
HighsMipSolverData::getInfeasibleRows(
    const std::vector<double>& solution) const {
  std::vector<std::tuple<HighsInt, HighsInt, double>> infeasibleRows;
  for (HighsInt i = 0; i != mipsolver.numRow(); ++i) {
    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];

    HighsCDouble row_activity_quad = 0.0;
    for (HighsInt j = start; j != end; ++j)
      row_activity_quad +=
          static_cast<HighsCDouble>(solution[ARindex_[j]]) * ARvalue_[j];

    double row_activity = static_cast<double>(row_activity_quad);
    if (row_activity > mipsolver.rowUpper(i) + feastol) {
      double difference = std::abs(row_activity - mipsolver.rowUpper(i));
      infeasibleRows.push_back({i, +1, difference});
    }
    if (row_activity < mipsolver.rowLower(i) - feastol) {
      double difference = std::abs(mipsolver.rowLower(i) - row_activity);
      infeasibleRows.push_back({i, -1, difference});
    }
  }
  return infeasibleRows;
}

bool HighsMipSolverData::trySolution(const std::vector<double>& solution,
                                     const int solution_source) {
  if (int(solution.size()) != mipsolver.numCol()) return false;

  HighsCDouble obj = 0;

  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    if (solution[i] < mipsolver.model_->col_lower_[i] - feastol) return false;
    if (solution[i] > mipsolver.model_->col_upper_[i] + feastol) return false;
    if (mipsolver.isColInteger(i) && fractionality(solution[i]) > feastol)
      return false;

    obj += static_cast<HighsCDouble>(mipsolver.colCost(i)) * solution[i];
  }

  for (HighsInt i = 0; i != mipsolver.numRow(); ++i) {
    double rowactivity = 0.0;

    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];

    for (HighsInt j = start; j != end; ++j)
      rowactivity += solution[ARindex_[j]] * ARvalue_[j];

    if (rowactivity > mipsolver.rowUpper(i) + feastol) return false;
    if (rowactivity < mipsolver.rowLower(i) - feastol) return false;
  }

  return addIncumbent(solution, double(obj), solution_source);
}

bool HighsMipSolverData::solutionRowFeasible(
    const std::vector<double>& solution) const {
  for (HighsInt i = 0; i != mipsolver.numRow(); ++i) {
    HighsCDouble c_double_rowactivity = HighsCDouble(0.0);

    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];

    for (HighsInt j = start; j != end; ++j)
      c_double_rowactivity +=
          static_cast<HighsCDouble>(solution[ARindex_[j]]) * ARvalue_[j];

    double rowactivity = double(c_double_rowactivity);
    if (rowactivity > mipsolver.rowUpper(i) + feastol) return false;
    if (rowactivity < mipsolver.rowLower(i) - feastol) return false;
  }
  return true;
}

HighsModelStatus HighsMipSolverData::trivialHeuristics() {
  //  printf("\nHighsMipSolverData::trivialHeuristics() Number of continuous
  //  columns is %d\n",
  //	 int(continuous_cols.size()));
  if (continuous_cols.size() > 0) return HighsModelStatus::kNotset;
  const HighsInt num_try_heuristic = 4;
  const std::vector<int> heuristic_source = {
      kSolutionSourceTrivialZ, kSolutionSourceTrivialL, kSolutionSourceTrivialU,
      kSolutionSourceTrivialP};

  std::vector<double> col_lower = mipsolver.model_->col_lower_;
  std::vector<double> col_upper = mipsolver.model_->col_upper_;
  const std::vector<double>& row_lower = mipsolver.model_->row_lower_;
  const std::vector<double>& row_upper = mipsolver.model_->row_upper_;
  const HighsSparseMatrix& matrix = mipsolver.model_->a_matrix_;
  // Determine the following properties, according to which some
  // trivial heuristics are duplicated or fail immediately
  bool all_integer_lower_non_positive = true;
  bool all_integer_lower_zero = true;
  bool all_integer_lower_finite = true;
  bool all_integer_upper_finite = true;
  for (HighsInt integer_col = 0; integer_col < numintegercols; integer_col++) {
    HighsInt iCol = integer_cols[integer_col];
    // Round bounds in to nearest integer
    col_lower[iCol] = std::ceil(col_lower[iCol]);
    col_upper[iCol] = std::floor(col_upper[iCol]);
    const bool legal_bounds =
        col_lower[iCol] <= col_upper[iCol] && col_lower[iCol] < kHighsInf &&
        col_upper[iCol] > -kHighsInf && !std::isnan(col_lower[iCol]) &&
        !std::isnan(col_upper[iCol]);
    if (!legal_bounds) {
      assert(legal_bounds);
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                   "HighsMipSolverData::trivialHeuristics() has detected "
                   "infeasible/illegal bounds [%g, %g] for column %d: MIP is "
                   "infeasible\n",
                   col_lower[iCol], col_upper[iCol], int(iCol));
      return HighsModelStatus::kInfeasible;
    }
    // If bounds are inconsistent then MIP is infeasible
    if (col_lower[iCol] > col_upper[iCol]) return HighsModelStatus::kInfeasible;

    if (col_lower[iCol] > 0) all_integer_lower_non_positive = false;
    if (col_lower[iCol]) all_integer_lower_zero = false;
    if (col_lower[iCol] <= -kHighsInf) all_integer_lower_finite = false;
    if (col_upper[iCol] >= kHighsInf) all_integer_upper_finite = false;
    // Only continue if one of the properties still holds
    if (!(all_integer_lower_non_positive || all_integer_lower_zero ||
          all_integer_upper_finite))
      break;
  }
  const bool all_integer_boxed =
      all_integer_lower_finite && all_integer_upper_finite;
  //  printf(
  //      "Trying trivial heuristics\n"
  //      "   all_integer_lower_non_positive = %d\n"
  //      "   all_integer_lower_zero = %d\n"
  //      "   all_integer_upper_finite = %d\n"
  //      "   all_integer_boxed = %d\n",
  //      all_integer_lower_non_positive, all_integer_lower_zero,
  //      all_integer_upper_finite, all_integer_boxed);
  const double feasibility_tolerance =
      mipsolver.options_mip_->mip_feasibility_tolerance;
  // Loop through the trivial heuristics
  std::vector<double> solution(mipsolver.numCol());
  for (HighsInt try_heuristic = 0; try_heuristic < num_try_heuristic;
       try_heuristic++) {
    if (try_heuristic == 0) {
      // First heuristic is to see whether all-zero for integer
      // variables is feasible
      //
      // If there is a positive lower bound then the heuristic fails
      if (!all_integer_lower_non_positive) continue;
      // Determine whether a zero row activity is feasible
      bool heuristic_failed = false;
      for (HighsInt iRow = 0; iRow < mipsolver.numRow(); iRow++) {
        if (row_lower[iRow] > feasibility_tolerance ||
            row_upper[iRow] < -feasibility_tolerance) {
          heuristic_failed = true;
          break;
        }
      }
      if (heuristic_failed) continue;
      solution.assign(mipsolver.numCol(), 0);
    } else if (try_heuristic == 1) {
      // Second heuristic is to see whether all-lower for integer
      // variables (if distinct from all-zero) is feasible
      if (all_integer_lower_zero) continue;
      // Trivially feasible for columns
      if (!solutionRowFeasible(col_lower)) continue;
      solution = col_lower;
    } else if (try_heuristic == 2) {
      // Third heuristic is to see whether all-upper for integer
      // variables is feasible
      //
      // If there is an infinite upper bound then the heuristic fails
      if (!all_integer_upper_finite) continue;
      // Trivially feasible for columns
      if (!solutionRowFeasible(col_upper)) continue;
      solution = col_upper;
    } else if (try_heuristic == 3) {
      // Fourth heuristic is to see whether the "lock point" is feasible
      if (!all_integer_boxed) continue;
      for (HighsInt integer_col = 0; integer_col < numintegercols;
           integer_col++) {
        HighsInt iCol = integer_cols[integer_col];
        HighsInt num_positive_values = 0;
        HighsInt num_negative_values = 0;
        for (HighsInt iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol + 1];
             iEl++) {
          if (matrix.value_[iEl] > 0)
            num_positive_values++;
          else
            num_negative_values++;
        }
        solution[iCol] = num_positive_values > num_negative_values
                             ? col_lower[iCol]
                             : col_upper[iCol];
      }
      // Trivially feasible for columns
      if (!solutionRowFeasible(solution)) continue;
    }

    HighsCDouble cdouble_obj = 0.0;
    for (HighsInt iCol = 0; iCol < mipsolver.numCol(); iCol++)
      cdouble_obj +=
          static_cast<HighsCDouble>(mipsolver.colCost(iCol)) * solution[iCol];
    double obj = static_cast<double>(cdouble_obj);
    const double save_upper_bound = upper_bound;
    const bool new_incumbent =
        addIncumbent(solution, obj, heuristic_source[try_heuristic]);
    const bool lc_report = false;
    if (lc_report) {
      printf("Trivial heuristic %d has succeeded: objective = %g",
             int(try_heuristic), obj);
      if (new_incumbent) {
        printf("; upper bound from %g to %g\n", save_upper_bound, upper_bound);
      } else {
        printf("\n");
      }
    }
  }
  return HighsModelStatus::kNotset;
}

void HighsMipSolverData::startAnalyticCenterComputation(
    const highs::parallel::TaskGroup& taskGroup) {
  taskGroup.spawn([&]() {
    // first check if the analytic centre computation should be cancelled, e.g.
    // due to early return in the root node evaluation
    //
    // Highs instantiation
    Highs ipm;
    ipm.setProfiling(mipsolver.profiling_);
    ipm.setOptionValue("output_flag", false);
    const std::vector<double>& sol = ipm.getSolution().col_value;
    // Don't use presolve - because this can lead to postsolve putting
    // integer variables onto bounds. This is not just a "less good"
    // AC. It can have implications leading to erroneous fixing of
    // variables and a suboptimal solution declared as optimal.
    ipm.setOptionValue("presolve", kHighsOffString);
    // Determine the solver
    const std::string mip_ipm_solver = mipsolver.options_mip_->mip_ipm_solver;
    // HiPO does not have centring yet, so use IPX
    const bool use_hipo = false;
    const std::string ipm_solver = use_hipo ? kHipoString : kIpxString;
    ipm.setOptionValue("solver", ipm_solver);
    ipm.setOptionValue("ipm_iteration_limit", 200);
    ipm.setOptionValue("run_crossover", kHighsOffString);
    ipm.setOptionValue("run_centring", true);
    HighsLp lpmodel(*mipsolver.model_);
    lpmodel.col_cost_.assign(lpmodel.num_col_, 0.0);
    lpmodel.integrality_.clear();
    ipm.passModel(std::move(lpmodel));
    const bool dump_ipm_lp = false;
    if (dump_ipm_lp && !mipsolver.submip) {
      const std::string file_name = mipsolver.model_->model_name_ + "_ac.mps";
      printf(
          "HighsMipSolverData::startAnalyticCenterComputation: Calling "
          "ipm.writeModel(%s)\n",
          file_name.c_str());
      ipm.writeModel(file_name);
      fflush(stdout);
      exit(1);
    }
    const bool ipm_logging = false;
    if (ipm_logging) {
      bool output_flag;
      ipm.getOptionValue("output_flag", output_flag);
      assert(output_flag == false);
      (void)output_flag;
      ipm.setOptionValue("output_flag", !mipsolver.submip);
    }
    const HighsInt profiling_clock =
        use_hipo ? kSubSolverHipoAc : kSubSolverIpxAc;
    if (mipsolver.profiling_) mipsolver.profiling_->start(profiling_clock);
    ipm.optimizeLp();
    if (mipsolver.profiling_) mipsolver.profiling_->stop(profiling_clock);
    if (ipm_logging) ipm.setOptionValue("output_flag", false);
    if (use_hipo && mip_ipm_solver == kHighsChooseString &&
        HighsInt(sol.size()) != mipsolver.numCol()) {
      printf(
          "In HighsMipSolverData::startAnalyticCenterComputation HiPO has "
          "failed to get a solution: status = %s Try IPX\n",
          ipm.modelStatusToString(ipm.getModelStatus()).c_str());
      // HiPO has failed to get a solution, so try IPX
      ipm.setOptionValue("solver", kIpxString);
      ipm.optimizeLp();
    }
    if (HighsInt(sol.size()) != mipsolver.numCol()) return;
    analyticCenterStatus = ipm.getModelStatus();
    // The model status only indicates that the LP is solved, which is
    // true for any feasible point since the costs are zero. Fixing
    // variables in finishAnalyticCenterComputation requires the point to
    // be centred, so only accept an optimal status if centring succeeded
    if (ipm.getInfo().centring_status != kCentringStatusSuccess)
      analyticCenterStatus = HighsModelStatus::kUnknown;
    analyticCenter = sol;
  });
}

void HighsMipSolverData::finishAnalyticCenterComputation(
    const highs::parallel::TaskGroup& taskGroup) {
  if (mipsolver.profiling_->mip_) {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - starting  analytic centre synch\n",
                 mipsolver.timer_.read());
    fflush(stdout);
  }
  taskGroup.sync();
  if (mipsolver.profiling_->mip_) {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - completed analytic centre synch\n",
                 mipsolver.timer_.read());
    fflush(stdout);
  }
  analyticCenterComputed = true;
  // Disable fixings if presolve=off as it possible to produce valid analytic
  // centers close to a bound when adding sets of equalities and redundant
  // constraints
  if (analyticCenterStatus == HighsModelStatus::kOptimal &&
      mipsolver.options_mip_->presolve != kHighsOffString) {
    HighsInt nfixed = 0;
    HighsInt nintfixed = 0;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      double boundRange = mipsolver.mipdata_->getDomain().col_upper_[i] -
                          mipsolver.mipdata_->getDomain().col_lower_[i];
      if (boundRange == 0.0) continue;

      const double tolerance =
          (mipsolver.mipdata_->feastol / 10) * std::min(boundRange, 1.0);

      if (analyticCenter[i] <= mipsolver.model_->col_lower_[i] + tolerance) {
        mipsolver.mipdata_->getDomain().changeBound(
            HighsBoundType::kUpper, i, mipsolver.model_->col_lower_[i],
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->getDomain().infeasible()) return;
        ++nfixed;
        if (mipsolver.isColInteger(i)) ++nintfixed;
      } else if (analyticCenter[i] >=
                 mipsolver.model_->col_upper_[i] - tolerance) {
        mipsolver.mipdata_->getDomain().changeBound(
            HighsBoundType::kLower, i, mipsolver.model_->col_upper_[i],
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->getDomain().infeasible()) return;
        ++nfixed;
        if (mipsolver.isColInteger(i)) ++nintfixed;
      }
    }
    if (nfixed > 0)
      highsLogDev(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                  "Fixing %d columns (%d integers) sitting at bound at "
                  "analytic center\n",
                  int(nfixed), int(nintfixed));
    mipsolver.mipdata_->getDomain().propagate();
    if (mipsolver.mipdata_->getDomain().infeasible()) return;
  }
}

void HighsMipSolverData::startSymmetryDetection(
    const highs::parallel::TaskGroup& taskGroup,
    std::unique_ptr<SymmetryDetectionData>& symData) {
  symData = std::unique_ptr<SymmetryDetectionData>(new SymmetryDetectionData());
  symData->symDetection.loadModelAsGraph(
      mipsolver.mipdata_->presolvedModel,
      mipsolver.options_mip_->small_matrix_value);
  detectSymmetries = symData->symDetection.initializeDetection();

  if (detectSymmetries) {
    taskGroup.spawn([&]() {
      double startTime = mipsolver.timer_.getWallTime();
      symData->symDetection.run(symData->symmetries, &mipsolver.timer_,
                                mipsolver.options_mip_->time_limit);
      symData->detectionTime = mipsolver.timer_.getWallTime() - startTime;
    });
  } else
    symData.reset();
}

void HighsMipSolverData::finishSymmetryDetection(
    const highs::parallel::TaskGroup& taskGroup,
    std::unique_ptr<SymmetryDetectionData>& symData) {
  taskGroup.sync();

  symmetries = std::move(symData->symmetries);
  std::string symmetry_time =
      mipsolver.options_mip_->timeless_log
          ? ""
          : highsFormatToString(" %.1fs", symData->detectionTime);
  highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
               "\nSymmetry detection completed in%s\n", symmetry_time.c_str());

  if (symmetries.numGenerators == 0) {
    detectSymmetries = false;
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "No symmetry present\n\n");
  } else if (symmetries.orbitopes.size() == 0) {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "Found %d generator(s)\n\n", int(symmetries.numGenerators));

  } else {
    if (symmetries.numPerms != 0) {
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                   "Found %d generator(s) and %d full orbitope(s) acting on %d "
                   "columns\n\n",
                   int(symmetries.numPerms), int(symmetries.orbitopes.size()),
                   int(symmetries.columnToOrbitope.size()));
    } else {
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                   "Found %d full orbitope(s) acting on %d columns\n\n",
                   int(symmetries.orbitopes.size()),
                   int(symmetries.columnToOrbitope.size()));
    }
  }
  symData.reset();

  for (HighsOrbitopeMatrix& orbitope : symmetries.orbitopes)
    orbitope.determineOrbitopeType(cliquetable);

  if (symmetries.numPerms != 0) {
    StabilizerOrbitWorkspace workspace;
    globalOrbits = symmetries.computeStabilizerOrbits(getDomain(), workspace);
  }
}

double HighsMipSolverData::limitsToGap(const double use_lower_bound,
                                       const double use_upper_bound, double& lb,
                                       double& ub) const {
  double offset = mipsolver.model_->offset_;
  lb = use_lower_bound + offset;
  if (std::abs(lb) <= epsilon) lb = 0;
  ub = kHighsInf;
  double gap = kHighsInf;
  if (use_upper_bound != kHighsInf) {
    ub = use_upper_bound + offset;
    if (std::fabs(ub) <= epsilon) ub = 0;
    lb = std::min(ub, lb);
    if (ub == 0.0)
      gap = lb == 0.0 ? 0.0 : kHighsInf;
    else
      gap = (ub - lb) / fabs(ub);
  }
  return gap;
}

double HighsMipSolverData::computeNewUpperLimit(double ub, double mip_abs_gap,
                                                double mip_rel_gap) const {
  double new_upper_limit;
  if (objectiveFunction.isIntegral()) {
    new_upper_limit =
        (std::floor(objectiveFunction.integralScale() * ub - 0.5) /
         objectiveFunction.integralScale());

    if (mip_rel_gap != 0.0)
      new_upper_limit = std::min(
          new_upper_limit,
          ub - std::ceil(mip_rel_gap * fabs(ub + mipsolver.model_->offset_) *
                             objectiveFunction.integralScale() -
                         mipsolver.mipdata_->epsilon) /
                   objectiveFunction.integralScale());

    if (mip_abs_gap != 0.0)
      new_upper_limit = std::min(
          new_upper_limit,
          ub - std::ceil(mip_abs_gap * objectiveFunction.integralScale() -
                         mipsolver.mipdata_->epsilon) /
                   objectiveFunction.integralScale());

    // add feasibility tolerance so that the next best integer feasible solution
    // is definitely included in the remaining search
    new_upper_limit += feastol;
  } else {
    new_upper_limit = std::min(ub - feastol, std::nextafter(ub, -kHighsInf));

    if (mip_rel_gap != 0.0)
      new_upper_limit =
          std::min(new_upper_limit,
                   ub - mip_rel_gap * fabs(ub + mipsolver.model_->offset_));

    if (mip_abs_gap != 0.0)
      new_upper_limit = std::min(new_upper_limit, ub - mip_abs_gap);
  }

  return new_upper_limit;
}

bool HighsMipSolverData::moreHeuristicsAllowed() const {
  // in the beginning of the search and in sub-MIP heuristics we only allow
  // what is proportionally for the currently spent effort plus an initial
  // offset. This is because in a sub-MIP we usually do a truncated search and
  // therefore should not extrapolate the time we spent for heuristics as in
  // the other case. Moreover, since we estimate the total effort for
  // exploring the tree based on the weight of the already pruned nodes, the
  // estimated effort the is not expected to be a good prediction in the
  // beginning.
  if (mipsolver.submip) {
    return heuristic_lp_iterations < total_lp_iterations * heuristic_effort;
  } else if (pruned_treeweight < 1e-3 &&
             num_leaves - num_leaves_before_run < 10 &&
             num_nodes - num_nodes_before_run < 1000) {
    // in the main MIP solver allow an initial offset of 10000 heuristic LP
    // iterations
    if (heuristic_lp_iterations <
        total_lp_iterations * heuristic_effort + 10000)
      return true;
  } else if (heuristic_lp_iterations <
             100000 + ((total_lp_iterations - heuristic_lp_iterations -
                        sb_lp_iterations) >>
                       1)) {
    // compute the node LP iterations in the current run as only those should be
    // used when estimating the total required LP iterations to complete the
    // search
    int64_t heur_iters_curr_run =
        heuristic_lp_iterations - heuristic_lp_iterations_before_run;
    int64_t sb_iters_curr_run = sb_lp_iterations - sb_lp_iterations_before_run;
    int64_t node_iters_curr_run = total_lp_iterations -
                                  total_lp_iterations_before_run -
                                  heur_iters_curr_run - sb_iters_curr_run;
    // now estimate the total fraction of LP iterations that we have spent on
    // heuristics by assuming the node iterations of the current run will
    // grow proportional to the pruned weight of the current tree and the
    // iterations spent for anything else are just added as an offset
    double total_heuristic_effort_estim =
        heuristic_lp_iterations /
        ((total_lp_iterations - node_iters_curr_run) +
         node_iters_curr_run / std::max(0.01, double(pruned_treeweight)));
    // since heuristics help most in the beginning of the search, we want to
    // spent the time we have for heuristics in the first 80% of the tree
    // exploration. Additionally we want to spent the proportional effort
    // of heuristics that is allowed in the first 30% of tree exploration as
    // fast as possible, which is why we have the max(0.3/0.8,...).
    // Hence, in the first 30% of the tree exploration we allow to spent all
    // effort available for heuristics in that part of the search as early as
    // possible, whereas after that we allow the part that is proportionally
    // adequate when we want to spent all available time in the first 80%.
    if (total_heuristic_effort_estim <
        std::max(0.3 / 0.8, std::min(double(pruned_treeweight), 0.8) / 0.8) *
            heuristic_effort) {
      // printf(
      //     "heuristic lp iterations: %ld, total_lp_iterations: %ld, "
      //     "total_heur_effort_estim = %.3f%%\n",
      //     heuristic_lp_iterations, total_lp_iterations,
      //     total_heuristic_effort_estim);
      return true;
    }
  }

  return false;
}

void HighsMipSolverData::removeFixedIndices() {
  integral_cols.erase(
      std::remove_if(integral_cols.begin(), integral_cols.end(),
                     [&](HighsInt col) { return getDomain().isFixed(col); }),
      integral_cols.end());
  integer_cols.erase(
      std::remove_if(integer_cols.begin(), integer_cols.end(),
                     [&](HighsInt col) { return getDomain().isFixed(col); }),
      integer_cols.end());
  implint_cols.erase(
      std::remove_if(implint_cols.begin(), implint_cols.end(),
                     [&](HighsInt col) { return getDomain().isFixed(col); }),
      implint_cols.end());
  continuous_cols.erase(
      std::remove_if(continuous_cols.begin(), continuous_cols.end(),
                     [&](HighsInt col) { return getDomain().isFixed(col); }),
      continuous_cols.end());
}

void HighsMipSolverData::init() {
  postSolveStack.initializeIndexMaps(mipsolver.numRow(), mipsolver.numCol());
  mipsolver.orig_model_ = mipsolver.model_;
  feastol = mipsolver.options_mip_->mip_feasibility_tolerance;
  epsilon = mipsolver.options_mip_->small_matrix_value;
  if (mipsolver.clqtableinit)
    cliquetable.buildFrom(mipsolver.orig_model_, *mipsolver.clqtableinit);
  cliquetable.setMinEntriesForParallelism(
      highs::parallel::num_threads() > 1
          ? mipsolver.options_mip_->mip_min_cliquetable_entries_for_parallelism
          : kHighsIInf);
  if (mipsolver.implicinit) implications.buildFrom(*mipsolver.implicinit);
  heuristic_effort = mipsolver.options_mip_->mip_heuristic_effort;
  detectSymmetries = mipsolver.options_mip_->mip_detect_symmetry;

  firstlpsolobj = -kHighsInf;
  rootlpsolobj = -kHighsInf;
  analyticCenterComputed = false;
  analyticCenterStatus = HighsModelStatus::kNotset;
  maxTreeSizeLog2 = 0;
  numRestarts = 0;
  numRestartsRoot = 0;
  numImprovingSols = 0;
  pruned_treeweight = 0;
  avgrootlpiters = 0;
  num_nodes = 0;
  num_nodes_before_run = 0;
  num_leaves = 0;
  num_leaves_before_run = 0;
  total_repair_lp = 0;
  total_repair_lp_feasible = 0;
  total_repair_lp_iterations = 0;
  total_lp_iterations = 0;
  heuristic_lp_iterations = 0;
  sepa_lp_iterations = 0;
  sb_lp_iterations = 0;
  total_lp_iterations_before_run = 0;
  heuristic_lp_iterations_before_run = 0;
  sepa_lp_iterations_before_run = 0;
  sb_lp_iterations_before_run = 0;
  num_disp_lines = 0;
  numCliqueEntriesAfterPresolve = 0;
  numCliqueEntriesAfterFirstPresolve = 0;
  cliquesExtracted = false;
  rowMatrixSet = false;
  lower_bound = -kHighsInf;
  upper_bound = kHighsInf;
  upper_limit = mipsolver.options_mip_->objective_bound;
  optimality_limit = mipsolver.options_mip_->objective_bound;
  worker_lp_iterations_stop.store(std::numeric_limits<int64_t>::max(),
                                  std::memory_order_relaxed);
  primal_dual_integral.initialise();
  // The P-D integral is measured from the start of the MIP solver
  primal_dual_integral.prev_time = mipsolver.timer_.read();

  if (mipsolver.options_mip_->mip_report_level == 0)
    dispfreq = 0;
  else if (mipsolver.options_mip_->mip_report_level == 1)
    dispfreq = 2000;
  else
    dispfreq = 100;
}

void HighsMipSolverData::runMipPresolve(
    const HighsInt presolve_reduction_limit) {
  mipsolver.timer_.start(mipsolver.timer_.presolve_clock);
  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, presolve_reduction_limit);
  mipsolver.modelstatus_ = presolve.run(postSolveStack);
  presolve_status = presolve.getPresolveStatus();
  mipsolver.timer_.stop(mipsolver.timer_.presolve_clock);

  if (presolve.failed()) return;
  // Report the final presolve reductions unless this is a restart
  if (mipsolver.options_mip_->presolve != kHighsOffString && numRestarts == 0)
    reportPresolveReductions(mipsolver.options_mip_->log_options,
                             presolve_status, *mipsolver.orig_model_,
                             *mipsolver.model_);
}

void HighsMipSolverData::checkAddSolution() {
  // Run when (re-)starting the MIP solver
  if (mipsolver.solution_objective_ != kHighsInf) {
    // Assigning new incumbent
    incumbent = postSolveStack.getReducedPrimalSolution(mipsolver.solution_);
    // return the objective value in the transformed space
    double solobj =
        mipsolver.solution_objective_ * (int)mipsolver.orig_model_->sense_ -
        mipsolver.model_->offset_;
    const bool feasible =
        mipsolver.bound_violation_ <=
            mipsolver.options_mip_->mip_feasibility_tolerance &&
        mipsolver.integrality_violation_ <=
            mipsolver.options_mip_->mip_feasibility_tolerance &&
        mipsolver.row_violation_ <=
            mipsolver.options_mip_->mip_feasibility_tolerance;
    if (numRestarts == 0) {
      // Perform logging and checking of the MIP start solution
      if (feasible) {
        highsLogUser(
            mipsolver.options_mip_->log_options, HighsLogType::kInfo,
            "\nMIP start solution is feasible: objective value is %.12g\n",
            mipsolver.solution_objective_);
      } else {
        // The initial solution shouldn't be infeasible in the
        // original space, since it's tested externally
        highsLogUser(
            mipsolver.options_mip_->log_options, HighsLogType::kWarning,
            "\nMIP start solution is infeasible: objective value is "
            "%.12g with violations:"
            " variable bound = %.4g; constraint bound = %.4g; "
            "integrality = %.4g\n",
            mipsolver.solution_objective_, mipsolver.bound_violation_,
            mipsolver.row_violation_, mipsolver.integrality_violation_);
      }
      // Since the supplied solution (mipsolver.solution_) may not
      // correspond to postsolve applied to incumbent (#3510), compute
      // the objective locally using incumbent, saving the value of
      // mipsolver.solution_objective_ to identify changes
      const double previous_mipsolver_solution_objective =
          mipsolver.solution_objective_;
      double mipsolver_solution_objective =
          double(mipsolver.model_->objectiveCDoubleValue(incumbent));
      mipsolver.solution_objective_ =
          int(mipsolver.orig_model_->sense_) * mipsolver_solution_objective;
      const bool has_start_solution_objective_difference =
          std::fabs(mipsolver.solution_objective_ -
                    previous_mipsolver_solution_objective) /
              (1 + std::fabs(previous_mipsolver_solution_objective)) >
          mipsolver.options_mip_->mip_feasibility_tolerance;
      // This affects solobj - which is use used to determine bounds
      // etc (without the objective offset)
      solobj = mipsolver_solution_objective - mipsolver.model_->offset_;
      if (has_start_solution_objective_difference) {
        const bool reduced = postSolveStack.numReductions() > 0;
        highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                     "MIP start solution is %s%s: objective value is %.12g\n",
                     feasible ? "feasible" : "infeasible",
                     reduced ? " in the reduced space" : "",
                     mipsolver.solution_objective_);
      }
    }
    if (feasible && solobj < upper_bound) {
      double prev_upper_bound = upper_bound;

      upper_bound = solobj;

      bool bound_change = upper_bound != prev_upper_bound;
      if (!mipsolver.submip && bound_change)
        updatePrimalDualIntegral(lower_bound, lower_bound, prev_upper_bound,
                                 upper_bound);

      double new_upper_limit = computeNewUpperLimit(solobj, 0.0, 0.0);

      saveReportMipSolution(new_upper_limit);
      if (new_upper_limit < upper_limit) {
        upper_limit = new_upper_limit;
        optimality_limit =
            computeNewUpperLimit(solobj, mipsolver.options_mip_->mip_abs_gap,
                                 mipsolver.options_mip_->mip_rel_gap);
        nodequeue.setOptimalityLimit(optimality_limit);
      }
    }
    if (!mipsolver.submip && feasible && mipsolver.callback_->user_callback &&
        mipsolver.callback_->active[kCallbackMipSolution]) {
      assert(!mipsolver.submip);
      mipsolver.callback_->clearHighsCallbackOutput();
      mipsolver.callback_->data_out.mip_solution = mipsolver.solution_;
      const bool interrupt = interruptFromCallbackWithData(
          kCallbackMipSolution, mipsolver.solution_objective_,
          "Feasible solution");
      assert(!interrupt);
    }
  }
}

void HighsMipSolverData::runSetup() {
  const HighsLp& model = *mipsolver.model_;

  // Indicate that the first LP has not been solved
  this->getLp().setSolvedFirstLp(false);

  last_disptime = -kHighsInf;
  disptime = 0;

  // Transform the reference of the objective limit and lower/upper
  // bounds from the original model to the current model, undoing the
  // transformation done before restart so that the offset change due
  // to presolve is incorporated. Bound changes are transitory, so no
  // real gap change, and no update to P-D integral is necessary
  upper_limit -= mipsolver.model_->offset_;
  optimality_limit -= mipsolver.model_->offset_;

  lower_bound -= mipsolver.model_->offset_;
  upper_bound -= mipsolver.model_->offset_;

  checkAddSolution();

  if (mipsolver.numCol() == 0)
    addIncumbent(std::vector<double>(), 0, kSolutionSourceEmptyMip);

  redcostfixing = HighsRedcostFixing();
  getPseudoCost() = HighsPseudocost(mipsolver);
  nodequeue.setNumCol(mipsolver.numCol());
  nodequeue.setOptimalityLimit(optimality_limit);

  continuous_cols.clear();
  integer_cols.clear();
  implint_cols.clear();
  integral_cols.clear();

  rowMatrixSet = false;
  if (!rowMatrixSet) {
    rowMatrixSet = true;
    highsSparseTranspose(model.num_row_, model.num_col_, model.a_matrix_.start_,
                         model.a_matrix_.index_, model.a_matrix_.value_,
                         ARstart_, ARindex_, ARvalue_);
    // (re-)initialize number of uplocks and downlocks
    uplocks.assign(model.num_col_, 0);
    downlocks.assign(model.num_col_, 0);
    for (HighsInt i = 0; i != model.num_col_; ++i) {
      HighsInt start = model.a_matrix_.start_[i];
      HighsInt end = model.a_matrix_.start_[i + 1];
      for (HighsInt j = start; j != end; ++j) {
        HighsInt row = model.a_matrix_.index_[j];

        if (model.row_lower_[row] != -kHighsInf) {
          if (model.a_matrix_.value_[j] < 0)
            ++uplocks[i];
          else
            ++downlocks[i];
        }
        if (model.row_upper_[row] != kHighsInf) {
          if (model.a_matrix_.value_[j] < 0)
            ++downlocks[i];
          else
            ++uplocks[i];
        }
      }
    }
  }

  rowintegral.resize(mipsolver.numRow());

#ifdef HIGHS_DEBUGSOL
  if (debugSolution.debugSolActive) {
    debugSolution.debugSolution.clear();
    debugSolution.debugSolution = postSolveStack.getReducedPrimalSolution(
        debugSolution.debugOrigSolution);
    debugSolution.debugSolObjective = 0;
    HighsCDouble debugsolobj = 0.0;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i)
      debugsolobj += static_cast<HighsCDouble>(mipsolver.colCost(i)) *
                     debugSolution.debugSolution[i];
    debugSolution.debugSolObjective = static_cast<double>(debugsolobj);
    debugSolution.registerDomain(getDomain());
    assert(checkSolution(debugSolution.debugSolution));
  }
#endif

  // compute the maximal absolute coefficients to filter propagation
  maxAbsRowCoef.resize(mipsolver.numRow());
  for (HighsInt i = 0; i != mipsolver.numRow(); ++i) {
    double maxabsval = 0.0;

    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];
    bool integral = true;
    for (HighsInt j = start; j != end; ++j) {
      integral = integral && mipsolver.isColIntegral(ARindex_[j]) &&
                 fractionality(ARvalue_[j]) <= epsilon;

      maxabsval = std::max(maxabsval, std::abs(ARvalue_[j]));
    }

    if (integral) {
      if (presolvedModel.row_lower_[i] != -kHighsInf)
        presolvedModel.row_lower_[i] =
            std::ceil(presolvedModel.row_lower_[i] - feastol);

      if (presolvedModel.row_upper_[i] != kHighsInf)
        presolvedModel.row_upper_[i] =
            std::floor(presolvedModel.row_upper_[i] + feastol);
    }

    rowintegral[i] = integral;
    maxAbsRowCoef[i] = maxabsval;
  }

  // compute row activities and propagate all rows once
  objectiveFunction.setupCliquePartition(getDomain(), cliquetable);
  getDomain().setupObjectivePropagation();
  getDomain().computeRowActivities();
  getDomain().propagate();
  if (getDomain().infeasible()) {
    mipsolver.modelstatus_ = HighsModelStatus::kInfeasible;

    updateLowerBound(kHighsInf);

    pruned_treeweight = 1.0;
    return;
  }

  if (model.num_col_ == 0) {
    mipsolver.modelstatus_ = HighsModelStatus::kOptimal;
    return;
  }

  if (checkLimits()) return;
  // extract cliques if they have not been extracted before

  for (HighsInt col : getDomain().getChangedCols())
    implications.cleanupVarbounds(col);
  getDomain().clearChangedCols();

  getLp().getLpSolver().setOptionValue("presolve", kHighsOffString);

  checkObjIntegrality();
  rootlpsol.clear();
  firstlpsol.clear();
  HighsInt num_binary = 0;
  HighsInt num_domain_fixed = 0;
  maxTreeSizeLog2 = 0;
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    switch (mipsolver.variableType(i)) {
      case HighsVarType::kContinuous:
        if (getDomain().isFixed(i)) {
          num_domain_fixed++;
          continue;
        }
        continuous_cols.push_back(i);
        break;
      case HighsVarType::kImplicitInteger:
        if (getDomain().isFixed(i)) {
          num_domain_fixed++;
          continue;
        }
        implint_cols.push_back(i);
        integral_cols.push_back(i);
        break;
      case HighsVarType::kInteger:
        if (getDomain().isFixed(i)) {
          num_domain_fixed++;
          if (fractionality(getDomain().col_lower_[i]) > feastol) {
            // integer variable is fixed to a fractional value -> infeasible
            mipsolver.modelstatus_ = HighsModelStatus::kInfeasible;

            updateLowerBound(kHighsInf);

            pruned_treeweight = 1.0;
            return;
          }
          continue;
        }
        integer_cols.push_back(i);
        integral_cols.push_back(i);
        maxTreeSizeLog2 += (HighsInt)std::ceil(
            std::log2(std::min(1024.0, 1.0 + mipsolver.model_->col_upper_[i] -
                                           mipsolver.model_->col_lower_[i])));
        // NB Since this is for counting the number of times the
        // condition is true using the bitwise operator avoids having
        // any conditional branch whereas using the logical operator
        // would require a branch due to short circuit
        // evaluation. Semantically both is equivalent and correct. If
        // there was any code to be executed for the condition being
        // true then there would be a conditional branch in any case
        // and I would have used the logical to begin with.
        //
        // Hence any compiler warning can be ignored safely
        num_binary +=
            (static_cast<HighsInt>(mipsolver.model_->col_lower_[i] == 0.0) &
             static_cast<HighsInt>(mipsolver.model_->col_upper_[i] == 1.0));
        break;
      case HighsVarType::kSemiContinuous:
      case HighsVarType::kSemiInteger:
        highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kError,
                     "Semicontinuous or semiinteger variables should have been "
                     "reformulated away before HighsMipSolverData::runSetup() "
                     "is called.");
        throw std::logic_error("Unexpected variable type");
    }
  }

  basisTransfer();

  numintegercols = integer_cols.size();
  detectSymmetries = detectSymmetries && num_binary > 0;
  numCliqueEntriesAfterPresolve = cliquetable.getNumEntries();
  HighsInt num_col = mipsolver.numCol();
  HighsInt num_general_integer = numintegercols - num_binary;
  HighsInt num_implied_integer = implint_cols.size();
  HighsInt num_continuous = continuous_cols.size();
  assert(num_col == num_continuous + num_binary + num_general_integer +
                        num_implied_integer + num_domain_fixed);
  if (numRestarts == 0) {
    numCliqueEntriesAfterFirstPresolve = cliquetable.getNumEntries();
    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::kInfo,
        // clang-format off
		 "\nSolving MIP model with:\n"
		 "   %" HIGHSINT_FORMAT " row%s\n"
		 "   %" HIGHSINT_FORMAT " col%s ("
		 "%" HIGHSINT_FORMAT" binary, "
		 "%" HIGHSINT_FORMAT " integer, "
		 "%" HIGHSINT_FORMAT" implied int., "
		 "%" HIGHSINT_FORMAT " continuous, "
		 "%" HIGHSINT_FORMAT " domain fixed)\n"
		 "   %" HIGHSINT_FORMAT " nonzero%s\n"
		 "   Thread count: "
		 "%" HIGHSINT_FORMAT " available cores, "
		 "%" HIGHSINT_FORMAT " logical processors, "
		 "using up to %" HIGHSINT_FORMAT " workers "
		 "across %" HIGHSINT_FORMAT " threads.\n",
        // clang-format on
        mipsolver.numRow(), mipsolver.numRow() == 1 ? "" : "s", num_col,
        num_col == 1 ? "" : "s", num_binary, num_general_integer,
        num_implied_integer, num_continuous, num_domain_fixed,
        mipsolver.numNonzero(), mipsolver.numNonzero() == 1 ? "" : "s",
        HighsInt{static_cast<int>(highs::parallel::available_core_count())},
        HighsInt{static_cast<int>(std::thread::hardware_concurrency())},
        mipsolver.getMaxNumWorkers(), HighsInt{highs::parallel::num_threads()});
  } else {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "Model after restart has "
                 // clang-format off
		 "%" HIGHSINT_FORMAT " row%s, "
		 "%" HIGHSINT_FORMAT " col%s ("
		 "%" HIGHSINT_FORMAT " bin., "
		 "%" HIGHSINT_FORMAT " int., "
		 "%" HIGHSINT_FORMAT " impl., "
		 "%" HIGHSINT_FORMAT " cont., "
		 "%" HIGHSINT_FORMAT " dom.fix.), and "
		 "%" HIGHSINT_FORMAT " nonzero%s\n",
                 // clang-format on
                 mipsolver.numRow(), mipsolver.numRow() == 1 ? "" : "s",
                 num_col, num_col == 1 ? "" : "s", num_binary,
                 num_general_integer, num_implied_integer, num_continuous,
                 num_domain_fixed, mipsolver.numNonzero(),
                 mipsolver.numNonzero() == 1 ? "" : "s");
  }

  heuristics.setupIntCols();

  if (upper_limit == kHighsInf) analyticCenterComputed = false;
  analyticCenterStatus = HighsModelStatus::kNotset;
  analyticCenter.clear();

  symmetries.clear();

  if (numRestarts != 0)
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "\n");
}

double HighsMipSolverData::transformNewIntegerFeasibleSolution(
    const std::vector<double>& sol,
    const bool possibly_store_as_new_incumbent) {
  HighsSolution solution;
  solution.col_value = sol;
  solution.value_valid = true;
  // Perform primal postsolve to get the original column values
  const HighsInt report_col = -kHighsIInf;
  this->reportOriginalPresolvedCol(report_col, sol);
  postSolveStack.undoPrimal(*mipsolver.options_mip_, solution, report_col);
  // Determine the row values, as they aren't computed in primal
  // postsolve
  HighsStatus return_status =
      calculateRowValuesQuad(*mipsolver.orig_model_, solution);
  if (kAllowDeveloperAssert) assert(return_status == HighsStatus::kOk);
  bool allow_try_again = true;
try_again:

  // compute the objective value in the original space
  MipViolation violation;
  HighsCDouble mipsolver_quad_objective_value = 0;
  bool feasible = mipsolver.solutionFeasible(
      mipsolver.orig_model_, solution.col_value, &solution.row_value, violation,
      mipsolver_quad_objective_value);
  double bound_violation_ = 0;
  double row_violation_ = 0;
  double integrality_violation_ = 0;
  violation.copy(bound_violation_, row_violation_, integrality_violation_);
  double mipsolver_objective_value = double(mipsolver_quad_objective_value);
  if (!feasible && allow_try_again) {
    // printf(
    //     "trying to repair sol that is violated by %.12g bounds, %.12g rows,
    //     %.12g "integrality\n", bound_violation_, row_violation_,
    //     integrality_violation_);
    HighsLp fixedModel = *mipsolver.orig_model_;
    fixedModel.integrality_.clear();
    for (HighsInt i = 0; i != mipsolver.orig_model_->num_col_; ++i) {
      if (mipsolver.orig_model_->integrality_[i] == HighsVarType::kInteger) {
        double solval = std::round(solution.col_value[i]);
        fixedModel.col_lower_[i] = std::max(fixedModel.col_lower_[i], solval);
        fixedModel.col_upper_[i] = std::min(fixedModel.col_upper_[i], solval);
      }
    }
    this->total_repair_lp++;
    double time_available = std::max(
        mipsolver.options_mip_->time_limit - mipsolver.timer_.read(), 0.1);
    // Highs instantiation
    Highs tmpSolver;
    tmpSolver.setProfiling(mipsolver.profiling_);
    const bool debug_report = false;
    if (debug_report) {
      tmpSolver.setOptionValue("log_dev_level", 2);
      tmpSolver.setOptionValue("highs_analysis_level", 4);
    } else {
      tmpSolver.setOptionValue("output_flag", false);
    }
    // tmpSolver.setOptionValue("simplex_scale_strategy", 0);
    // tmpSolver.setOptionValue("presolve", kHighsOffString);
    tmpSolver.setOptionValue("time_limit", time_available);
    // Set primal feasibility tolerance for LP solves according to
    // mip_feasibility_tolerance. Interestingly, dual feasibility
    // tolerance not set to smaller tolerance as in
    // HighsLpRelaxationconstructor.
    double mip_primal_feasibility_tolerance =
        mipsolver.options_mip_->mip_feasibility_tolerance;
    tmpSolver.setOptionValue("primal_feasibility_tolerance",
                             mip_primal_feasibility_tolerance);
    // check if only root presolve is allowed
    const bool use_presolve = !mipsolver.options_mip_->mip_root_presolve_only;
    const std::string presolve =
        use_presolve ? kHighsChooseString : kHighsOffString;
    tmpSolver.setOptionValue("presolve", presolve);
    tmpSolver.passModel(std::move(fixedModel));
    // Until a good decision can be made on whether to use simplex,
    // HiPO or IPX to solve an LP without a basis, use simplex
    tmpSolver.setOptionValue("solver", kSimplexString);
    tmpSolver.optimizeLp();
    if (use_presolve &&
        tmpSolver.getModelStatus() == HighsModelStatus::kSolveError) {
      tmpSolver.setOptionValue("presolve", kHighsOffString);
      tmpSolver.optimizeLp();
    }
    this->total_repair_lp_iterations +=
        tmpSolver.getInfo().simplex_iteration_count;
    if (tmpSolver.getInfo().primal_solution_status == kSolutionStatusFeasible) {
      this->total_repair_lp_feasible++;
      solution = tmpSolver.getSolution();
      calculateRowValuesQuad(*mipsolver.orig_model_, solution.col_value,
                             solution.row_value);
      allow_try_again = false;
      goto try_again;
    }
  }

  const double transformed_solobj =
      static_cast<double>(static_cast<HighsInt>(mipsolver.orig_model_->sense_) *
                              mipsolver_quad_objective_value -
                          mipsolver.model_->offset_);

  // Possible MIP solution callback
  if (!mipsolver.submip && feasible && mipsolver.callback_->user_callback &&
      mipsolver.callback_->active[kCallbackMipSolution]) {
    mipsolver.callback_->clearHighsCallbackOutput();
    mipsolver.callback_->data_out.mip_solution = solution.col_value;
    const bool interrupt = interruptFromCallbackWithData(
        kCallbackMipSolution, mipsolver_objective_value, "Feasible solution");
    assert(!interrupt);
  }

  // Catch the case where the repaired solution now has worse objective
  // than the current stored solution
  if (transformed_solobj >= upper_bound && !sol.empty()) {
    return transformed_solobj;
  }

  if (possibly_store_as_new_incumbent) {
    // Store the solution as incumbent in the original space if there
    // is no solution or if it is feasible
    if (feasible) {
      // if (!allow_try_again)
      //   printf("repaired solution with value %g\n",
      //   mipsolver_objective_value);
      // store
      mipsolver.row_violation_ = row_violation_;
      mipsolver.bound_violation_ = bound_violation_;
      mipsolver.integrality_violation_ = integrality_violation_;
      mipsolver.solution_ = std::move(solution.col_value);
      mipsolver.solution_objective_ = mipsolver_objective_value;
    } else {
      bool currentFeasible =
          mipsolver.solution_objective_ != kHighsInf &&
          mipsolver.bound_violation_ <=
              mipsolver.options_mip_->mip_feasibility_tolerance &&
          mipsolver.integrality_violation_ <=
              mipsolver.options_mip_->mip_feasibility_tolerance &&
          mipsolver.row_violation_ <=
              mipsolver.options_mip_->mip_feasibility_tolerance;
      violation.log(mipsolver.options_mip_->log_options,
                    mipsolver_objective_value, "Solution");
      if (!currentFeasible) {
        // if the current incumbent is non existent or also not feasible we
        // still store the new one
        mipsolver.row_violation_ = row_violation_;
        mipsolver.bound_violation_ = bound_violation_;
        mipsolver.integrality_violation_ = integrality_violation_;
        mipsolver.solution_ = std::move(solution.col_value);
        mipsolver.solution_objective_ = mipsolver_objective_value;
      }

      // return infinity so that it is not used for bounding
      return kHighsInf;
    }
  }

  // return the objective value in the transformed space
  return transformed_solobj;
}

double HighsMipSolverData::percentageInactiveIntegers() const {
  return 100.0 *
         (1.0 - static_cast<double>(integer_cols.size() -
                                    cliquetable.getSubstitutions().size()) /
                    numintegercols);
}

void HighsMipSolverData::performRestart() {
  HighsBasis root_basis;
  HighsPseudocostInitialization pscostinit(
      getPseudoCost(), mipsolver.options_mip_->mip_pscost_minreliable,
      postSolveStack);

  mipsolver.pscostinit = &pscostinit;
  ++numRestarts;
  num_leaves_before_run = num_leaves;
  num_nodes_before_run = num_nodes;
  total_lp_iterations_before_run = total_lp_iterations;
  heuristic_lp_iterations_before_run = heuristic_lp_iterations;
  sepa_lp_iterations_before_run = sepa_lp_iterations;
  sb_lp_iterations_before_run = sb_lp_iterations;
  HighsInt numLpRows = getLp().getLp().num_row_;
  HighsInt numModelRows = mipsolver.numRow();
  HighsInt numCuts = numLpRows - numModelRows;
  postSolveStack.appendCutsToModel(numCuts);
  auto integrality = std::move(presolvedModel.integrality_);
  double offset = presolvedModel.offset_;
  presolvedModel = getLp().getLp();
  presolvedModel.offset_ = offset;
  presolvedModel.integrality_ = std::move(integrality);
#ifdef HIGHS_DEBUGSOL
  bool debugSolActive = false;
  std::swap(debugSolution.debugSolActive, debugSolActive);
#endif

  const HighsBasis& basis = firstrootbasis;
  if (basis.valid) {
    // if we have a basis after solving the root LP, we expand it to the
    // original space so that it can be used for constructing a starting basis
    // for the presolved model after the restart
    root_basis.col_status.resize(postSolveStack.getNextColIndex());
    root_basis.row_status.resize(postSolveStack.getNextRowIndex(),
                                 HighsBasisStatus::kBasic);
    root_basis.valid = true;
    root_basis.useful = true;

    for (HighsInt i = 0; i < static_cast<HighsInt>(basis.col_status.size());
         ++i)
      root_basis.col_status[postSolveStack.getOrigColIndex(i)] =
          basis.col_status[i];

    for (HighsInt i = 0; i < static_cast<HighsInt>(basis.row_status.size());
         ++i)
      root_basis.row_status[postSolveStack.getOrigRowIndex(i)] =
          basis.row_status[i];

    mipsolver.rootbasis = &root_basis;
  }

  // Transform the reference of the objective limit and lower/upper
  // bounds to the original model, since offset will generally change
  // in presolve. Bound changes are transitory, so no real gap change,
  // and no update to P-D integral is necessary
  upper_limit += mipsolver.model_->offset_;
  optimality_limit += mipsolver.model_->offset_;

  upper_bound += mipsolver.model_->offset_;
  lower_bound += mipsolver.model_->offset_;

  // remove the current incumbent. Any incumbent is already transformed into the
  // original space and kept there
  incumbent.clear();
  pruned_treeweight = 0;
  nodequeue.clear();
  globalOrbits.reset();

  // Need to be able to set presolve reduction limit separately when
  // restarting - so that bugs in presolve restart can be investigated
  // independently (see #1553)
  //
  // However, when restarting, presolve is (naturally) applied to the
  // presolved problem, so have to control the number of _further_
  // presolve reductions
  //
  // The number of further presolve reductions must be positive,
  // otherwise the MIP solver cycles, hence
  // restart_presolve_reduction_limit cannot be zero
  //
  // Although postSolveStack.numReductions() is size_t, it makes no
  // sense to use presolve_reduction_limit when the number of
  // reductions is vast
  HighsInt num_reductions = HighsInt(postSolveStack.numReductions());
  HighsInt restart_presolve_reduction_limit =
      mipsolver.options_mip_->restart_presolve_reduction_limit;
  assert(restart_presolve_reduction_limit);
  HighsInt further_presolve_reduction_limit =
      restart_presolve_reduction_limit >= 0
          ? num_reductions + restart_presolve_reduction_limit
          : -1;
  runMipPresolve(further_presolve_reduction_limit);

  if (mipsolver.modelstatus_ != HighsModelStatus::kNotset) {
    // transform the objective limit to the current model
    upper_limit -= mipsolver.model_->offset_;
    optimality_limit -= mipsolver.model_->offset_;

    if (mipsolver.modelstatus_ == HighsModelStatus::kOptimal) {
      mipsolver.mipdata_->upper_bound = 0;
      mipsolver.mipdata_->transformNewIntegerFeasibleSolution(
          std::vector<double>());
    } else {
      upper_bound -= mipsolver.model_->offset_;
    }

    // lower_bound still relates to the original model, and the offset
    // is never applied, since MIP solving is complete, and
    // lower_bound is set to upper_bound, so apply the offset now, so
    // that housekeeping in updatePrimalDualIntegral is correct
    lower_bound -= mipsolver.model_->offset_;

    // There must be a gap change, since it's now zero, so always call
    // updatePrimalDualIntegral (unless solving a sub-MIP)
    //
    // Surely there must be a lower bound change
    updateLowerBound(upper_bound, true,
                     mipsolver.modelstatus_ != HighsModelStatus::kOptimal);
    if (mipsolver.solution_objective_ != kHighsInf &&
        mipsolver.modelstatus_ == HighsModelStatus::kInfeasible)
      mipsolver.modelstatus_ = HighsModelStatus::kOptimal;
    return;
  }
  // Bounds are currently in the original space since presolve will have
  // changed offset_
#ifdef HIGHS_DEBUGSOL
  debugSolution.debugSolActive = debugSolActive;
#endif
  runSetup();
  if (mipsolver.terminate()) return;

  postSolveStack.removeCutsFromModel(numCuts);

  // HighsNodeQueue oldNodeQueue;
  // std::swap(nodequeue, oldNodeQueue);

  // Ensure master worker is pointing to the correct cut and conflict pools
  if (!workers.empty()) {
    workers[0].setCutPool(&getCutPool());
    workers[0].setConflictPool(&getConflictPool());
    workers[0].setGlobalDomain(&getDomain());
    workers[0].setPseudocost(&getPseudoCost());
    workers[0].upper_bound = upper_bound;
    workers[0].upper_limit = upper_limit;
    workers[0].optimality_limit = optimality_limit;
  }

  // remove the pointer into the stack-space of this function
  if (mipsolver.rootbasis == &root_basis) mipsolver.rootbasis = nullptr;
  mipsolver.pscostinit = nullptr;
}

void HighsMipSolverData::basisTransfer() {
  // if a root basis is given, construct a basis for the root LP from
  // in the reduced problem space after presolving
  if (mipsolver.rootbasis) {
    const HighsInt numRow = mipsolver.numRow();
    const HighsInt numCol = mipsolver.numCol();
    firstrootbasis.col_status.assign(numCol, HighsBasisStatus::kNonbasic);
    firstrootbasis.row_status.assign(numRow, HighsBasisStatus::kNonbasic);
    firstrootbasis.valid = true;
    firstrootbasis.alien = true;
    firstrootbasis.useful = true;

    for (HighsInt i : postSolveStack.getOrigRows()) {
      HighsBasisStatus status =
          mipsolver.rootbasis->row_status[postSolveStack.getOrigRowIndex(i)];
      firstrootbasis.row_status[i] = status;
    }

    for (HighsInt i : postSolveStack.getOrigCols()) {
      HighsBasisStatus status =
          mipsolver.rootbasis->col_status[postSolveStack.getOrigColIndex(i)];
      firstrootbasis.col_status[i] = status;
    }
  }
}

const std::vector<double>& HighsMipSolverData::getSolution() const {
  return incumbent;
}

bool HighsMipSolverData::addIncumbent(const std::vector<double>& sol,
                                      double solobj, const int solution_source,
                                      const bool print_display_line,
                                      const bool is_user_solution) {
  assert(!parallelLockActive());
  const bool execute_mip_solution_callback =
      !is_user_solution && !mipsolver.submip &&
      (mipsolver.callback_->user_callback
           ? mipsolver.callback_->active[kCallbackMipSolution]
           : false);
  // Determine whether the potential new incumbent should be
  // transformed
  //
  // Happens if solobj improves on the upper bound or the MIP solution
  // callback is active
  const bool possibly_store_as_new_incumbent = solobj < upper_bound;
  const bool get_transformed_solution =
      possibly_store_as_new_incumbent || execute_mip_solution_callback;
  // Get the transformed objective and solution if required
  const double transformed_solobj =
      get_transformed_solution ? transformNewIntegerFeasibleSolution(
                                     sol, possibly_store_as_new_incumbent)
                               : 0;
  const bool highs_solution_report = false;
  if (solution_source == kSolutionSourceHighsSolution && highs_solution_report
      //&& possibly_store_as_new_incumbent
  ) {
    std::stringstream ss;
    ss.str(std::string());
    ss << highsFormatToString(
        "HighsMipSolverData::addIncumbent HiGHS solution Obj "
        "= %15.8g; UB = %15.8g; Obj-UB = %11.4g; PossAdd = %s",
        solobj, upper_bound, solobj - upper_bound,
        possibly_store_as_new_incumbent ? "T" : "F");
    if (possibly_store_as_new_incumbent)
      ss << highsFormatToString(
          "; TransObj = %15.8g; TransObj-UB = %11.4g; TransSolobj < UB %s",
          transformed_solobj, transformed_solobj - upper_bound,
          transformed_solobj < upper_bound ? "T" : "F");
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "%s\n", ss.str().c_str());
    fflush(stdout);
  }
  if (possibly_store_as_new_incumbent) {
    solobj = transformed_solobj;
    if (solobj >= upper_bound) return false;

    double prev_upper_bound = upper_bound;

    upper_bound = solobj;
    for (HighsMipWorker& worker : workers) {
      worker.upper_bound = upper_bound;
    }

    bool bound_change = upper_bound != prev_upper_bound;
    if (!mipsolver.submip && bound_change)
      updatePrimalDualIntegral(lower_bound, lower_bound, prev_upper_bound,
                               upper_bound);

    // Assigning new incumbent
    incumbent = sol;
    double new_upper_limit = computeNewUpperLimit(solobj, 0.0, 0.0);

    if (!is_user_solution && !mipsolver.submip)
      saveReportMipSolution(new_upper_limit);
    if (new_upper_limit < upper_limit) {
      ++numImprovingSols;
      upper_limit = new_upper_limit;
      optimality_limit =
          computeNewUpperLimit(solobj, mipsolver.options_mip_->mip_abs_gap,
                               mipsolver.options_mip_->mip_rel_gap);
      nodequeue.setOptimalityLimit(optimality_limit);
      for (HighsMipWorker& worker : workers) {
        worker.upper_limit = upper_limit;
        worker.optimality_limit = optimality_limit;
      }
      debugSolution.newIncumbentFound();
      getDomain().propagate();
      if (!getDomain().infeasible())
        redcostfixing.propagateRootRedcost(mipsolver);

      // Two calls to printDisplayLine added for completeness,
      // ensuring that when the root node has an integer solution, a
      // logging line is issued

      if (getDomain().infeasible()) {
        pruned_treeweight = 1.0;
        nodequeue.clear();
        if (print_display_line)
          printDisplayLine(solution_source);  // Added for completeness
        return true;
      }
      cliquetable.extractObjCliques(mipsolver);
      if (getDomain().infeasible()) {
        pruned_treeweight = 1.0;
        nodequeue.clear();
        if (print_display_line)
          printDisplayLine(solution_source);  // Added for completeness
        return true;
      }
      pruned_treeweight += nodequeue.performBounding(upper_limit);
      printDisplayLine(solution_source);
    }
  } else if (incumbent.empty())
    // Assigning new incumbent
    incumbent = sol;

  return true;
}

static std::array<char, 22> convertToPrintString(int64_t val) {
  decltype(convertToPrintString(std::declval<int64_t>())) printString = {};
  double l = std::log10(std::max(1.0, double(val)));
  switch (int(l)) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
      std::snprintf(printString.data(), printString.size(), "%" PRId64, val);
      break;
    case 6:
    case 7:
    case 8:
      std::snprintf(printString.data(), printString.size(), "%" PRId64 "k",
                    val / 1000);
      break;
    default:
      std::snprintf(printString.data(), printString.size(), "%" PRId64 "m",
                    val / 1000000);
  }

  return printString;
}

static std::array<char, 22> convertToPrintString(double val,
                                                 const char* trailingStr = "") {
  decltype(convertToPrintString(std::declval<double>(),
                                std::declval<char*>())) printString = {};
  double l = std::abs(val) == kHighsInf
                 ? 0.0
                 : std::log10(std::max(1e-6, std::abs(val)));
  switch (int(l)) {
    case 0:
    case 1:
    case 2:
    case 3:
      std::snprintf(printString.data(), printString.size(), "%.10g%s", val,
                    trailingStr);
      break;
    case 4:
      std::snprintf(printString.data(), printString.size(), "%.11g%s", val,
                    trailingStr);
      break;
    case 5:
      std::snprintf(printString.data(), printString.size(), "%.12g%s", val,
                    trailingStr);
      break;
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
      std::snprintf(printString.data(), printString.size(), "%.13g%s", val,
                    trailingStr);
      break;
    default:
      std::snprintf(printString.data(), printString.size(), "%.9g%s", val,
                    trailingStr);
  }

  return printString;
}

void HighsMipSolverData::printSolutionSourceKey() const {
  std::stringstream ss;
  // Last MipSolutionSource enum is kSolutionSourceCleanup - which is
  // not a solution source, but used to force the last logging line to
  // be printed
  const int last_enum = kSolutionSourceCount - 1;
  // Set the index of the last solution source to be printed in each
  // line of the key. Four or five can be printed, depending on the
  // lengths of the solution source strings in that line
  std::vector<int> limits = {4, 9, 14, last_enum};
  assert(last_enum > limits[limits.size() - 2]);

  ss.str(std::string());
  for (int k = 0; k < limits[0]; k++) {
    if (k == 0) {
      ss << "\nSrc: ";
    } else {
      ss << "; ";
    }
    ss << solutionSourceToString(k) << " => "
       << solutionSourceToString(k, false);
  }
  highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
               "%s;\n", ss.str().c_str());
  int to_line = limits.size() - 1;
  for (int line = 0; line < to_line; line++) {
    ss.str(std::string());
    for (int k = limits[line]; k < limits[line + 1]; k++) {
      if (k == limits[line]) {
        ss << "     ";
      } else {
        ss << "; ";
      }
      ss << solutionSourceToString(k) << " => "
         << solutionSourceToString(k, false);
    }
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "%s%s\n", ss.str().c_str(), line < to_line - 1 ? ";" : "");
  }
}

void HighsMipSolverData::printDisplayLine(const int solution_source) {
  // MIP logging method
  //
  // Note that if the original problem is a maximization, the cost
  // coefficients are negated so that the MIP solver only solves a
  // minimization. Hence, in preparing to print the display line, the
  // dual bound (lb) is always less than the primal bound (ub). When
  // printed, the sense of the optimization is applied so that the
  // values printed correspond to the original objective.

  // No point in computing all the logging values if logging is off
  bool output_flag = *mipsolver.options_mip_->log_options.output_flag;
  if (!output_flag) return;

  bool timeless_log = mipsolver.options_mip_->timeless_log;
  disptime = timeless_log ? disptime + 1 : mipsolver.timer_.read();
  if (solution_source == kSolutionSourceNone &&
      disptime - last_disptime <
          mipsolver.options_mip_->mip_min_logging_interval)
    return;
  last_disptime = disptime;
  std::string time_string =
      timeless_log ? "" : highsFormatToString(" %7.1fs", disptime);

  if (num_disp_lines % 20 == 0) {
    if (num_disp_lines == 0) printSolutionSourceKey();
    std::string work_string0 = timeless_log ? "   Work" : "      Work      ";
    std::string work_string1 = timeless_log ? "LpIters" : "LpIters     Time";
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 // clang-format off
	"\n        Nodes      |    B&B Tree     |            Objective Bounds              |  Dynamic Constraints | %s\n"
	  "Src  Proc. InQueue |  Leaves   Expl. | BestBound       BestSol              Gap |   Cuts   InLp Confl. | %s\n\n",
                 // clang-format on
                 work_string0.c_str(), work_string1.c_str());

    //"   %7s | %10s | %10s | %10s | %10s | %-15s | %-15s | %7s | %7s "
    //"| %8s | %8s\n",
    //"time", "open nodes", "nodes", "leaves", "lpiters", "dual bound",
    //"primal bound", "cutpool", "confl.", "gap", "explored");
  }

  ++num_disp_lines;

  auto print_nodes = convertToPrintString(num_nodes);
  auto queue_nodes = convertToPrintString(nodequeue.numActiveNodes());
  auto print_leaves = convertToPrintString(num_leaves - num_leaves_before_run);

  double explored = 100 * double(pruned_treeweight);

  double lb;
  double ub;
  double gap = limitsToGap(lower_bound, upper_bound, lb, ub);
  gap *= 1e2;
  if (mipsolver.options_mip_->objective_bound < ub)
    ub = mipsolver.options_mip_->objective_bound;

  auto print_lp_iters = convertToPrintString(total_lp_iterations);
  HighsInt dynamic_constraints_in_lp =
      getLp().numRows() > 0 ? getLp().numRows() - getLp().getNumModelRows() : 0;
  if (upper_bound != kHighsInf) {
    std::array<char, 22> gap_string = {};
    if (gap >= 9999.)
      std::strcpy(gap_string.data(), "Large");
    else
      std::snprintf(gap_string.data(), gap_string.size(), "%.2f%%", gap);

    std::array<char, 22> ub_string;
    if (mipsolver.options_mip_->objective_bound < ub) {
      ub_string =
          convertToPrintString((int)mipsolver.orig_model_->sense_ * ub, "*");
    } else
      ub_string = convertToPrintString((int)mipsolver.orig_model_->sense_ * ub);

    auto lb_string =
        convertToPrintString((int)mipsolver.orig_model_->sense_ * lb);

    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::kInfo,
        // clang-format off
                 " %s %7s %7s   %7s %6.2f%%   %-15s %-15s %8s   %6" HIGHSINT_FORMAT " %6" HIGHSINT_FORMAT " %6" HIGHSINT_FORMAT "   %7s%s\n",
        // clang-format on
        solutionSourceToString(solution_source).c_str(), print_nodes.data(),
        queue_nodes.data(), print_leaves.data(), explored, lb_string.data(),
        ub_string.data(), gap_string.data(), getCutPool().getNumCuts(),
        dynamic_constraints_in_lp, getConflictPool().getNumConflicts(),
        print_lp_iters.data(), time_string.c_str());
  } else {
    std::array<char, 22> ub_string;
    if (mipsolver.options_mip_->objective_bound < ub) {
      ub_string =
          convertToPrintString((int)mipsolver.orig_model_->sense_ * ub, "*");
    } else
      ub_string = convertToPrintString((int)mipsolver.orig_model_->sense_ * ub);

    auto lb_string =
        convertToPrintString((int)mipsolver.orig_model_->sense_ * lb);

    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::kInfo,
        // clang-format off
        " %s %7s %7s   %7s %6.2f%%   %-15s %-15s %8.2f   %6" HIGHSINT_FORMAT " %6" HIGHSINT_FORMAT " %6" HIGHSINT_FORMAT "   %7s%s\n",
        // clang-format on
        solutionSourceToString(solution_source).c_str(), print_nodes.data(),
        queue_nodes.data(), print_leaves.data(), explored, lb_string.data(),
        ub_string.data(), gap, getCutPool().getNumCuts(),
        dynamic_constraints_in_lp, getConflictPool().getNumConflicts(),
        print_lp_iters.data(), time_string.c_str());
  }
  // Check that limitsToBounds yields the same values for the
  // dual_bound, primal_bound (modulo optimization sense) and
  // mip_rel_gap
  double dual_bound;
  double primal_bound;
  double mip_rel_gap;
  limitsToBounds(dual_bound, primal_bound, mip_rel_gap);
  mip_rel_gap *= 1e2;
  assert(dual_bound == (int)mipsolver.orig_model_->sense_ * lb);
  assert(primal_bound == (int)mipsolver.orig_model_->sense_ * ub);
  assert(gap == mip_rel_gap);

  // Possibly interrupt from MIP logging callback
  mipsolver.callback_->clearHighsCallbackOutput();
  const bool interrupt = interruptFromCallbackWithData(
      kCallbackMipLogging, mipsolver.solution_objective_, "MIP logging");
  assert(!interrupt);
}

bool HighsMipSolverData::rootSeparationRound(
    HighsMipWorker& worker, HighsSeparation& sepa, HighsInt& ncuts,
    HighsLpRelaxation::Status& status) {
  int64_t tmpLpIters = -getLp().getNumLpIterations();
  ncuts = sepa.separationRound(getDomain(), status);
  tmpLpIters += getLp().getNumLpIterations();
  avgrootlpiters = getLp().getAvgSolveIters();
  total_lp_iterations += tmpLpIters;
  sepa_lp_iterations += tmpLpIters;

  status = evaluateRootLp(worker);
  if (status == HighsLpRelaxation::Status::kInfeasible) return true;

  const std::vector<double>& solvals =
      getLp().getLpSolver().getSolution().col_value;

  if (mipsolver.submip || incumbent.empty()) {
    heuristics.randomizedRounding(worker, solvals);
    if (mipsolver.options_mip_->mip_heuristic_run_shifting)
      heuristics.shifting(worker, solvals);
    heuristics.flushStatistics(mipsolver, worker);
    status = evaluateRootLp(worker);
    if (status == HighsLpRelaxation::Status::kInfeasible) return true;
  }

  return false;
}

HighsLpRelaxation::Status HighsMipSolverData::evaluateRootLp(
    HighsMipWorker& worker) {
  do {
    getDomain().propagate();

    if (globalOrbits && !getDomain().infeasible())
      globalOrbits->orbitalFixing(getDomain());

    if (getDomain().infeasible()) {
      updateLowerBound(std::min(kHighsInf, upper_bound));
      pruned_treeweight = 1.0;
      num_nodes += 1;
      num_leaves += 1;
      return HighsLpRelaxation::Status::kInfeasible;
    }

    bool lpBoundsChanged = false;
    if (!getDomain().getChangedCols().empty()) {
      lpBoundsChanged = true;
      removeFixedIndices();
      getLp().flushDomain(getDomain());
    }

    bool lpWasSolved = false;
    HighsLpRelaxation::Status status;
    if (lpBoundsChanged ||
        getLp().getLpSolver().getModelStatus() == HighsModelStatus::kNotset) {
      int64_t lpIters = -getLp().getNumLpIterations();
      status = getLp().resolveLp(&getDomain());
      lpIters += getLp().getNumLpIterations();
      total_lp_iterations += lpIters;
      avgrootlpiters = getLp().getAvgSolveIters();
      lpWasSolved = true;

      if (status == HighsLpRelaxation::Status::kUnbounded) {
        if (mipsolver.solution_.empty())
          mipsolver.modelstatus_ = HighsModelStatus::kUnboundedOrInfeasible;
        else
          mipsolver.modelstatus_ = HighsModelStatus::kUnbounded;

        pruned_treeweight = 1.0;
        num_nodes += 1;
        num_leaves += 1;
        return status;
      }

      if (status == HighsLpRelaxation::Status::kOptimal &&
          getLp().getFractionalIntegers().empty() &&
          addIncumbent(getLp().getLpSolver().getSolution().col_value,
                       getLp().getObjective(), kSolutionSourceEvaluateNode)) {
        mipsolver.modelstatus_ = HighsModelStatus::kOptimal;
        updateLowerBound(upper_bound);
        pruned_treeweight = 1.0;
        num_nodes += 1;
        num_leaves += 1;
        return HighsLpRelaxation::Status::kInfeasible;
      }

      if (status == HighsLpRelaxation::Status::kOptimal &&
          mipsolver.options_mip_->mip_heuristic_run_zi_round)
        heuristics.ziRound(worker,
                           getLp().getLpSolver().getSolution().col_value);

    } else
      status = getLp().getStatus();

    if (status == HighsLpRelaxation::Status::kInfeasible) {
      updateLowerBound(std::min(kHighsInf, upper_bound));
      pruned_treeweight = 1.0;
      num_nodes += 1;
      num_leaves += 1;
      return status;
    }

    if (getLp().unscaledDualFeasible(getLp().getStatus())) {
      updateLowerBound(std::max(getLp().getObjective(), lower_bound));

      if (lpWasSolved) {
        redcostfixing.addRootRedcost(
            mipsolver, getLp().getLpSolver().getSolution().col_dual,
            getLp().getObjective());
        if (upper_limit != kHighsInf)
          redcostfixing.propagateRootRedcost(mipsolver);
      }
    }

    if (lower_bound > optimality_limit) {
      pruned_treeweight = 1.0;
      num_nodes += 1;
      num_leaves += 1;
      return HighsLpRelaxation::Status::kInfeasible;
    }

    if (getDomain().getChangedCols().empty()) return status;
  } while (true);
}

static void clockOff(HighsProfiling* profiling) {
  if (!profiling->mip_) return;
  if (profiling->isSubMip()) return;
  // Make sure that exactly one of the following clocks is running
  const int clock0_running =
      profiling->running(kMipClockEvaluateRootNode0) ? 1 : 0;
  const int clock1_running =
      profiling->running(kMipClockEvaluateRootNode1) ? 1 : 0;
  const int clock2_running =
      profiling->running(kMipClockEvaluateRootNode2) ? 1 : 0;
  const bool one_running = clock0_running + clock1_running + clock2_running;
  if (!one_running)
    printf("HighsMipSolverData::clockOff Clocks running are (%d; %d; %d)\n",
           clock0_running, clock1_running, clock2_running);
  assert(one_running);
  if (clock0_running) profiling->stop(kMipClockEvaluateRootNode0);
  if (clock1_running) profiling->stop(kMipClockEvaluateRootNode1);
  if (clock2_running) profiling->stop(kMipClockEvaluateRootNode2);
}

void HighsMipSolverData::evaluateRootNode(HighsMipWorker& worker) {
  const bool compute_analytic_centre = true;
  if (!compute_analytic_centre) printf("NOT COMPUTING ANALYTIC CENTRE!\n");
  HighsInt maxSepaRounds = mipsolver.submip ? 5 : kHighsIInf;
  if (numRestarts == 0)
    maxSepaRounds =
        std::min(HighsInt(2 * std::sqrt(maxTreeSizeLog2)), maxSepaRounds);
  std::unique_ptr<SymmetryDetectionData> symData;
  highs::parallel::TaskGroup tg;
  HighsProfiling* profiling = mipsolver.profiling_;
restart:
  profiling->start(kMipClockEvaluateRootNode0);

  if (detectSymmetries) {
    profiling->start(kMipClockStartSymmetryDetection);
    startSymmetryDetection(tg, symData);
    profiling->stop(kMipClockStartSymmetryDetection);
  }
  if (compute_analytic_centre && !analyticCenterComputed) {
    if (profiling->mip_)
      highsLogUser(
          mipsolver.options_mip_->log_options, HighsLogType::kInfo,
          "MIP-Timing: %11.2g - starting analytic centre calculation\n",
          mipsolver.timer_.read());
    profiling->start(kMipClockStartAnalyticCentreComputation);
    startAnalyticCenterComputation(tg);
    profiling->stop(kMipClockStartAnalyticCentreComputation);
  }

  // lp.getLpSolver().setOptionValue(
  //     "dual_simplex_cost_perturbation_multiplier", 10.0);
  getLp().setIterationLimit();
  getLp().loadModel();
  getDomain().clearChangedCols();
  getLp().setObjectiveLimit(upper_limit);

  updateLowerBound(std::max(lower_bound, getDomain().getObjectiveLowerBound()));

  printDisplayLine();

  // Possibly query existence of an external solution
  if (!mipsolver.submip)
    mipsolver.mipdata_->queryExternalSolution(
        mipsolver.solution_objective_,
        kExternalMipSolutionQueryOriginEvaluateRootNode0);

  // check if only root presolve is allowed
  if (firstrootbasis.valid)
    getLp().getLpSolver().setBasis(firstrootbasis,
                                   "HighsMipSolverData::evaluateRootNode");
  else if (mipsolver.options_mip_->mip_root_presolve_only)
    getLp().getLpSolver().setOptionValue("presolve", kHighsOffString);
  else
    getLp().getLpSolver().setOptionValue("presolve", kHighsOnString);
  if (mipsolver.options_mip_->highs_debug_level)
    getLp().getLpSolver().setOptionValue("output_flag",
                                         mipsolver.options_mip_->output_flag);
  //  lp.getLpSolver().setOptionValue("log_dev_level", kHighsLogDevLevelInfo);
  //  lp.getLpSolver().setOptionValue("log_file",
  //  mipsolver.options_mip_->log_file);

  profiling->start(kMipClockEvaluateRootLp);
  HighsLpRelaxation::Status status = evaluateRootLp(worker);
  profiling->stop(kMipClockEvaluateRootLp);
  if (numRestarts == 0) firstrootlpiters = total_lp_iterations;

  getLp().getLpSolver().setOptionValue("output_flag", false);
  getLp().getLpSolver().setOptionValue("presolve", kHighsOffString);
  getLp().getLpSolver().setOptionValue("parallel", kHighsOffString);

  if (status == HighsLpRelaxation::Status::kInfeasible ||
      status == HighsLpRelaxation::Status::kUnbounded)
    return clockOff(profiling);

  firstlpsol = getLp().getSolution().col_value;
  firstlpsolobj = getLp().getObjective();
  rootlpsolobj = firstlpsolobj;

  if (getLp().getLpSolver().getBasis().valid &&
      getLp().numRows() == mipsolver.numRow())
    firstrootbasis = getLp().getLpSolver().getBasis();
  else {
    // the root basis is later expected to be consistent for the model without
    // cuts so set it to the slack basis if the current basis already includes
    // cuts, e.g. due to a restart
    firstrootbasis.col_status.assign(mipsolver.numCol(),
                                     HighsBasisStatus::kNonbasic);
    firstrootbasis.row_status.assign(mipsolver.numRow(),
                                     HighsBasisStatus::kBasic);
    firstrootbasis.valid = true;
    firstrootbasis.useful = true;
  }

  if (getCutPool().getNumCuts() != 0) {
    assert(numRestarts != 0);
    HighsCutSet cutset;
    profiling->start(kMipClockSeparateLpCuts);
    getCutPool().separateLpCutsAfterRestart(cutset);
    profiling->stop(kMipClockSeparateLpCuts);
#ifdef HIGHS_DEBUGSOL
    for (HighsInt i = 0; i < cutset.numCuts(); ++i) {
      debugSolution.checkCut(cutset.ARindex_.data() + cutset.ARstart_[i],
                             cutset.ARvalue_.data() + cutset.ARstart_[i],
                             cutset.ARstart_[i + 1] - cutset.ARstart_[i],
                             cutset.upper_[i]);
    }
#endif
    getLp().addCuts(cutset);
    profiling->start(kMipClockEvaluateRootLp);
    status = evaluateRootLp(worker);
    profiling->stop(kMipClockEvaluateRootLp);
    getLp().removeObsoleteRows();
    if (status == HighsLpRelaxation::Status::kInfeasible)
      return clockOff(profiling);
  }

  getLp().setIterationLimit(std::max(10000, int(10 * avgrootlpiters)));

  // make sure first line after solving root LP is printed
  last_disptime = -kHighsInf;
  disptime = 0;

  if (mipsolver.options_mip_->mip_heuristic_run_zi_round)
    heuristics.ziRound(worker, firstlpsol);
  profiling->start(kMipClockRandomizedRounding);
  heuristics.randomizedRounding(worker, firstlpsol);
  profiling->stop(kMipClockRandomizedRounding);
  if (mipsolver.options_mip_->mip_heuristic_run_shifting)
    heuristics.shifting(worker, firstlpsol);

  heuristics.flushStatistics(mipsolver, worker);

  profiling->start(kMipClockEvaluateRootLp);
  status = evaluateRootLp(worker);
  profiling->stop(kMipClockEvaluateRootLp);
  if (status == HighsLpRelaxation::Status::kInfeasible)
    return clockOff(profiling);

  rootlpsolobj = firstlpsolobj;
  removeFixedIndices();
  if (mipsolver.options_mip_->mip_allow_restart &&
      mipsolver.options_mip_->presolve != kHighsOffString) {
    double fixingRate = percentageInactiveIntegers();
    if (fixingRate >= 10.0) {
      tg.cancel();
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                   "\n%.1f%% inactive integer columns, restarting\n",
                   fixingRate);
      tg.taskWait();
      profiling->start(kMipClockPerformRestart);
      performRestart();
      profiling->stop(kMipClockPerformRestart);
      ++numRestartsRoot;
      if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
        clockOff(profiling);
        goto restart;
      }

      return clockOff(profiling);
    }
  }

  // begin separation
  if (profiling->mip_) {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - starting  separation\n",
                 mipsolver.timer_.read());
    fflush(stdout);
  }
  profiling->start(kMipClockRootSeparation);
  std::vector<double> avgdirection;
  std::vector<double> curdirection;
  avgdirection.resize(mipsolver.numCol());
  curdirection.resize(mipsolver.numCol());

  HighsInt stall = 0;
  double smoothprogress = 0.0;
  HighsInt nseparounds = 0;
  HighsSeparation sepa(worker);
  sepa.setLpRelaxation(&getLp());

  while (getLp().scaledOptimal(status) &&
         !getLp().getFractionalIntegers().empty() && stall < 3) {
    printDisplayLine();

    if (checkLimits()) {
      profiling->stop(kMipClockRootSeparation);
      return clockOff(profiling);
    }

    if (nseparounds == maxSepaRounds) break;

    removeFixedIndices();

    if (!mipsolver.submip &&
        mipsolver.options_mip_->presolve != kHighsOffString) {
      double fixingRate = percentageInactiveIntegers();
      if (fixingRate >= 10.0) {
        stall = -1;
        break;
      }
    }

    ++nseparounds;

    HighsInt ncuts;

    profiling->start(kMipClockRootSeparationRound);
    const bool root_separation_round_result =
        rootSeparationRound(worker, sepa, ncuts, status);
    profiling->stop(kMipClockRootSeparationRound);
    if (root_separation_round_result) {
      profiling->stop(kMipClockRootSeparation);
      return clockOff(profiling);
    }
    if (nseparounds >= 5 && !mipsolver.submip && !analyticCenterComputed &&
        compute_analytic_centre) {
      if (checkLimits()) {
        profiling->stop(kMipClockRootSeparation);
        return clockOff(profiling);
      }
      profiling->start(kMipClockRootSeparationFinishAnalyticCentreComputation);
      finishAnalyticCenterComputation(tg);
      profiling->stop(kMipClockRootSeparationFinishAnalyticCentreComputation);

      profiling->start(kMipClockRootSeparationCentralRounding);
      heuristics.centralRounding(worker);
      profiling->stop(kMipClockRootSeparationCentralRounding);

      heuristics.flushStatistics(mipsolver, worker);

      if (checkLimits()) {
        profiling->stop(kMipClockRootSeparation);
        return clockOff(profiling);
      }
      profiling->start(kMipClockRootSeparationEvaluateRootLp);
      status = evaluateRootLp(worker);
      profiling->stop(kMipClockRootSeparationEvaluateRootLp);
      if (status == HighsLpRelaxation::Status::kInfeasible) {
        profiling->stop(kMipClockRootSeparation);
        return clockOff(profiling);
      }
    }

    HighsCDouble sqrnorm = 0.0;
    const auto& solvals = getLp().getSolution().col_value;

    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      curdirection[i] = firstlpsol[i] - solvals[i];

      // if (mip.integrality_[i] == 2 && lp.getObjective() > firstobj &&
      //    std::abs(curdirection[i]) > 1e-6)
      //  pseudocost.addObservation(i, -curdirection[i],
      //                            lp.getObjective() - firstobj);

      sqrnorm += curdirection[i] * curdirection[i];
    }
#if 1
    double scale = double(1.0 / sqrt(sqrnorm));
    sqrnorm = 0.0;
    HighsCDouble dotproduct = 0.0;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      avgdirection[i] =
          (scale * curdirection[i] - avgdirection[i]) / nseparounds;
      sqrnorm += avgdirection[i] * avgdirection[i];
      dotproduct += avgdirection[i] * curdirection[i];
    }
#endif

    double progress = double(dotproduct / sqrt(sqrnorm));

    if (nseparounds == 1) {
      smoothprogress = progress;
    } else {
      double alpha = 1.0 / 3.0;
      double nextprogress = (1.0 - alpha) * smoothprogress + alpha * progress;

      if (nextprogress < smoothprogress * 1.01 &&
          (getLp().getObjective() - firstlpsolobj) <=
              (rootlpsolobj - firstlpsolobj) * 1.001)
        ++stall;
      else {
        stall = 0;
      }
      smoothprogress = nextprogress;
    }

    rootlpsolobj = getLp().getObjective();
    getLp().setIterationLimit(std::max(10000, int(10 * avgrootlpiters)));
    if (ncuts == 0) break;

    // Possibly query existence of an external solution
    if (!mipsolver.submip)
      mipsolver.mipdata_->queryExternalSolution(
          mipsolver.solution_objective_,
          kExternalMipSolutionQueryOriginEvaluateRootNode1);
  }
  profiling->stop(kMipClockRootSeparation);
  if (profiling->mip_) {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - completed separation\n",
                 mipsolver.timer_.read());
    fflush(stdout);
  }

  getLp().setIterationLimit();
  profiling->start(kMipClockEvaluateRootLp);
  status = evaluateRootLp(worker);
  profiling->stop(kMipClockEvaluateRootLp);
  if (status == HighsLpRelaxation::Status::kInfeasible)
    return clockOff(profiling);

  rootlpsol = getLp().getLpSolver().getSolution().col_value;
  rootlpsolobj = getLp().getObjective();
  getLp().setIterationLimit(std::max(10000, int(10 * avgrootlpiters)));

  if (mipsolver.options_mip_->mip_heuristic_run_zi_round) {
    heuristics.ziRound(worker, firstlpsol);
    heuristics.flushStatistics(mipsolver, worker);
  }
  if (mipsolver.options_mip_->mip_heuristic_run_shifting) {
    heuristics.shifting(worker, rootlpsol);
    heuristics.flushStatistics(mipsolver, worker);
  }

  if (!analyticCenterComputed && compute_analytic_centre) {
    if (checkLimits()) return clockOff(profiling);

    profiling->start(kMipClockFinishAnalyticCentreComputation);
    finishAnalyticCenterComputation(tg);
    profiling->stop(kMipClockFinishAnalyticCentreComputation);

    profiling->start(kMipClockRootCentralRounding);
    heuristics.centralRounding(worker);
    profiling->stop(kMipClockRootCentralRounding);

    heuristics.flushStatistics(mipsolver, worker);

    // if there are new global bound changes we re-evaluate the LP and do one
    // more separation round
    if (checkLimits()) return clockOff(profiling);
    bool separate = !getDomain().getChangedCols().empty();
    profiling->start(kMipClockEvaluateRootLp);
    status = evaluateRootLp(worker);
    profiling->stop(kMipClockEvaluateRootLp);
    if (status == HighsLpRelaxation::Status::kInfeasible)
      return clockOff(profiling);
    if (separate && getLp().scaledOptimal(status)) {
      HighsInt ncuts;
      profiling->start(kMipClockRootSeparationRound0);
      const bool root_separation_round_result =
          rootSeparationRound(worker, sepa, ncuts, status);
      profiling->stop(kMipClockRootSeparationRound0);
      if (root_separation_round_result) return clockOff(profiling);
      ++nseparounds;
      printDisplayLine();
    }
  }

  printDisplayLine();
  // Possibly query existence of an external solution
  if (!mipsolver.submip)
    mipsolver.mipdata_->queryExternalSolution(
        mipsolver.solution_objective_,
        kExternalMipSolutionQueryOriginEvaluateRootNode2);

  // Possible cut extraction callback
  if (!mipsolver.submip && mipsolver.callback_->user_callback &&
      mipsolver.callback_->callbackActive(kCallbackMipGetCutPool))
    mipsolver.callbackGetCutPool();
  if (checkLimits()) return clockOff(profiling);

  profiling->stop(kMipClockEvaluateRootNode0);
  profiling->start(kMipClockEvaluateRootNode1);
  do {
    if (rootlpsol.empty()) break;
    if (upper_limit != kHighsInf && !moreHeuristicsAllowed()) break;

    if (mipsolver.options_mip_->mip_heuristic_run_root_reduced_cost) {
      profiling->start(kMipClockRootHeuristicsReducedCost);
      heuristics.rootReducedCost(worker);
      profiling->stop(kMipClockRootHeuristicsReducedCost);
      heuristics.flushStatistics(mipsolver, worker);
    }

    if (checkLimits()) return clockOff(profiling);

    // if there are new global bound changes we re-evaluate the LP and do one
    // more separation round
    bool separate = !getDomain().getChangedCols().empty();
    profiling->start(kMipClockEvaluateRootLp);
    status = evaluateRootLp(worker);
    profiling->stop(kMipClockEvaluateRootLp);
    if (status == HighsLpRelaxation::Status::kInfeasible)
      return clockOff(profiling);
    if (separate && getLp().scaledOptimal(status)) {
      HighsInt ncuts;
      profiling->start(kMipClockRootSeparationRound1);
      const bool root_separation_round_result =
          rootSeparationRound(worker, sepa, ncuts, status);
      profiling->stop(kMipClockRootSeparationRound1);
      if (root_separation_round_result) return clockOff(profiling);
      ++nseparounds;
      printDisplayLine();
    }

    if (upper_limit != kHighsInf && !moreHeuristicsAllowed()) break;

    if (checkLimits()) return clockOff(profiling);
    if (mipsolver.options_mip_->mip_heuristic_run_rens) {
      profiling->start(kMipClockRootHeuristicsRens);
      heuristics.RENS(worker, rootlpsol);
      profiling->stop(kMipClockRootHeuristicsRens);
      heuristics.flushStatistics(mipsolver, worker);
    }

    if (checkLimits()) return clockOff(profiling);
    // if there are new global bound changes we re-evaluate the LP and do one
    // more separation round
    separate = !getDomain().getChangedCols().empty();
    profiling->start(kMipClockEvaluateRootLp);
    status = evaluateRootLp(worker);
    profiling->stop(kMipClockEvaluateRootLp);
    if (status == HighsLpRelaxation::Status::kInfeasible)
      return clockOff(profiling);
    if (separate && getLp().scaledOptimal(status)) {
      HighsInt ncuts;
      profiling->start(kMipClockRootSeparationRound2);
      const bool root_separation_round_result =
          rootSeparationRound(worker, sepa, ncuts, status);
      profiling->stop(kMipClockRootSeparationRound2);
      if (root_separation_round_result) return clockOff(profiling);
      ++nseparounds;

      printDisplayLine();
      // Possibly query existence of an external solution
      if (!mipsolver.submip)
        mipsolver.mipdata_->queryExternalSolution(
            mipsolver.solution_objective_,
            kExternalMipSolutionQueryOriginEvaluateRootNode3);
    }

    if (upper_limit != kHighsInf || mipsolver.submip) break;

    if (checkLimits()) return clockOff(profiling);
    profiling->start(kMipClockRootFeasibilityPump);
    heuristics.feasibilityPump(worker);
    profiling->stop(kMipClockRootFeasibilityPump);
    heuristics.flushStatistics(mipsolver, worker);

    if (checkLimits()) return clockOff(profiling);
    profiling->start(kMipClockEvaluateRootLp);
    status = evaluateRootLp(worker);
    profiling->stop(kMipClockEvaluateRootLp);
    if (status == HighsLpRelaxation::Status::kInfeasible)
      return clockOff(profiling);
  } while (false);

  profiling->stop(kMipClockEvaluateRootNode1);
  profiling->start(kMipClockEvaluateRootNode2);
  if (lower_bound > upper_limit) {
    mipsolver.modelstatus_ = HighsModelStatus::kOptimal;
    pruned_treeweight = 1.0;
    num_nodes += 1;
    num_leaves += 1;
    return clockOff(profiling);
  }

  // if there are new global bound changes we re-evaluate the LP and do one
  // more separation round
  bool separate = !getDomain().getChangedCols().empty();
  profiling->start(kMipClockEvaluateRootLp);
  status = evaluateRootLp(worker);
  profiling->stop(kMipClockEvaluateRootLp);
  if (status == HighsLpRelaxation::Status::kInfeasible)
    return clockOff(profiling);
  if (separate && getLp().scaledOptimal(status)) {
    HighsInt ncuts;
    profiling->start(kMipClockRootSeparationRound3);
    const bool root_separation_round_result =
        rootSeparationRound(worker, sepa, ncuts, status);
    profiling->stop(kMipClockRootSeparationRound3);
    if (root_separation_round_result) return clockOff(profiling);
    ++nseparounds;
    printDisplayLine();
  }

  // Possibly query existence of an external solution
  if (!mipsolver.submip)
    mipsolver.mipdata_->queryExternalSolution(
        mipsolver.solution_objective_,
        kExternalMipSolutionQueryOriginEvaluateRootNode4);

  removeFixedIndices();
  if (getLp().getLpSolver().getBasis().valid) getLp().removeObsoleteRows();
  rootlpsolobj = getLp().getObjective();

  printDisplayLine();

  if (lower_bound <= upper_limit) {
    if (!mipsolver.submip && mipsolver.options_mip_->mip_allow_restart &&
        mipsolver.options_mip_->presolve != kHighsOffString) {
      if (!analyticCenterComputed && compute_analytic_centre) {
        profiling->start(kMipClockFinishAnalyticCentreComputation);
        finishAnalyticCenterComputation(tg);
        profiling->stop(kMipClockFinishAnalyticCentreComputation);
      }
      double fixingRate = percentageInactiveIntegers();
      if (fixingRate >= 2.5 + 7.5 * mipsolver.submip ||
          (!mipsolver.submip && fixingRate > 0 && numRestarts == 0)) {
        tg.cancel();
        highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                     "\n%.1f%% inactive integer columns, restarting\n",
                     fixingRate);
        if (stall != -1) maxSepaRounds = std::min(maxSepaRounds, nseparounds);
        tg.taskWait();
        profiling->start(kMipClockPerformRestart);
        performRestart();
        profiling->stop(kMipClockPerformRestart);
        if (mipsolver.terminate()) return;
        ++numRestartsRoot;
        if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
          clockOff(profiling);
          goto restart;
        }
        return clockOff(profiling);
      }
    }

    if (detectSymmetries) {
      finishSymmetryDetection(tg, symData);
      profiling->start(kMipClockEvaluateRootLp);
      status = evaluateRootLp(worker);
      profiling->stop(kMipClockEvaluateRootLp);
      if (status == HighsLpRelaxation::Status::kInfeasible)
        return clockOff(profiling);
    }

    // add the root node to the nodequeue to initialize the search
    nodequeue.emplaceNode(
        std::vector<HighsDomainChange>(), std::vector<HighsInt>(), lower_bound,
        getLp().computeBestEstimate(worker.getPseudocost()), 1);
  }
  // End of HighsMipSolverData::evaluateRootNode()
  clockOff(profiling);
}

bool HighsMipSolverData::checkLimits(int64_t nodeOffset) const {
  const HighsOptions& options = *mipsolver.options_mip_;

  // This MIP instance may have been terminated
  if (terminatorActive())
    if (this->terminatorTerminated()) return true;

  // Possible user interrupt
  if (!mipsolver.submip && !parallelLockActive() &&
      mipsolver.callback_->user_callback) {
    mipsolver.callback_->clearHighsCallbackOutput();
    if (interruptFromCallbackWithData(kCallbackMipInterrupt,
                                      mipsolver.solution_objective_,
                                      "MIP check limits")) {
      if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
        highsLogDev(options.log_options, HighsLogType::kInfo,
                    "User interrupt\n");
        mipsolver.modelstatus_ = HighsModelStatus::kInterrupt;
      }
      return true;
    }
  }
  // Possible termination due to objective being at least as good as
  // the target value
  if (!mipsolver.submip && mipsolver.solution_objective_ < kHighsInf &&
      options.objective_target > -kHighsInf) {
    // Note:
    //
    // Whether the sense is ObjSense::kMinimize or
    // ObjSense::kMaximize, the undefined value of
    // mipsolver.solution_objective_ is kHighsInf, and the default
    // target value is -kHighsInf, so had to rule out these cases in
    // the conditional statement above.
    //
    // mipsolver.solution_objective_ is the actual objective of the
    // MIP - including the offset, and independent of objective sense
    //
    // The target is reached if the objective is below (above) the
    // target value when minimizing (maximizing).
    const int int_sense = int(this->mipsolver.orig_model_->sense_);
    const bool reached_objective_target =
        int_sense * mipsolver.solution_objective_ <
        int_sense * options.objective_target;
    if (reached_objective_target) {
      if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
        highsLogDev(options.log_options, HighsLogType::kInfo,
                    "Reached objective target\n");
        mipsolver.modelstatus_ = HighsModelStatus::kObjectiveTarget;
      }
      return true;
    }
  }

  if (options.mip_max_nodes != kHighsIInf &&
      num_nodes + nodeOffset >= options.mip_max_nodes) {
    if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Reached node limit\n");
      mipsolver.modelstatus_ = HighsModelStatus::kSolutionLimit;
    }
    return true;
  }

  if (options.mip_max_leaves != kHighsIInf &&
      num_leaves >= options.mip_max_leaves) {
    if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Reached leaf node limit\n");
      mipsolver.modelstatus_ = HighsModelStatus::kSolutionLimit;
    }
    return true;
  }

  if (options.mip_max_improving_sols != kHighsIInf &&
      numImprovingSols >= options.mip_max_improving_sols) {
    if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Reached improving solution limit\n");
      mipsolver.modelstatus_ = HighsModelStatus::kSolutionLimit;
    }
    return true;
  }

  //  const double time = mipsolver.timer_.read();
  //  printf("checkLimits: time = %g\n", time);
  if (options.time_limit < kHighsInf &&
      mipsolver.timer_.read() >= options.time_limit) {
    if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Reached time limit\n");
      mipsolver.modelstatus_ = HighsModelStatus::kTimeLimit;
    }
    return true;
  }

  return false;
}

void HighsMipSolverData::checkObjIntegrality() {
  objectiveFunction.checkIntegrality(epsilon);
  if (objectiveFunction.isIntegral() && numRestarts == 0) {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "Objective function is integral with scale %g\n",
                 objectiveFunction.integralScale());
  }
}

void HighsMipSolverData::setupDomainPropagation() {
  const HighsLp& model = *mipsolver.model_;
  highsSparseTranspose(model.num_row_, model.num_col_, model.a_matrix_.start_,
                       model.a_matrix_.index_, model.a_matrix_.value_, ARstart_,
                       ARindex_, ARvalue_);

  getPseudoCost() = HighsPseudocost(mipsolver);

  // compute the maximal absolute coefficients to filter propagation
  maxAbsRowCoef.resize(mipsolver.numRow());
  for (HighsInt i = 0; i != mipsolver.numRow(); ++i) {
    double maxabsval = 0.0;

    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];
    for (HighsInt j = start; j != end; ++j)
      maxabsval = std::max(maxabsval, std::abs(ARvalue_[j]));

    maxAbsRowCoef[i] = maxabsval;
  }

  getDomain() = HighsDomain(mipsolver);
  getDomain().computeRowActivities();
}

void HighsMipSolverData::saveReportMipSolution(const double new_upper_limit) {
  const bool non_improving = new_upper_limit >= upper_limit;
  if (mipsolver.submip) return;
  if (non_improving) return;

  if (mipsolver.callback_->user_callback) {
    if (mipsolver.callback_->active[kCallbackMipImprovingSolution]) {
      mipsolver.callback_->clearHighsCallbackOutput();
      mipsolver.callback_->data_out.mip_solution = mipsolver.solution_;
      const bool interrupt = interruptFromCallbackWithData(
          kCallbackMipImprovingSolution, mipsolver.solution_objective_,
          "Improving solution");
      assert(!interrupt);
    }
  }

  if (mipsolver.options_mip_->mip_improving_solution_save) {
    HighsObjectiveSolution record;
    record.objective = mipsolver.solution_objective_;
    record.col_value = mipsolver.solution_;
    mipsolver.saved_objective_and_solution_.push_back(record);
  }
  FILE* file = mipsolver.improving_solution_file_;
  if (file) {
    writeLpObjective(file, mipsolver.options_mip_->log_options,
                     *(mipsolver.orig_model_), mipsolver.solution_);
    writePrimalSolution(
        file, mipsolver.options_mip_->log_options, *(mipsolver.orig_model_),
        mipsolver.solution_,
        mipsolver.options_mip_->mip_improving_solution_report_sparse);
  }
}

void HighsMipSolverData::limitsToBounds(double& dual_bound,
                                        double& primal_bound,
                                        double& mip_rel_gap) const {
  mip_rel_gap = limitsToGap(lower_bound, upper_bound, dual_bound, primal_bound);
  primal_bound =
      std::min(mipsolver.options_mip_->objective_bound, primal_bound);
  // Adjust objective sense in case of maximization problem
  if (this->mipsolver.orig_model_->sense_ == ObjSense::kMaximize) {
    dual_bound = -dual_bound;
    primal_bound = -primal_bound;
  }
}

void HighsMipSolverData::updateLowerBound(double new_lower_bound,
                                          const bool check_bound_change,
                                          const bool check_prev_data) {
  // Update lower bound
  double prev_lower_bound = lower_bound;
  lower_bound = new_lower_bound;
  if (!mipsolver.submip && lower_bound != prev_lower_bound)
    updatePrimalDualIntegral(prev_lower_bound, lower_bound, upper_bound,
                             upper_bound, check_bound_change, check_prev_data);
}

// Interface to callbackAction, with mipsolver_objective_value since
// incumbent value (mipsolver.solution_objective_) is not right for
// callback_type = kCallbackMipSolution

void HighsMipSolverData::setCallbackDataOut(
    const double mipsolver_objective_value) const {
  double dual_bound;
  double primal_bound;
  double mip_rel_gap;
  limitsToBounds(dual_bound, primal_bound, mip_rel_gap);
  mipsolver.callback_->data_out.running_time = mipsolver.timer_.read();
  mipsolver.callback_->data_out.objective_function_value =
      mipsolver_objective_value;
  mipsolver.callback_->data_out.mip_node_count = mipsolver.mipdata_->num_nodes;
  mipsolver.callback_->data_out.mip_total_lp_iterations =
      mipsolver.mipdata_->total_lp_iterations;
  mipsolver.callback_->data_out.mip_primal_bound = primal_bound;
  mipsolver.callback_->data_out.mip_dual_bound = dual_bound;
  mipsolver.callback_->data_out.mip_gap = mip_rel_gap;
}

bool HighsMipSolverData::interruptFromCallbackWithData(
    const int callback_type, const double mipsolver_objective_value,
    const std::string& message) const {
  if (!mipsolver.callback_->callbackActive(callback_type)) return false;
  assert(!mipsolver.submip);
  setCallbackDataOut(mipsolver_objective_value);
  return mipsolver.callback_->callbackAction(callback_type, message);
}

void HighsMipSolverData::queryExternalSolution(
    const double mipsolver_objective_value,
    const ExternalMipSolutionQueryOrigin external_solution_query_origin) {
  assert(!mipsolver.submip);
  HighsCallback* callback = mipsolver.callback_;
  const bool use_callback =
      callback->user_callback && callback->active[kCallbackMipUserSolution];
  if (use_callback) {
    setCallbackDataOut(mipsolver_objective_value);
    callback->data_out.external_solution_query_origin =
        external_solution_query_origin;
    callback->clearHighsCallbackInput();

    const bool interrupt =
        callback->callbackAction(kCallbackMipUserSolution, "MIP User solution");
    assert(!interrupt);
    if (callback->data_in.user_has_solution) {
      // Objective is assumed to be original_offset +
      // (original_c)^T(original_x), but MIP solver bounds are based on the
      // reduced objective (reduced_c)^T(reduced_x)
      //
      // Now, original_sense*[reduced_offset + (reduced_c)^T(reduced_x)] is an
      // objective in the original space, so
      //
      // f0 + c0^Tx0 = s*(f1 + c1^Tx1)
      //
      // where 0 => original; 1 => reduced
      //
      // This allows the reduced objective value to be deduced as
      //
      // c1^Tx1 = s*(f0 + c0^Tx0) - f1
      //
      // (reduced_c)^T(reduced_x) = original_sense*[original_offset +
      // (original_c)^T(original_x) - reduced_offset]
      const auto& user_solution = callback->data_in.user_solution;
      MipViolation violation;
      HighsCDouble user_solution_quad_objective_value = 0;
      const bool feasible = mipsolver.solutionFeasible(
          mipsolver.orig_model_, user_solution, nullptr, violation,
          user_solution_quad_objective_value);
      double bound_violation_ = 0;
      double row_violation_ = 0;
      double integrality_violation_ = 0;
      violation.copy(bound_violation_, row_violation_, integrality_violation_);
      if (!feasible) {
        violation.log(mipsolver.options_mip_->log_options,
                      double(user_solution_quad_objective_value),
                      "User-supplied solution");
        return;
      }
      std::vector<double> reduced_user_solution;
      reduced_user_solution =
          postSolveStack.getReducedPrimalSolution(user_solution);
      const double reduced_user_solution_objective_value =
          double(HighsInt(mipsolver.orig_model_->sense_) *
                     user_solution_quad_objective_value -
                 mipsolver.model_->offset_);
      const bool print_display_line = true;
      const bool is_user_solution = true;
      addIncumbent(reduced_user_solution, reduced_user_solution_objective_value,
                   kSolutionSourceUserSolution, print_display_line,
                   is_user_solution);
    }
  }
}

HighsInt HighsMipSolverData::terminatorConcurrency() const {
  return mipsolver.terminator_.num_instance;
}

HighsInt HighsMipSolverData::terminatorMyInstance() const {
  return mipsolver.terminator_.my_instance;
}

void HighsMipSolverData::terminatorTerminate() {
  assert(terminatorActive());
  mipsolver.terminator_.terminate();
}

bool HighsMipSolverData::terminatorTerminated() const {
  if (this->terminatorActive())
    mipsolver.termination_status_ = mipsolver.terminator_.terminationStatus();
  return mipsolver.termination_status_ != HighsModelStatus::kNotset;
}

void HighsMipSolverData::terminatorReport() const {
  if (this->terminatorActive())
    mipsolver.terminator_.report(mipsolver.options_mip_->log_options);
}

void HighsMipSolverData::reportOriginalPresolvedCol(
    const HighsInt original_col,
    const std::vector<double>& presolved_solution) {
  if (original_col < 0 || original_col >= mipsolver.orig_model_->num_col_)
    return;
  // Find this column in the presolved model
  HighsInt presolved_col = postSolveStack.getPresolvedColumnIndex(original_col);
  printf("Original col = %d (%s) [%g, %g] %s", int(original_col),
         mipsolver.orig_model_->col_names_[original_col].c_str(),
         mipsolver.orig_model_->col_lower_[original_col],
         mipsolver.orig_model_->col_upper_[original_col],
         mipsolver.orig_model_->integrality_[original_col] ==
                 HighsVarType::kContinuous
             ? "Continuous"
             : "Discrete");
  if (presolved_col > 0) {
    printf("; Presolved col = %d (%s) [%g, %g] %s has value %g\n",
           int(presolved_col),
           mipsolver.model_->col_names_[presolved_col].c_str(),
           mipsolver.model_->col_lower_[presolved_col],
           mipsolver.model_->col_upper_[presolved_col],
           mipsolver.model_->integrality_[presolved_col] ==
                   HighsVarType::kContinuous
               ? "Continuous"
               : "Discrete",
           presolved_solution[presolved_col]);
  } else {
    printf("; Not in presolved model\n");
  }
}

static double possInfRelDiff(const double v0, const double v1,
                             const double den) {
  double rel_diff;
  if (std::fabs(v0) == kHighsInf) {
    if (std::fabs(v1) == kHighsInf) {
      rel_diff = 0;
    } else {
      rel_diff = kHighsInf;
    }
  } else {
    if (std::fabs(v1) == kHighsInf) {
      rel_diff = kHighsInf;
    } else {
      rel_diff = std::fabs(v1 - v0) / std::max(1.0, std::fabs(den));
    }
  }
  return rel_diff;
}

void HighsMipSolverData::updatePrimalDualIntegral(const double from_lower_bound,
                                                  const double to_lower_bound,
                                                  const double from_upper_bound,
                                                  const double to_upper_bound,
                                                  const bool check_bound_change,
                                                  const bool check_prev_data) {
  // Parameters to updatePrimalDualIntegral are lower and upper bounds
  // before/after a change
  //
  // updatePrimalDualIntegral should only be called when there is a
  // change in one of the bounds, except when the final update is
  // made, in which case the bounds must NOT have changed. By default,
  // a check for some bound change is made, unless check_bound_change
  // is false, in which case there is a check for unchanged bounds.
  //
  HighsPrimaDualIntegral& pdi = this->primal_dual_integral;
  // HighsPrimaDualIntegral struct contains the following data
  //
  // * value: Current value of the P-D integral
  //
  // * prev_lb: Value of lb that was computed from to_lower_bound in
  //   the previous call. Used as a check that the value of lb
  //   computed from from_lower_bound in this call is equal - to
  //   within bound_change_tolerance. If not true, then a change in lb
  //   has been missed. Only for checking/debugging
  //
  // * prev_ub: Ditto for upper_bound. Only for checking/debugging
  //
  // * prev_gap: Ditto for gap. Only for checking/debugging
  //
  // * prev_time: Used to determine the time spent at the previous gap

  double from_lb;
  double from_ub;
  const double from_gap =
      this->limitsToGap(from_lower_bound, from_upper_bound, from_lb, from_ub);
  double to_lb;
  double to_ub;
  const double to_gap =
      this->limitsToGap(to_lower_bound, to_upper_bound, to_lb, to_ub);

  const double lb_difference = possInfRelDiff(from_lb, to_lb, to_lb);
  const double ub_difference = possInfRelDiff(from_ub, to_ub, to_ub);
  const double bound_change_tolerance = 0;
  const bool bound_change = lb_difference > bound_change_tolerance ||
                            ub_difference > bound_change_tolerance;

  if (check_bound_change) {
    if (!bound_change) {
      if (from_lower_bound == to_lower_bound &&
          from_upper_bound == to_upper_bound) {
        const double lower_bound_difference =
            possInfRelDiff(from_lower_bound, to_lower_bound, to_lower_bound);
        const double upper_bound_difference =
            possInfRelDiff(from_upper_bound, to_upper_bound, to_upper_bound);
        assert(bound_change);
      }
    }
  } else {
    if (bound_change) {
      if (from_lower_bound != to_lower_bound ||
          from_upper_bound != to_upper_bound) {
        const double lower_bound_difference =
            possInfRelDiff(from_lower_bound, to_lower_bound, to_lower_bound);
        const double upper_bound_difference =
            possInfRelDiff(from_upper_bound, to_upper_bound, to_upper_bound);
        assert(!bound_change);
      }
    }
  }
  if (pdi.value <= -kHighsInf) {
    // First call, so zero the P-D integral and no housekeeping to
    // test
    pdi.value = 0;
    pdi.prev_gap = 1.0;
  } else if (check_prev_data) {
    // updatePrimalDualIntegral has been called previously, so can
    // test housekeeping, except after restart
    //
    // These housekeeping tests check that the previous saved
    // lower/upper bounds and gap are very close to the "from"
    // lower/upper bounds and corresponding gap. They are usually
    // identical, but rounding error can occur when passing through
    // reset, when the old/new offsets are added/subtracted from the
    // bounds due to changes in offset during presolve.
    const double lb_inconsistency =
        possInfRelDiff(from_lb, pdi.prev_lb, pdi.prev_lb);
    const bool lb_consistent = lb_inconsistency < 1e-12;
    const double ub_inconsistency =
        possInfRelDiff(from_ub, pdi.prev_ub, pdi.prev_ub);
    const bool ub_consistent = ub_inconsistency < 1e-12;
    const double gap_inconsistency =
        possInfRelDiff(from_gap, pdi.prev_gap, 1.0);
    const bool gap_consistent = gap_inconsistency < 1e-12;
    assert(lb_consistent);
    assert(ub_consistent);
    assert(gap_consistent);
  }
  // Update the P-D integral using the gap capped at 1, so an infinite
  // gap (no incumbent, no lower bound, or gap bigger than 1)
  // contributes 1
  double time = mipsolver.timer_.read();
  double time_diff = time - pdi.prev_time;
  assert(time_diff >= 0);
  pdi.value += time_diff * std::min(pdi.prev_gap, 1.0);
  pdi.prev_time = time;
  pdi.prev_lb = to_lb;
  pdi.prev_ub = to_ub;
  pdi.prev_gap = to_gap;
}

void HighsPrimaDualIntegral::initialise() { this->value = -kHighsInf; }

void HighsTerminator::clear() {
  this->num_instance = 0;
  this->my_instance = kNoThreadInstance;
  this->record = nullptr;
}

void HighsTerminator::initialise(HighsInt num_instance_, HighsInt my_instance_,
                                 HighsModelStatus* record_) {
  this->clear();
  this->num_instance = num_instance_;
  this->my_instance = my_instance_;
  this->record = record_;
}

HighsInt HighsTerminator::concurrency() const { return this->num_instance; }

void HighsTerminator::terminate() {
  assert(this->record);
  assert(this->my_instance < this->num_instance);
  this->record[this->my_instance] = HighsModelStatus::kHighsInterrupt;
}

HighsModelStatus HighsTerminator::terminationStatus() const {
  assert(this->record);
  for (HighsInt instance = 0; instance < this->num_instance; instance++) {
    if (this->record[instance] != HighsModelStatus::kNotset)
      return this->record[instance];
  }
  return HighsModelStatus::kNotset;
}

void HighsTerminator::report(const HighsLogOptions& log_options) const {
  highsLogUser(log_options, HighsLogType::kInfo, "\nTerminator:        ");
  for (HighsInt instance = 0; instance < this->num_instance; instance++)
    highsLogUser(log_options, HighsLogType::kInfo, " %20d",
                 int(this->record[instance]));
  highsLogUser(log_options, HighsLogType::kInfo, "\n");
}

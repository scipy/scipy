/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkPrimal.h
 * @brief Phase 2 primal simplex solver for HiGHS
 */
#ifndef SIMPLEX_HEKKPRIMAL_H_
#define SIMPLEX_HEKKPRIMAL_H_

#include <utility>

#include "simplex/HEkk.h"
#include "util/HSet.h"

using std::pair;

const SimplexAlgorithm algorithm = SimplexAlgorithm::kPrimal;

/**
 * @brief Primal simplex solver for HiGHS
 */

class HEkkPrimal {
 public:
  HEkkPrimal(HEkk& simplex) : ekk_instance_(simplex) { initialiseInstance(); }
  /**
   * @brief Solve a model instance
   */
  HighsStatus solve(const bool force_phase2 = false);

 private:
  /**
   * @brief Initialise a primal simplex instance
   */
  void initialiseInstance();
  /**
   * @brief Initialise a primal simplex solve
   */
  void initialiseSolve();
  void solvePhase1();
  void solvePhase2();
  void cleanup();
  void rebuild();

  void iterate();
  void chuzc();
  void chooseColumn(const bool hyper_sparse = false);
  bool useVariableIn();
  void phase1ChooseRow();
  void chooseRow();

  void considerBoundSwap();
  void assessPivot();

  void update();

  void updateDual();

  void hyperChooseColumn();
  void hyperChooseColumnStart();
  void hyperChooseColumnClear();
  void hyperChooseColumnChangedInfeasibility(const double infeasibility,
                                             const HighsInt iCol);
  void hyperChooseColumnBasicFeasibilityChange();
  void hyperChooseColumnDualChange();

  void phase1ComputeDual();
  void phase1UpdatePrimal();
  void basicFeasibilityChangeBtran();
  void basicFeasibilityChangePrice();
  void basicFeasibilityChangeUpdateDual();

  void phase2UpdatePrimal(const bool initialise = false);

  void considerInfeasibleValueIn();

  void initialiseDevexFramework();
  void updateDevex();
  void computePrimalSteepestEdgeWeights();
  double computePrimalSteepestEdgeWeight(const HighsInt iVar,
                                         HVector& local_col_aq);
  void updatePrimalSteepestEdgeWeights();
  void updateDualSteepestEdgeWeights();
  void updateFtranDSE(HVector& col_steepest_edge);
  void updateBtranPSE(HVector& col_steepest_edge);

  void updateVerify();

  void iterationAnalysisData();
  void iterationAnalysis();
  void localReportIterHeader();
  void localReportIter(const bool header = false);
  void reportRebuild(const HighsInt reason_for_rebuild = -1);
  void getNonbasicFreeColumnSet();
  void removeNonbasicFreeColumn();
  void adjustPerturbedEquationOut();
  void getBasicPrimalInfeasibility();
  bool correctPrimal(const bool initialise = false);
  void shiftBound(const bool lower, const HighsInt iVar, const double value,
                  const double random_value, double& bound, double& shift);
  void savePrimalRay();
  HighsDebugStatus debugPrimalSimplex(const std::string& message,
                                      const bool initialise = false);
  HighsDebugStatus debugPrimalSteepestEdgeWeights(const std::string& message);
  HighsDebugStatus debugPrimalSteepestEdgeWeights(
      const HighsInt alt_debug_level = -1);

  bool isBadBasisChange();

  // References:
  HEkk& ekk_instance_;

  // Pointers:
  HighsSimplexAnalysis* analysis;

  // Class data members
  HighsInt num_col;
  HighsInt num_row;
  HighsInt num_tot;
  HighsInt solve_phase;
  EdgeWeightMode edge_weight_mode;
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  double objective_target;
  HighsInt rebuild_reason;
  // Pivot related
  HighsInt variable_in;
  HighsInt move_in;
  HighsInt row_out;
  HighsInt variable_out;
  HighsInt move_out;
  double theta_dual;
  double theta_primal;
  double value_in;
  double alpha_col;
  double alpha_row;
  double numericalTrouble;

  HighsInt num_flip_since_rebuild;
  // Primal phase 1 tools
  vector<pair<double, HighsInt>> ph1SorterR;
  vector<pair<double, HighsInt>> ph1SorterT;
  // Edge weights
  // Edge weight
  vector<double> edge_weight_;
  HighsInt num_devex_iterations_;
  HighsInt num_bad_devex_weight_;
  vector<HighsInt> devex_index_;
  // Nonbasic free column data.
  HighsInt num_free_col;
  HSet nonbasic_free_col_set;
  // Hyper-sparse CHUZC data
  bool use_hyper_chuzc = false;
  bool initialise_hyper_chuzc;
  bool done_next_chuzc;
  const HighsInt max_num_hyper_chuzc_candidates = 50;
  HighsInt num_hyper_chuzc_candidates;
  vector<HighsInt> hyper_chuzc_candidate;
  vector<double> hyper_chuzc_measure;
  HSet hyper_chuzc_candidate_set;
  double max_hyper_chuzc_non_candidate_measure;
  double max_changed_measure_value;
  HighsInt max_changed_measure_column;
  const bool report_hyper_chuzc = false;
  // Solve buffer
  HVector row_ep;
  HVector row_ap;
  HVector col_aq;
  HVector col_basic_feasibility_change;
  HVector row_basic_feasibility_change;
  HVector col_steepest_edge;
  HighsRandom random_;  // Just for checking PSE weights

  double max_max_local_primal_infeasibility_;
  double max_max_ignored_violation_;
  double max_max_primal_correction_;
  HighsInt last_header_iteration_count_;

  const HighsInt primal_correction_strategy =
      kSimplexPrimalCorrectionStrategyAlways;
  double debug_max_relative_primal_steepest_edge_weight_error = 0;

  const HighsInt check_iter = 9999999;
  const HighsInt check_column = -2133;
};

#endif /* SIMPLEX_HEKKPRIMAL_H_ */

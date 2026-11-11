/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkDual.h
 * @brief Dual simplex solver for HiGHS
 */
#ifndef SIMPLEX_HEKKDUAL_H_
#define SIMPLEX_HEKKDUAL_H_

#include <set>
#include <string>
#include <vector>

#include "simplex/HEkk.h"
#include "simplex/HEkkDualRHS.h"
#include "simplex/HEkkDualRow.h"
#include "simplex/HSimplex.h"
#include "util/HVector.h"
#include "util/HVectorBase.h"

// Limit on the number of column slices for parallel calculations. SIP
// uses num_threads-2 slices; PAMI uses num_threads-1 slices
const HighsInt kHighsSlicedLimit = kSimplexConcurrencyLimit;
// Was 100, but can't see why this should be higher than
// kSimplexConcurrencyLimit; const double kMaxOkGrowth = 1e4;

/**
 * @brief Dual simplex solver for HiGHS
 */
class HEkkDual {
 public:
  HEkkDual(HEkk& simplex)
      : ekk_instance_(simplex), dualRow(simplex), dualRHS(simplex) {
    initialiseInstance();
    dualRow.setup();
    dualRHS.setup();
    if (!(ekk_instance_.info_.simplex_strategy == kSimplexStrategyDualPlain))
      initialiseInstanceParallel(simplex);
  }

  /**
   * @brief Solve a model instance
   */
  HighsStatus solve(const bool force_phase2 = false);

  const SimplexAlgorithm algorithm = SimplexAlgorithm::kDual;

 public:
  /**
   * @brief Initialise a dual simplex instance
   *
   * Copy dimensions and pointers to matrix, factor and solver-related
   * model data, plus tolerances. Sets up local std::vectors
   * (columnDSE, columnBFRT, column, row_ep and row_ap), and buffers
   * for dualRow and dualRHS.
   */
  void initialiseInstance();

  /**
   * @brief Initialise parallel aspects of a dual simplex instance
   *
   * Sets up data structures for SIP or PAMI
   */
  void initialiseInstanceParallel(HEkk& simplex);

  /**
   * @brief Initialise matrix slices and slices of row_ap or dualRow for SIP or
   * PAMI
   *
   * TODO generalise call slice_matrix[i].setup_lgBs so slice can be
   * used with non-logical initial basis
   */
  void initSlice(
      const HighsInt init_sliced_num  //!< Ideal number of slices - true number
                                      //!< is modified in light of limits
  );

  /**
   * @brief Initialise a dual simplex solve
   */
  void initialiseSolve();

  /**
   * @brief Perform Phase 1 dual simplex iterations
   */
  void solvePhase1();

  /**
   * @brief Perform Phase 2 dual simplex iterations
   */
  void solvePhase2();

  /**
   * @brief Reinvert if INVERT not fresh, then recompute dual and primal values
   *
   * Also collects primal infeasibilities and computes the dual objective value
   */

  void rebuild();

  /**
   * @brief Remove perturbation and recompute the dual solution
   *
   * Also collects primal infeasibilities and computes the dual objective value
   */
  void cleanup();

  /**
   * @brief Perform a single serial dual simplex iteration
   *
   * All the methods it calls have as their first line "if (rebuild_reason)
   * return;", where rebuild_reason is, for example, set to 1 when CHUZR
   * finds no candidate. This causes a break from the inner loop of
   * solvePhase% and, hence, a call to rebuild().
   */
  void iterate();

  /**
   * @brief Perform a single SIP dual simplex iteration
   */
  void iterateTasks();

  /**
   * @brief Perform a single PAMI dual simplex iteration - source code in
   * HEkkDualMulti.cpp
   */
  void iterateMulti();  // in HEkkDualMulti.cpp

  /**
   * @brief Pass the data for the serial iteration analysis, report and rebuild
   * report
   */
  void iterationAnalysisData();

  /**
   * @brief Perform the serial iteration analysis
   */
  void iterationAnalysis();

  /**
   * @brief Pass the data for the PAMI iteration analysis for a minor iteration,
   * report and rebuild report
   */
  void iterationAnalysisMinorData();

  /**
   * @brief Perform the PAMI iteration analysis for a minor iteration
   */
  void iterationAnalysisMinor();

  /**
   * @brief Pass the data for the PAMI iteration analysis for a major iteration
   */
  void iterationAnalysisMajorData();

  /**
   * @brief Perform the PAMI iteration analysis for a major iteration
   */
  void iterationAnalysisMajor();

  /**
   * @brief Single line report after rebuild or cleanup
   */
  void reportRebuild(const HighsInt reason_for_rebuild);

  /**
   * @brief Choose the index of a good row to leave the basis (CHUZR)
   */
  void chooseRow();

  /**
   * @brief Determine whether the updated_edge_weight is accurate enough to
   * be accepted, and update the analysis of weight errors
   */
  bool acceptDualSteepestEdgeWeight(const double updated_edge_weight);

  /**
   * @brief Determine whether the updated_edge_weight error should trigger a new
   * Devex framework
   */
  bool newDevexFramework(const double updated_edge_weight);

  /**
   * @brief Compute pivot row (PRICE) and choose the index of a good column to
   * enter the basis (CHUZC)
   */
  void chooseColumn(HVector* row_ep);
  void improveChooseColumnRow(HVector* row_ep);

  /**
   * @brief Choose the index of a good column to enter the basis (CHUZC) by
   * exploiting slices of the pivotal row - for SIP and PAMI
   */
  void chooseColumnSlice(HVector* row_ep);

  /**
   * @brief Compute the pivotal column (FTRAN)
   */
  void updateFtran();

  /**
   * @brief Compute the RHS changes corresponding to the BFRT
   * (FTRAN-BFRT)
   */
  void updateFtranBFRT();

  /**
   * @brief Compute the std::vector required to update DSE weights - being
   * FTRAN applied to the pivotal column (FTRAN-DSE)
   */
  void updateFtranDSE(HVector* DSE_Vector  //!< Pivotal column as RHS for FTRAN
  );
  /**
   * @brief Compare the pivot value computed row-wise and column-wise
   * and determine whether reinversion is advisable
   */
  void updateVerify();

  /**
   * @brief Update the dual values
   */
  void updateDual();

  /**
   * @brief Update the primal values and any edge weights
   */
  void updatePrimal(HVector* DSE_Vector  //!< FTRANned pivotal column
  );

  /**
   * @brief Update the basic and nonbasic variables, iteration count,
   * invertible representation of the basis matrix and row-wise
   * representation of the nonbasic columns, delete the Freelist entry
   * for the entering column, update the primal value for the row
   * where the basis change has occurred, and set the corresponding
   * primal infeasibility value in dualRHS.work_infeasibility, and
   * then determine whether to reinvert according to the synthetic
   * clock
   */
  void updatePivots();

  void shiftCost(const HighsInt iCol, const double amount);
  void shiftBack(const HighsInt iCol);

  /**
   * @brief Initialise a Devex framework: reference set is all basic
   * variables
   */
  void initialiseDevexFramework();

  /**
   * @brief Interpret the dual edge weight strategy as setting of a mode and
   * other actions
   */
  void interpretDualEdgeWeightStrategy(
      const HighsInt simplex_dual_edge_weight_strategy);

  bool reachedExactObjectiveBound();
  double computeExactDualObjectiveValue(HVector& dual_col, HVector& dual_row);

  /**
   * @brief PAMI: Choose the indices of a good set of rows to leave the
   * basis (CHUZR)
   */
  void majorChooseRow();

  /**
   * @brief PAMI: Perform multiple BTRAN
   */
  void majorChooseRowBtran();

  /**
   * @brief PAMI: Choose the index (from the set of indices) of a good
   * row to leave the basis (CHUZR-MI)
   */
  void minorChooseRow();

  /**
   * @brief PAMI: Update the data during minor iterations
   */
  void minorUpdate();

  /**
   * @brief PAMI: Update the dual values during minor iterations
   */
  void minorUpdateDual();

  /**
   * @brief PAMI: Update the primal values during minor iterations
   */
  void minorUpdatePrimal();

  /**
   * @brief PAMI: Perform a basis change during minor iterations
   */
  void minorUpdatePivots();

  /**
   * @brief PAMI: Update the tableau rows during minor iterations
   */
  void minorUpdateRows();

  /**
   * @brief PAMI: Initialise a new Devex framework during minor iterations
   */
  void minorInitialiseDevexFramework();

  /**
   * @brief PAMI: Perform updates after a set of minor iterations
   */
  void majorUpdate();

  /**
   * @brief PAMI: Prepare for the FTRANs after a set of minor iterations
   */
  void majorUpdateFtranPrepare();

  /**
   * @brief PAMI: Perform the parallel part of multiple FTRANs after a
   * set of minor iterations
   */
  void majorUpdateFtranParallel();

  /**
   * @brief PAMI: Perform the final part of multiple FTRANs after a set
   * of minor iterations
   */
  void majorUpdateFtranFinal();

  /**
   * @brief PAMI: Update the primal values after a set of minor
   * iterations
   */
  void majorUpdatePrimal();

  /**
   * @brief PAMI: Update the invertible representation of the basis
   * matrix after a set of minor iterations
   */
  void majorUpdateFactor();

  /**
   * @brief PAMI: Roll back some iterations if numerical trouble
   * detected when updating the invertible representation of the basis
   * matrix after a set of minor iterations
   */
  void majorRollback();

  // private:
  void possiblyUseLiDualSteepestEdge();
  void computeDualInfeasibilitiesWithFixedVariableFlips();
  void correctDualInfeasibilities(HighsInt& free_infeasibility_count);

  bool proofOfPrimalInfeasibility();
  void saveDualRay();
  void assessPhase1Optimality();
  void assessPhase1OptimalityUnperturbed();
  void exitPhase1ResetDuals();
  void reportOnPossibleLpDualInfeasibility();

  bool checkNonUnitWeightError(const std::string& message) const;
  bool dualInfoOk(const HighsLp& lp) const;
  bool bailoutOnDualObjective();
  HighsDebugStatus debugDualSimplex(const std::string& message,
                                    const bool initialise = false);

  bool isBadBasisChange();
  void assessPossiblyDualUnbounded();

  // Devex scalars
  HighsInt num_devex_iterations =
      0;  //!< Number of Devex iterations with the current framework
  bool new_devex_framework = false;  //!< Set a new Devex framework
  bool minor_new_devex_framework =
      false;  //!< Set a new Devex framework in PAMI minor iterations

  // References:
  HEkk& ekk_instance_;

  // Model
  HighsInt solver_num_row;
  HighsInt solver_num_col;
  HighsInt solver_num_tot;
  double inv_solver_num_row;  // 1.0 / solver_num_row

  const HighsSparseMatrix* a_matrix;
  const HSimplexNla* simplex_nla;
  HighsSimplexAnalysis* analysis;

  const int8_t* jMove;
  const double* workRange;
  const double* baseLower;
  const double* baseUpper;
  double* baseValue;
  double* workDual;
  double* workValue;
  double* colLower;
  double* colUpper;
  double* rowLower;
  double* rowUpper;
  int8_t* nonbasicFlag;

  // Retained value
  bool initial_basis_is_logical_;

  // Options
  EdgeWeightMode edge_weight_mode;
  bool allow_dual_steepest_edge_to_devex_switch;

  double Tp;  // Tolerance for primal
  double primal_feasibility_tolerance;

  double Td;  // Tolerance for dual
  double dual_feasibility_tolerance;
  double objective_bound;

  bool force_phase2;
  HighsInt solve_phase;
  HighsInt rebuild_reason;

  HVector row_ep;
  HVector row_ap;
  HVector col_aq;
  HVector col_BFRT;
  HVector col_DSE;

  HVector dev_row_ep;
  HVector dev_col_DSE;

  HEkkDualRow dualRow;

  // Solving related buffers
  HighsInt dualInfeasCount;

  HEkkDualRHS dualRHS;

  // Simplex pivotal information
  HighsInt row_out;
  HighsInt variable_out;
  HighsInt move_out;  // -1 from small to lower, +1 to upper
  HighsInt variable_in;
  double delta_primal;
  double theta_dual;
  double theta_primal;
  double alpha_col;
  double alpha_row;
  double numericalTrouble;
  // (Local) value of computed weight
  double computed_edge_weight;

  bool check_invert_condition = false;

  // Partitioned coefficient matrix
  HighsInt slice_num;
  HighsInt slice_PRICE;
  HighsInt slice_start[kHighsSlicedLimit + 1];
  HighsSparseMatrix slice_a_matrix[kHighsSlicedLimit];
  HighsSparseMatrix slice_ar_matrix[kHighsSlicedLimit];
  HVector slice_row_ap[kHighsSlicedLimit];
  std::vector<HEkkDualRow> slice_dualRow;

  /**
   * @brief Multiple CHUZR data
   */
  struct MChoice {
    HighsInt row_out;
    double baseValue;
    double baseLower;
    double baseUpper;
    double infeasValue;
    double infeasEdWt;
    double infeasLimit;
    HVector row_ep;
    HVector col_aq;
    HVector col_BFRT;
  };

  /**
   * @brief Multiple minor iteration data
   */
  struct MFinish {
    int8_t move_in;
    double shiftOut;
    std::vector<HighsInt> flipList;

    HighsInt row_out;
    HighsInt variable_out;
    HighsInt variable_in;
    double alpha_row;
    double theta_primal;
    double basicBound;
    double basicValue;
    double EdWt;
    HVector_ptr row_ep;
    HVector_ptr col_aq;
    HVector_ptr col_BFRT;
  };

  HighsInt multi_num;
  HighsInt multi_chosen;
  HighsInt multi_iChoice;
  HighsInt multi_nFinish;
  HighsInt multi_iteration;
  HighsInt multi_chooseAgain;
  MChoice multi_choice[kSimplexConcurrencyLimit];
  MFinish multi_finish[kSimplexConcurrencyLimit];
};

#endif /* SIMPLEX_HEKKDUAL_H_ */

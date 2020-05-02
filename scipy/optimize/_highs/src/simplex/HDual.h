/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDual.h
 * @brief Dual simplex solver for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HDUAL_H_
#define SIMPLEX_HDUAL_H_

#include <set>
#include <string>
#include <vector>

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HCrash.h"
#include "simplex/HDualRHS.h"
#include "simplex/HDualRow.h"
#include "simplex/HSimplex.h"
#include "simplex/HVector.h"

class HFactor;

/**
 * Limit on the number of column slices for parallel calculations. SIP uses
 * num_threads-2 slices; PAMI uses num_threads-1 slices
 */
const int HIGHS_SLICED_LIMIT =
    HIGHS_THREAD_LIMIT;  // Was 100, but can't see why this should be higher
                         // than HIGHS_THREAD_LIMIT;

/**
 * Parameters controlling number of Devex iterations.
 *
 * There is a new Devex framework if either
 *
 * 1) The weight inaccuracy ratio exceeds maxAllowedDevexWeightRatio
 *
 * 2) There have been max(minAbsNumberDevexIterations,
 * numRow/minRlvNumberDevexIterations) Devex iterations
 */
const int minAbsNumberDevexIterations = 25;
const double minRlvNumberDevexIterations = 1e-2;
const double maxAllowedDevexWeightRatio = 3.0;

/**
 * Multiplier used in running average calculations
 */
const double runningAverageMu = 0.05;

/**
 * Candidate persistence cut-off in PAMI
 */
const double pami_cutoff = 0.95;

/**
 * @brief Dual simplex solver for HiGHS
 */
class HDual {
 public:
  HDual(HighsModelObject& model_object)
      : workHMO(model_object), dualRow(model_object), dualRHS(model_object) {
    dualRow.setup();
    for (int i = 0; i < HIGHS_SLICED_LIMIT; i++)
      slice_dualRow.push_back(HDualRow(model_object));
    dualRHS.setup();
  }

  /**
   * @brief Solve a model instance
   */
  HighsStatus solve();

  const SimplexAlgorithm algorithm = SimplexAlgorithm::DUAL;

 public:
  /**
   * @brief Set solver options from simplex options
   */
  void options();
  /**
   * @brief Initialise a dual simplex instance
   *
   * Copy dimensions and pointers to matrix, factor and solver-related
   * model data, plus tolerances. Sets up local std::vectors (columnDSE,
   * columnBFRT, column, row_ep and row_ap), scalars for their average
   * density and buffers for dualRow and dualRHS.
   */
  void init();

  /**
   * @brief Initialise parallel aspects of a dual simplex instance
   *
   * Sets up data structures for SIP or PAMI
   */
  void initParallel();

  /**
   * @brief Initialise matrix slices and slices of row_ap or dualRow for SIP or
   * PAMI
   *
   * TODO generalise call slice_matrix[i].setup_lgBs so slice can be
   * used with non-logical initial basis
   */
  void initSlice(
      const int init_sliced_num  //!< Ideal number of slices - true number
                                 //!< is modified in light of limits
  );

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
   * All the methods it calls have as their first line "if (invertHint)
   * return;", where invertHint is, for example, set to 1 when CHUZR
   * finds no candidate. This causes a break from the inner loop of
   * solve_phase% and, hence, a call to rebuild().
   */
  void iterate();

  /**
   * @brief Perform a single SIP dual simplex iteration
   */
  void iterateTasks();

  /**
   * @brief Perform a single PAMI dual simplex iteration - source code in
   * HDualMulti.cpp
   */
  void iterateMulti();  // in HDualMulti.cpp

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
   * @brief Single line report after rebuild
   */
  void reportRebuild(const int rebuild_invert_hint = -1);

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

  /**
   * @brief Initialise a Devex framework: reference set is all basic
   * variables
   */
  void initialiseDevexFramework(const bool parallel = false);

  /**
   * @brief Interpret the dual edge weight strategy as setting of a mode and
   * other actions
   */
  void interpretDualEdgeWeightStrategy(
      const int simplex_dual_edge_weight_strategy);

  /**
   * @brief Interpret the PRICE strategy as setting of a mode and other actions
   */
  /*
  void interpretPriceStrategy(
                              const int simplex_price_strategy
                              );
  */

  bool reachedExactDualObjectiveValueUpperBound();
  double computeExactDualObjectiveValue();

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

  void assessPhase1Optimality();
  void exitPhase1ResetDuals();
  void reportOnPossibleLpDualInfeasibility();

  bool checkNonUnitWeightError(std::string message);
  bool dualInfoOk(const HighsLp& lp);
  bool bailoutReturn();
  bool bailoutOnTimeIterations();
  bool bailoutOnDualObjective();

  bool solve_bailout;  //!< Set true if control is to be returned immediately to
                       //!< calling function

  // Devex scalars
  int num_devex_iterations =
      0;  //!< Number of Devex iterations with the current framework
  bool new_devex_framework = false;  //!< Set a new Devex framework
  bool minor_new_devex_framework =
      false;  //!< Set a new Devex framework in PAMI minor iterations

  // Model
  HighsModelObject& workHMO;
  int solver_num_row;
  int solver_num_col;
  int solver_num_tot;

  const HMatrix* matrix;
  const HFactor* factor;
  HighsSimplexAnalysis* analysis;

  const int* jMove;
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
  int* nonbasicFlag;

  // Options
  DualEdgeWeightMode dual_edge_weight_mode;
  bool initialise_dual_steepest_edge_weights;
  bool allow_dual_steepest_edge_to_devex_switch;

  const double min_dual_steepest_edge_weight = 1e-4;

  double Tp;  // Tolerance for primal
  double primal_feasibility_tolerance;

  double Td;  // Tolerance for dual
  double dual_feasibility_tolerance;
  double dual_objective_value_upper_bound;

  int solvePhase;
  int invertHint;

  HVector row_ep;
  HVector row_ap;
  HVector col_aq;
  HVector col_BFRT;
  HVector col_DSE;

  HDualRow dualRow;

  // Solving related buffers
  int dualInfeasCount;

  HDualRHS dualRHS;

  // Simplex pivotal information
  int rowOut;
  int columnOut;
  int sourceOut;  // -1 from small to lower, +1 to upper
  int columnIn;
  double deltaPrimal;
  double thetaDual;
  double thetaPrimal;
  double alpha;
  double alphaRow;
  double numericalTrouble;
  // (Local) value of computed weight
  double computed_edge_weight;

  // Partitioned coefficient matrix
  int slice_num;
  int slice_PRICE;
  int slice_start[HIGHS_SLICED_LIMIT + 1];
  HMatrix slice_matrix[HIGHS_SLICED_LIMIT];
  HVector slice_row_ap[HIGHS_SLICED_LIMIT];
  std::vector<HDualRow> slice_dualRow;

  /**
   * @brief Multiple CHUZR data
   */
  struct MChoice {
    int rowOut;
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
    int moveIn;
    double shiftOut;
    std::vector<int> flipList;

    int rowOut;
    int columnOut;
    int columnIn;
    double alphaRow;
    double thetaPrimal;
    double basicBound;
    double basicValue;
    double EdWt;
    HVector_ptr row_ep;
    HVector_ptr col_aq;
    HVector_ptr col_BFRT;
  };

  int multi_num;
  int multi_chosen;
  int multi_iChoice;
  int multi_nFinish;
  int multi_iteration;
  int multi_chooseAgain;
  MChoice multi_choice[HIGHS_THREAD_LIMIT];
  MFinish multi_finish[HIGHS_THREAD_LIMIT];

#ifdef HiGHSDEV
  const bool rp_iter_da = false;                  // true;//
  const bool rp_reinvert_syntheticClock = false;  // true;//
  const bool rp_numericalTrouble = false;         // true;//
#endif
  const double original_multi_build_syntheticTick_mu = 1.5;
  const double multi_build_syntheticTick_mu = 1.0;
  // original_multi_build_syntheticTick_mu;//
  const double numerical_trouble_tolerance = 1e-7;
  const double original_multi_numerical_trouble_tolerance = 1e-8;
  const double multi_numerical_trouble_tolerance = 1e-7;
  // original_multi_numerical_trouble_tolerance;

  const int synthetic_tick_reinversion_min_update_count = 50;
  const int original_multi_synthetic_tick_reinversion_min_update_count = 201;
  const int multi_synthetic_tick_reinversion_min_update_count =
      synthetic_tick_reinversion_min_update_count;
  // original_multi_synthetic_tick_reinversion_min_update_count;

  double build_syntheticTick;
  double total_syntheticTick;
};

#endif /* SIMPLEX_HDUAL_H_ */

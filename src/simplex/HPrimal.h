/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HPrimal.h
 * @brief Phase 2 primal simplex solver for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HPRIMAL_H_
#define SIMPLEX_HPRIMAL_H_

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HSimplex.h"
#include "simplex/HVector.h"

/**
 * @brief Phase 2 primal simplex solver for HiGHS
 *
 * Not an efficient primal simplex solver: just a way of tidying up
 * dual infeasibilities when dual optimality (primal feasibility) has
 * been acheived with the dual simplex method
 */

class HPrimal {
 public:
  HPrimal(HighsModelObject& model_object) : workHMO(model_object) {}
  /**
   * @brief Solve a model instance
   */
  HighsStatus solve();

  /**
   * @brief Perform Phase 2 primal simplex iterations
   */
  void solvePhase2();

  const SimplexAlgorithm algorithm = SimplexAlgorithm::PRIMAL;

 private:
  void primalRebuild();
  void primalChooseColumn();
  void primalChooseRow();
  void primalUpdate();

  void iterationAnalysisData();
  void iterationAnalysis();
  void reportRebuild(const int rebuild_invert_hint = -1);
  bool bailout();
  bool solve_bailout;  //!< Set true if control is to be returned immediately to
                       //!< calling function

  // Model pointer
  HighsModelObject& workHMO;

  int solver_num_col;
  int solver_num_row;
  int solver_num_tot;
  HighsSimplexAnalysis* analysis;

  bool no_free_columns;

  int solvePhase;
  // Pivot related
  int invertHint;
  int columnIn;
  int rowOut;
  int columnOut;
  double thetaDual;
  double thetaPrimal;
  double alpha;
  //  double alphaRow;
  double numericalTrouble;
  int num_flip_since_rebuild;

  // Solve buffer
  HVector row_ep;
  HVector row_ap;
  HVector col_aq;
};

#endif /* SIMPLEX_HPRIMAL_H_ */

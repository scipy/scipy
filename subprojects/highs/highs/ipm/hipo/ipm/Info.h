#ifndef HIPO_INFO_H
#define HIPO_INFO_H

#include "Options.h"
#include "Status.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/ipx/info.h"

namespace hipo {

enum TimeProfileItem {
  // matrix
  kMatrixStructureTime_NE,
  kMatrixStructureTime_AS,
  kMatrixValuesTime,
  // analyse
  kAnalyseTime,
  kAnalyseTime_NE,
  kAnalyseTime_AS,
  // factorise
  kFactoriseTime,
  // solve
  kSolveTime,
  kInsertSplitTime,
  kRefinementResTime,
  kRefinementOmegaTime,
  // ipm phases
  kInitialiseTime,
  kPrepareTime,
  kPredictorTime,
  kCorrectorsTime,
  kStepTime,
  //
  kSolve2x2Time,
  kRecoverTime,
  kRefineTime,
  kResidualsTime,
  //
  KTimeProfileItemCount
};

struct Info {
  // Size of problem, as seen by the solver
  Int m_solver, n_solver;

  // Size of original problem
  Int m_original, n_original;

  // Status of solver, see IpmStatus.h
  Status status = kStatusNotSet;
  Int error = kOk;

  // residuals and objectives of final solution
  double p_res_rel, p_res_abs, d_res_rel, d_res_abs, p_obj, d_obj, pd_gap;

  // Number of ipm iterations performed
  Int ipm_iter = 0;

  // True if ipx was invoked, whether to refine solution or for crossover
  bool ipx_used = false;

  // Info from ipx
  ipx::Info ipx_info;

  // Number of correctors used
  Int correctors;

  // Counters
  Int factor_number{};
  Int solve_number{};

  double times[KTimeProfileItemCount];
};

}  // namespace hipo

#endif
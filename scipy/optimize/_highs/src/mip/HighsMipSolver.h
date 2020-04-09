/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef MIP_HIGHS_MIP_SOLVER_H_
#define MIP_HIGHS_MIP_SOLVER_H_

#include "Highs.h"
#include "lp_data/HighsOptions.h"
#include "mip/SolveMip.h"

enum class HighsMipStatus {
  kOptimal,
  kTimeout,
  kReachedSimplexIterationLimit,
  kError,
  kNodeOptimal,
  kNodeInfeasible,
  kNodeUnbounded,
  kNodeNotOptimal,
  kNodeError,
  kRootNodeOptimal,
  kRootNodeNotOptimal,
  kRootNodeError,
  kMaxNodeReached,
  kUnderDevelopment,
  kTreeExhausted
};

const double unscaled_primal_feasibility_tolerance = 1e-4;
const double unscaled_dual_feasibility_tolerance = 1e-4;

class HighsMipSolver : Highs {
 public:
  HighsMipSolver(const HighsOptions& options, const HighsLp& lp)
      : options_mip_(options), mip_(lp) {}

  HighsMipStatus runMipSolver();

 private:
#ifdef HiGHSDEV
  void writeSolutionForIntegerVariables(Node& node);
#endif
  HighsMipStatus solveRootNode();
  HighsMipStatus solveNode(Node& node, bool hotstart = true);
  HighsMipStatus solveTree(Node& root);
  void reportMipSolverProgress(const HighsMipStatus mip_status);
  void reportMipSolverProgressLine(std::string message,
                                   const bool header = false);
  std::string highsMipStatusToString(const HighsMipStatus mip_status);

  Tree tree_;
  const HighsOptions options_mip_;
  const HighsLp mip_;

  int num_nodes_solved = 0;
  int num_nodes_pruned = 0;
  double root_objective_ = HIGHS_CONST_INF;
};

#endif

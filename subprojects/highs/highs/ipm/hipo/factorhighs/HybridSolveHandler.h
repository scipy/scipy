#ifndef FACTORHIGHS_HYBRID_SOLVE_HANDLER_H
#define FACTORHIGHS_HYBRID_SOLVE_HANDLER_H

#include "SolveHandler.h"

namespace hipo {

class HybridSolveHandler : public SolveHandler {
  const std::vector<std::vector<Int>>& swaps_;
  const std::vector<std::vector<char>>& any_swaps_;
  const std::vector<std::vector<double>>& pivot_2x2_;

  mutable std::vector<std::vector<double>> parallel_gemv_workspace_;
  mutable std::vector<double> serial_gemv_workspace_;

  std::vector<Int> first_child_, next_child_;
  mutable std::vector<std::vector<Int>> task_rows_;
  mutable std::vector<std::vector<double>> task_vals_;

  void processForwardTask(Int task, double* x) const;
  void processBackwardTask(Int task, double* x) const;

  void processForwardSn(Int sn, double* x, std::vector<double>& work,
                        Int task = -1, Int end_col_in_task = -1) const;
  void processBackwardSn(Int sn, double* x, std::vector<double>& work) const;
  void processDiagSn(Int sn, double* x) const;

 public:
  void forwardSolve(double* x) const override;
  void backwardSolve(double* x) const override;
  void diagSolve(double* x) const override;
  void inertia(Int& pos, Int& neg, Int& zero, double tol) const override;

  HybridSolveHandler(const Symbolic& S,
                     const std::vector<std::vector<double>>& sn_columns,
                     const std::vector<std::vector<Int>>& swaps,
                     const std::vector<std::vector<char>>& any_swap,
                     const std::vector<std::vector<double>>& pivot_2x2,
                     DataCollector& data, const FHoptions& options);
};

}  // namespace hipo

#endif
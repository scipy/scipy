#ifndef FACTORHIGHS_SOLVE_HANDLER_H
#define FACTORHIGHS_SOLVE_HANDLER_H

#include <vector>

#include "DataCollector.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "FactorHighsOptions.h"

namespace hipo {

// Interface class to handle different formats of dense matrices during the
// solve phase.
// Any implementation of a specific format needs to define:
// - forwardSolve  : to perform the solve with matrix L
// - backwardSolve : to perform the solve with matrix L^T
// - diagSolve     : to perform the solve with matrix D

class SolveHandler {
 protected:
  const Symbolic& S_;
  const std::vector<std::vector<double>>& sn_columns_;
  DataCollector& data_;

  const FHoptions& options_;

 public:
  SolveHandler(const Symbolic& S,
               const std::vector<std::vector<double>>& sn_columns,
               DataCollector& data, const FHoptions& options);

  // avoid copies
  SolveHandler(const SolveHandler&) = delete;
  SolveHandler& operator=(const SolveHandler&) = delete;

  // virtual destructor
  virtual ~SolveHandler() = default;

  // =================================================================
  // Pure virtual functions.
  // These need to be defined by any derived class.
  // =================================================================
  virtual void forwardSolve(double* x) const = 0;
  virtual void backwardSolve(double* x) const = 0;
  virtual void diagSolve(double* x) const = 0;

  virtual void inertia(Int& pos, Int& neg, Int& zero, double tol) const = 0;
};

}  // namespace hipo

#endif
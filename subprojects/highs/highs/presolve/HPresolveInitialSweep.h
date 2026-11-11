/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef PRESOLVE_HIGHS_PRESOLVE_INITIAL_SWEEP_H_
#define PRESOLVE_HIGHS_PRESOLVE_INITIAL_SWEEP_H_

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"

namespace presolve {

class HighsPostsolveStack;

class HPresolveInitialSweep {
 public:
  enum class Result {
    kOk,
    kPrimalInfeasible,
    kDualInfeasible,
  };

  HPresolveInitialSweep(HighsLp& lp, const HighsOptions& options,
                        const std::vector<HighsBool>& allow_rule,
                        const double primal_feastol);

  Result run(HighsPostsolveStack& postsolve_stack);

  HighsInt numDeletedRows() const { return num_deleted_rows_; }
  HighsInt numDeletedCols() const { return num_deleted_cols_; }

 private:
  HighsLp& lp_;
  const HighsOptions& options_;
  const std::vector<HighsBool>& allow_rule_;
  const double primal_feastol_;
  HighsInt num_deleted_rows_;
  HighsInt num_deleted_cols_;

  Result checkColBounds(HighsInt col, bool& isFixed);
  Result emptyCol(HighsPostsolveStack& postsolve_stack, HighsInt col);
  void removeFixedCol(HighsInt col);
  Result emptyRow(HighsPostsolveStack& postsolve_stack, HighsInt row);
  Result singletonRow(HighsPostsolveStack& postsolve_stack, HighsInt row,
                      HighsInt col, double val);
  double getMaxAbsColVal(HighsInt col) const;
};

}  // namespace presolve
#endif

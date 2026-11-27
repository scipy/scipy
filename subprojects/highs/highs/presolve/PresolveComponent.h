/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file PresolveComponent.h
 * @brief The HiGHS class
 */
#ifndef PRESOLVE_PRESOLVE_COMPONENT_H_
#define PRESOLVE_PRESOLVE_COMPONENT_H_

// Not all necessary, but copied from Presolve.h to avoid non-Linux
// failures
#include <list>
#include <map>
#include <stack>
#include <string>
#include <utility>

#include "lp_data/HighsLp.h"
#include "presolve/HighsPostsolveStack.h"
#include "util/HighsComponent.h"
#include "util/HighsTimer.h"

// Class defining the Presolve Component to be used in HiGHS.
// What used to be in Presolve.h but allowing for further testing and dev.

// The structure of component is general, of the presolve component - presolve
// specific.

struct PresolveComponentData : public HighsComponentData {
  HighsLp reduced_lp_;
  presolve::HighsPostsolveStack postSolveStack;
  HighsSolution recovered_solution_;
  HighsBasis recovered_basis_;
  HighsPresolveLog presolve_log_;

  void clear() {
    is_valid = false;

    postSolveStack = presolve::HighsPostsolveStack();

    reduced_lp_.clear();
    recovered_solution_.clear();
    recovered_basis_.clear();
  }

  virtual ~PresolveComponentData() = default;
};

// HighsComponentInfo is a placeholder for details we want to query from outside
// of HiGHS like execution information. Times are recorded at the end of
// Highs::run()
struct PresolveComponentInfo : public HighsComponentInfo {
  HighsInt n_rows_removed = 0;
  HighsInt n_cols_removed = 0;
  HighsInt n_nnz_removed = 0;

  double presolve_time = 0;
  double postsolve_time = 0;

  virtual ~PresolveComponentInfo() = default;
};

class PresolveComponent : public HighsComponent {
 public:
  void clear() override;

  HighsStatus init(const HighsLp& lp, HighsTimer& timer, bool mip = false);

  HighsPresolveStatus run();

  HighsLp& getReducedProblem() { return data_.reduced_lp_; }
  HighsPresolveLog& getPresolveLog() { return data_.presolve_log_; }

  PresolveComponentInfo info_;
  PresolveComponentData data_;
  const HighsOptions* options_;
  HighsTimer* timer;

  HighsPresolveStatus presolve_status_ = HighsPresolveStatus::kNotPresolved;
  HighsPostsolveStatus postsolve_status_ = HighsPostsolveStatus::kNotPresolved;

  virtual ~PresolveComponent() = default;
};
#endif

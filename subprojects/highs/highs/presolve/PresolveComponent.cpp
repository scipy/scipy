/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file PresolveComponent.cpp
 * @brief The HiGHS class
 */

#include "presolve/PresolveComponent.h"

#include "presolve/HPresolve.h"

HighsStatus PresolveComponent::init(const HighsLp& lp, HighsTimer& timer,
                                    bool mip) {
  data_.postSolveStack.initializeIndexMaps(lp.num_row_, lp.num_col_);
  data_.reduced_lp_ = lp;
  this->timer = &timer;
  return HighsStatus::kOk;
}

HighsPresolveStatus PresolveComponent::run() {
  presolve::HPresolve presolve;
  presolve.setInput(data_.reduced_lp_, *options_,
                    options_->presolve_reduction_limit, timer);
  presolve.run(data_.postSolveStack);
  presolve_status_ = presolve.getPresolveStatus();
  if (!presolve.failed()) data_.presolve_log_ = presolve.getPresolveLog();
  return presolve_status_;
}

void PresolveComponent::clear() { data_.clear(); }

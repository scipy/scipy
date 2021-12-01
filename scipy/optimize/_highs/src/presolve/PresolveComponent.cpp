/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file PresolveComponent.cpp
 * @brief The HiGHS class
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "presolve/PresolveComponent.h"

HighsStatus PresolveComponent::init(const HighsLp& lp, HighsTimer& timer) {
  assert(options_.presolve_on);
  data_.presolve_.push_back(presolve::Presolve(timer));
  data_.presolve_[0].load(lp);
  return HighsStatus::OK;
}

HighsStatus PresolveComponent::setOptions(const HighsOptions& options) {
  if (options.presolve == off_string) {
    options_.presolve_on = false;
    return HighsStatus::OK;
  }

  if (options.presolve != on_string) return HighsStatus::Error;

  assert(options_.presolve_on);
  return HighsStatus::OK;
}

void PresolveComponent::negateReducedLpColDuals(bool reduced) {
  if (reduced)
    for (unsigned int col = 0; col < data_.reduced_solution_.col_dual.size();
         col++)
      data_.reduced_solution_.col_dual[col] =
          data_.reduced_solution_.col_dual[col];
  else
    for (unsigned int col = 0; col < data_.recovered_solution_.col_dual.size();
         col++)
      data_.recovered_solution_.col_dual[col] =
          data_.recovered_solution_.col_dual[col];
  return;
}

void PresolveComponent::negateReducedLpCost() {
  for (unsigned int col = 0; col < data_.reduced_lp_.colCost_.size(); col++)
    data_.reduced_lp_.colCost_[col] = -data_.reduced_lp_.colCost_[col];
  return;
}

HighsPresolveStatus PresolveComponent::run() {
  has_run_ = true;
  assert(data_.presolve_.size() > 0);
  // Set options.
  bool options_ok = presolve::checkOptions(options_);
  if (options_ok) {
    if (options_.order.size() > 0) data_.presolve_[0].order = options_.order;

    // max iterations
    if (options_.iteration_strategy == "num_limit")
      data_.presolve_[0].max_iterations = options_.max_iterations;

    // time limit
    if (options_.time_limit < presolve::inf && options_.time_limit > 0)
      data_.presolve_[0].setTimeLimit(options_.time_limit);

    // order and selection of presolvers
    if (options_.order.size() > 0) data_.presolve_[0].order = options_.order;

    // printing
    if (options_.dev) data_.presolve_[0].iPrint = -1;

    data_.presolve_[0].setNumericalTolerances();

    // Run presolve.
    presolve_status_ = data_.presolve_[0].presolve();
  } else {
    presolve_status_ = HighsPresolveStatus::OptionsError;
  }

  // else: Run default.

  if (presolve_status_ == HighsPresolveStatus::Reduced ||
      presolve_status_ == HighsPresolveStatus::ReducedToEmpty) {
    // Move vectors so no copying happens. presolve does not need that lp
    // any more.
    data_.reduced_lp_.numCol_ = data_.presolve_[0].numCol;
    data_.reduced_lp_.numRow_ = data_.presolve_[0].numRow;
    data_.reduced_lp_.Astart_ = std::move(data_.presolve_[0].Astart);
    data_.reduced_lp_.Aindex_ = std::move(data_.presolve_[0].Aindex);
    data_.reduced_lp_.Avalue_ = std::move(data_.presolve_[0].Avalue);
    data_.reduced_lp_.colCost_ = std::move(data_.presolve_[0].colCost);
    data_.reduced_lp_.colLower_ = std::move(data_.presolve_[0].colLower);
    data_.reduced_lp_.colUpper_ = std::move(data_.presolve_[0].colUpper);
    data_.reduced_lp_.rowLower_ = std::move(data_.presolve_[0].rowLower);
    data_.reduced_lp_.rowUpper_ = std::move(data_.presolve_[0].rowUpper);

    data_.reduced_lp_.sense_ = ObjSense::MINIMIZE;
    data_.reduced_lp_.offset_ = 0;
    data_.reduced_lp_.model_name_ =
        std::move(data_.presolve_[0].modelName);  //"Presolved model";
  }

  return presolve_status_;
}

void PresolveComponent::clear() {
  has_run_ = false;
  data_.clear();
}
namespace presolve {

bool checkOptions(const PresolveComponentOptions& options) {
  // todo: check options in a smart way
  if (options.dev) std::cout << "Checking presolve options... ";

  if (!(options.iteration_strategy == "smart" ||
        options.iteration_strategy == "off" ||
        options.iteration_strategy == "num_limit")) {
    if (options.dev)
      std::cout << "error: iteration strategy unknown: "
                << options.iteration_strategy << "." << std::endl;
    return false;
  }

  if (options.iteration_strategy == "num_limit" && options.max_iterations < 0) {
    if (options.dev)
      std::cout << "warning: negative iteration limit: "
                << options.max_iterations
                << ". Presolve will be run with no limit on iterations."
                << std::endl;
    return false;
  }

  return true;
}

}  // namespace presolve

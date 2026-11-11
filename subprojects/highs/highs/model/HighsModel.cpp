/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsModel.cpp
 * @brief
 */
#include "model/HighsModel.h"

#include <cassert>

bool HighsModel::operator==(const HighsModel& model) const {
  bool equal = equalButForNames(model);
  equal = this->lp_.equalNames(model.lp_) && equal;
  return equal;
}

bool HighsModel::equalButForNames(const HighsModel& model) const {
  bool equal = this->lp_.equalButForNames(model.lp_);
  equal = this->hessian_ == model.hessian_ && equal;
  return equal;
}

HighsStatus HighsModel::lpHessianDimensionsCompatible(
    const HighsLogOptions& log_options) const {
  if (!this->hessian_.isEmpty() && this->lp_.num_col_ < this->hessian_.dim()) {
    highsLogUser(log_options, HighsLogType::kError,
                 "HighsModel: LP column dimension = %d is less than Hessian "
                 "dimension = %d\n",
                 int(this->lp_.num_col_), int(this->hessian_.dim()));
    return HighsStatus::kError;
  }
  return HighsStatus::kOk;
}

HighsStatus HighsModel::assess(const HighsOptions& options) {
  HighsStatus lp_status = this->lp_.assess(options);
  HighsStatus hessian_status = this->hessian_.assess(options);
  if (lp_status == HighsStatus::kError ||
      hessian_status == HighsStatus::kError ||
      this->lpHessianDimensionsCompatible(options.log_options) ==
          HighsStatus::kError)
    return HighsStatus::kError;
  if (lp_status == HighsStatus::kWarning ||
      hessian_status == HighsStatus::kWarning)
    return HighsStatus::kWarning;
  return HighsStatus::kOk;
}

void HighsModel::clear() {
  this->lp_.clear();
  this->hessian_.clear();
}

double HighsModel::objectiveValue(const std::vector<double>& solution) const {
  return this->hessian_.objectiveValue(solution) +
         this->lp_.objectiveValue(solution);
}

void HighsModel::objectiveGradient(const std::vector<double>& solution,
                                   std::vector<double>& gradient) const {
  if (this->isQp()) {
    this->hessian_.product(solution, gradient);
  } else {
    gradient.assign(this->lp_.num_col_, 0);
  }
  for (HighsInt iCol = 0; iCol < this->lp_.num_col_; iCol++)
    gradient[iCol] += this->lp_.col_cost_[iCol];
}

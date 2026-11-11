/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file model/HighsModel.h
 * @brief
 */
#ifndef MODEL_HIGHS_MODEL_H_
#define MODEL_HIGHS_MODEL_H_

#include <vector>

#include "lp_data/HighsLp.h"
#include "model/HighsHessian.h"

class HighsModel;

class HighsModel {
 public:
  HighsLp lp_;
  HighsHessian hessian_;
  bool operator==(const HighsModel& model) const;
  bool equalButForNames(const HighsModel& model) const;
  bool isQp() const {
    return (this->hessian_.dim_ != 0 || this->hessian_.isOracle());
  }
  bool isMip() const { return this->lp_.isMip(); }
  bool isEmpty() const {
    return (this->lp_.num_col_ == 0 && this->lp_.num_row_ == 0);
  }
  bool needsMods(const double infinite_cost) const {
    return this->lp_.needsMods(infinite_cost);
  }
  bool hasMods() const { return this->lp_.hasMods(); }

  HighsStatus lpHessianDimensionsCompatible(
      const HighsLogOptions& log_options) const;

  HighsStatus assess(const HighsOptions& options);
  void clear();
  double objectiveValue(const std::vector<double>& solution) const;
  void objectiveGradient(const std::vector<double>& solution,
                         std::vector<double>& gradient) const;
};

#endif

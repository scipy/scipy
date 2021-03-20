/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HStruct.h
 * @brief Structs for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HSTRUCT_H_
#define LP_DATA_HSTRUCT_H_

#include <vector>

#include "lp_data/HConst.h"

struct HighsIterationCounts {
  int simplex = 0;
  int ipm = 0;
  int crossover = 0;
};

struct HighsScale {
  bool is_scaled_ = false;
  double cost_;
  std::vector<double> col_;
  std::vector<double> row_;
};

struct HighsSolution {
  std::vector<double> col_value;
  std::vector<double> col_dual;
  std::vector<double> row_value;
  std::vector<double> row_dual;
};

struct HighsBasis {
  bool valid_ = false;
  std::vector<HighsBasisStatus> col_status;
  std::vector<HighsBasisStatus> row_status;
};

struct HighsSolutionParams {
  // Input to solution analysis method
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  int primal_status = PrimalDualStatus::STATUS_NOTSET;
  int dual_status = PrimalDualStatus::STATUS_NOTSET;
  // Output from solution analysis method
  double objective_function_value;
  int num_primal_infeasibilities;
  double sum_primal_infeasibilities;
  double max_primal_infeasibility;
  int num_dual_infeasibilities;
  double sum_dual_infeasibilities;
  double max_dual_infeasibility;
};

#endif /* LP_DATA_HSTRUCT_H_ */

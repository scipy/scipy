/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_STATISTICS_HPP__
#define __SRC_LIB_STATISTICS_HPP__

#include <chrono>
#include <vector>

struct Statistics {
  HighsInt phase1_iterations = 0;
  HighsInt num_iterations = 0;
  std::chrono::high_resolution_clock::time_point time_start;
  std::chrono::high_resolution_clock::time_point time_end;

  std::vector<HighsInt> iteration;
  std::vector<HighsInt> nullspacedimension;
  std::vector<double> objval;
  std::vector<double> time;
  std::vector<double> sum_primal_infeasibilities;
  std::vector<HighsInt> num_primal_infeasibilities;
  std::vector<double> density_nullspace;
  std::vector<double> density_factor;
};

#endif

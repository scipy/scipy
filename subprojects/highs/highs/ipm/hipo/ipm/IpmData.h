#ifndef HIPO_IPM_DATA_H
#define HIPO_IPM_DATA_H

#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// structs to collect data during the ipm iterations

struct IpmIterData {
  double sigma_aff = 0.0;
  double sigma = 0.0;
  Int correctors = 0;
  Int num_solves = 0;
  double min_theta = 0.0;
  double max_theta = 0.0;
  double min_prod = 0.0;
  double max_prod = 0.0;
  Int num_small_prod = 0;
  Int num_large_prod = 0;
  double omega = 0.0;
  Int factorisation_used = 0;
};

struct IpmData {
  std::vector<IpmIterData> record;
  void append();
  IpmIterData& back();
};

}  // namespace hipo

#endif
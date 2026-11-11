#ifndef FACTOR_HIGHS_OPTIONS_H
#define FACTOR_HIGHS_OPTIONS_H

#include "FactorHighsSettings.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

struct FHoptions {
  double reg_p = 0.0;
  double reg_d = 0.0;
  Int nb = kBlockSize;
  bool pivoting = true;
  bool one_indexing = false;
  bool parallel_tree = false;
  bool parallel_node = false;
  bool parallel_forward = false;
  bool parallel_backward = false;
  bool parallel_diag = false;
};

}  // namespace hipo

#endif
#ifndef FACTOR_HIGHS_OPTIONS_H
#define FACTOR_HIGHS_OPTIONS_H

#include "FactorHighsSettings.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

struct FHoptions {
  // C-style callback that returns nonzero to stop a long computation early.
  // It must not throw and may be called concurrently.
  using StopCallback = int (*)(void* data);
  StopCallback should_stop = nullptr;
  void* stop_data = nullptr;
  bool stopRequested() const { return should_stop && should_stop(stop_data); }
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
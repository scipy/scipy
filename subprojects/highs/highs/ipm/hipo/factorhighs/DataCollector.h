#ifndef FACTORHIGHS_DATA_COLLECTOR_H
#define FACTORHIGHS_DATA_COLLECTOR_H

#include <atomic>
#include <limits>
#include <mutex>
#include <vector>

#include "FactorHighsSettings.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

struct IterData {
// data of one factorisation
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  double minD = kHighsInf;
  double maxD = 0.0;
  double minL = kHighsInf;
  double maxL = 0.0;
  double max_reg = 0.0;
  Int n_reg_piv = 0;
  Int n_swap = 0;
  Int n_2x2 = 0;
  Int n_wrong_sign = 0;
  double max_wrong_sign = 0.0;
  double M_norm1;
  double M_maxdiag;
#endif
};

// DataCollector is used to collect debug data during the factorisations.
//
// Expensive data related to the factorisation is only collected if
// HIPO_COLLECT_EXPENSIVE_DATA is defined.
//
// HIPO_COLLECT_EXPENSIVE_DATA is selected at compile time, because some data is
// collected many times, during many steps of the factorisation, using locks; if
// HIPO_COLLECT_EXPENSIVE_DATA is off at compile time, then there is no
// performance hit, since the functions are empty and are completely ignored by
// the compiler.
//
// Times are only collected if HIPO_TIMING_LEVEL is defined at compile time, for
// the same reason.
//
// All functions that collect data have a mutex, so that concurrent data
// collection is possible. However, data collection should be used only for
// debugging.
//

class DataCollector {
  // Record of times and BLAS calls
  std::vector<double> times_{};
  std::vector<Int> blas_calls_{};

  // record of data of ipm iterations
  std::vector<IterData> iter_data_record_{};

  // Mutex for concurrent access
  std::mutex mutex_;

 public:
  DataCollector();

  IterData& back();
  void append();

  void sumTime(TimeItems i, double t);
  void countRegPiv();
  void countSwap();
  void count2x2();
  void setWrongSign(double p);
  void setMaxReg(double new_reg);
  void setExtremeEntries(double minD, double maxD, double minoffD,
                         double maxoffD);
  void setNorms(double norm1, double maxdiag);

  void printTimes(const Logger& logger) const;
  void printIter(const Logger& logger) const;
};

}  // namespace hipo

#endif
#ifndef HIPO_OPTION_H
#define HIPO_OPTION_H

#include "Parameters.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsOptions.h"

namespace hipo {

enum class ParallelType { kOff, kChoose, kOn };

struct Options {
  // Solver options
  std::string nla = kHighsChooseString;
  std::string crossover = kHighsOffString;
  std::string ordering = kHighsChooseString;
  std::string factor = kHighsChooseString;

  // Ipm parameters
  Int max_iter = kMaxIterDefault;
  double feasibility_tol = kIpmTolDefault;
  double optimality_tol = kIpmTolDefault;
  double crossover_tol = kIpmTolDefault;
  bool refine_with_ipx = true;
  double time_limit = -1.0;
  Int block_size = 0;
  Int random_seed = 0;
  ParallelType parallel[static_cast<Int>(ParallelTechnique::kCount)];

  // Logging
  bool display = true;
  bool timeless_log = false;
  const HighsLogOptions* log_options = nullptr;

  inline bool getParallel(ParallelTechnique bit) const {
    return static_cast<bool>(parallel[static_cast<Int>(bit)]);
  }
  inline void setParallel(ParallelTechnique bit, ParallelType type) {
    parallel[static_cast<Int>(bit)] = type;
  }
  inline void chooseParallel(ParallelTechnique bit, bool default_behaviour) {
    ParallelType type_default =
        default_behaviour ? ParallelType::kOn : ParallelType::kOff;
    if (parallel[static_cast<Int>(bit)] == ParallelType::kChoose)
      setParallel(bit, type_default);
  }
};

inline bool testParallelBit(Int option, Int bit) { return option & (1 << bit); }

}  // namespace hipo

#endif
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_RUNTIME_HPP__
#define __SRC_LIB_RUNTIME_HPP__

#include "instance.hpp"
#include "qpsolver/qpconst.hpp"
#include "settings.hpp"
#include "statistics.hpp"
#include "util/HighsTimer.h"

struct Runtime {
  Instance instance;
  Instance relaxed_for_ratiotest;
  Instance scaled;
  Instance perturbed;
  Settings settings;
  Statistics& statistics;

  QpVector primal;
  QpVector rowactivity;
  QpVector dualvar;
  QpVector dualcon;
  QpModelStatus status = QpModelStatus::kUndetermined;

  std::vector<BasisStatus> status_var;
  std::vector<BasisStatus> status_con;

  Runtime(Instance& inst, Statistics& stats)
      : instance(inst),
        statistics(stats),
        primal(QpVector(instance.num_var)),
        rowactivity(QpVector(instance.num_con)),
        dualvar(instance.num_var),
        dualcon(instance.num_con),
        status_var(instance.num_var),
        status_con(instance.num_con) {}
};

#endif

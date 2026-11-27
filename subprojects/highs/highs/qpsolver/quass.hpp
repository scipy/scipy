/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_QUASS_HPP__
#define __SRC_LIB_QUASS_HPP__

#include "lp_data/HighsCallback.h"
#include "qpsolver/basis.hpp"
#include "qpsolver/eventhandler.hpp"
#include "qpsolver/factor.hpp"
#include "qpsolver/instance.hpp"
#include "qpsolver/runtime.hpp"

struct Quass {
  Quass(Runtime& rt);

  void solve(const QpVector& x0, const QpVector& ra, Basis& b0,
             HighsTimer& timer, HighsCallback& callback);

 private:
  Runtime& runtime;
};

#endif

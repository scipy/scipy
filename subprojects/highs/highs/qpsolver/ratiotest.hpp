/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_RATIOTEST_HPP__
#define __SRC_LIB_RATIOTEST_HPP__

#include <limits>

#include "runtime.hpp"

struct RatiotestResult {
  double alpha;
  HighsInt limitingconstraint;
  bool nowactiveatlower;
};

RatiotestResult ratiotest(Runtime& runtime, const QpVector& p,
                          const QpVector& rowmove, double alphastart);

Instance ratiotest_relax_instance(Runtime& runtime);

#endif

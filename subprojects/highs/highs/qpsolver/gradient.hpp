/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_GRADIENT_HPP__
#define __SRC_LIB_GRADIENT_HPP__

#include "qpvector.hpp"
#include "runtime.hpp"

class Gradient {
  Runtime& runtime;

  QpVector gradient;
  bool uptodate;
  HighsInt numupdates = 0;

 public:
  Gradient(Runtime& rt)
      : runtime(rt), gradient(QpVector(rt.instance.num_var)), uptodate(false) {}

  void recompute() {
    runtime.instance.Q.vec_mat(runtime.primal, gradient);
    gradient += runtime.instance.c;
    uptodate = true;
    numupdates = 0;
  }

  QpVector& getGradient() {
    if (!uptodate ||
        numupdates >= runtime.settings.gradientrecomputefrequency) {
      recompute();
    }
    return gradient;
  }

  void update(QpVector& buffer_Qp, double stepsize) {
    gradient.saxpy(stepsize, buffer_Qp);
    numupdates++;
  }
};

#endif

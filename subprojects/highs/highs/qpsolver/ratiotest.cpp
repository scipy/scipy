/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "qpsolver/ratiotest.hpp"

static double step(double x, double p, double l, double u, double t) {
  if (p < -t && l > -kHighsInf) {
    return (l - x) / p;
  } else if (p > t && u < kHighsInf) {
    return (u - x) / p;
  } else {
    return kHighsInf;
  }
}

static RatiotestResult ratiotest_textbook(Runtime& rt, const QpVector& p,
                                          const QpVector& rowmove,
                                          Instance& instance,
                                          const double alphastart) {
  RatiotestResult result;
  result.limitingconstraint = -1;
  result.alpha = alphastart;
  result.nowactiveatlower = false;

  // check ratio towards variable bounds
  for (HighsInt j = 0; j < p.num_nz; j++) {
    HighsInt i = p.index[j];
    double alpha_i = step(rt.primal.value[i], p.value[i], instance.var_lo[i],
                          instance.var_up[i], rt.settings.ratiotest_t);
    if (alpha_i < result.alpha) {
      result.alpha = alpha_i;
      result.limitingconstraint = instance.num_con + i;
      result.nowactiveatlower = p.value[i] < 0;
    }
  }

  // check ratio towards constraint bounds
  for (HighsInt j = 0; j < rowmove.num_nz; j++) {
    HighsInt i = rowmove.index[j];
    double alpha_i =
        step(rt.rowactivity.value[i], rowmove.value[i], instance.con_lo[i],
             instance.con_up[i], rt.settings.ratiotest_t);
    if (alpha_i < result.alpha) {
      result.alpha = alpha_i;
      result.limitingconstraint = i;
      result.nowactiveatlower = rowmove.value[i] < 0;
    }
  }

  return result;
}

static RatiotestResult ratiotest_twopass(Runtime& runtime, const QpVector& p,
                                         const QpVector& rowmove,
                                         Instance& relaxed,
                                         const double alphastart) {
  RatiotestResult res1 =
      ratiotest_textbook(runtime, p, rowmove, relaxed, alphastart);

  if (res1.limitingconstraint == -1) {
    return res1;
  }

  RatiotestResult result = res1;

  double max_pivot = 0;
  if (res1.limitingconstraint != -1) {
    if ((int)result.limitingconstraint < runtime.instance.num_con) {
      max_pivot = rowmove.value[result.limitingconstraint];
    } else {
      max_pivot = p.value[result.limitingconstraint - runtime.instance.num_con];
    }
  }

  for (HighsInt i = 0; i < runtime.instance.num_con; i++) {
    double step_i = step(runtime.rowactivity.value[i], rowmove.value[i],
                         runtime.instance.con_lo[i], runtime.instance.con_up[i],
                         runtime.settings.ratiotest_t);
    if (fabs(rowmove.value[i]) >= fabs(max_pivot) && step_i <= res1.alpha) {
      max_pivot = rowmove.value[i];
      result.limitingconstraint = i;
      result.alpha = step_i;
      result.nowactiveatlower = rowmove.value[i] < 0;
    }
  }

  for (HighsInt i = 0; i < runtime.instance.num_var; i++) {
    double step_i =
        step(runtime.primal.value[i], p.value[i], runtime.instance.var_lo[i],
             runtime.instance.var_up[i], runtime.settings.ratiotest_t);
    if (fabs(p.value[i]) >= fabs(max_pivot) && step_i <= res1.alpha) {
      max_pivot = p.value[i];
      result.limitingconstraint = runtime.instance.num_con + i;
      result.alpha = step_i;
      result.nowactiveatlower = p.value[i] < 0;
    }
  }

  result.alpha = fmax(result.alpha, 0.0);
  return result;
}

Instance ratiotest_relax_instance(Runtime& runtime) {
  Instance relaxed_instance = runtime.instance;
  for (double& bound : relaxed_instance.con_lo) {
    if (bound != -kHighsInf) {
      bound -= runtime.settings.ratiotest_d;
    }
  }

  for (double& bound : relaxed_instance.con_up) {
    if (bound != kHighsInf) {
      bound += runtime.settings.ratiotest_d;
    }
  }

  for (double& bound : relaxed_instance.var_lo) {
    if (bound != -kHighsInf) {
      bound -= runtime.settings.ratiotest_d;
    }
  }

  for (double& bound : relaxed_instance.var_up) {
    if (bound != kHighsInf) {
      bound += runtime.settings.ratiotest_d;
    }
  }
  return relaxed_instance;
}

RatiotestResult ratiotest(Runtime& runtime, const QpVector& p,
                          const QpVector& rowmove, double alphastart) {
  switch (runtime.settings.ratiotest) {
    case RatiotestStrategy::Textbook:
      return ratiotest_textbook(runtime, p, rowmove, runtime.instance,
                                alphastart);
    case RatiotestStrategy::TwoPass:
    default:  // to fix -Wreturn-type warning
      return ratiotest_twopass(runtime, p, rowmove,
                               runtime.relaxed_for_ratiotest, alphastart);
  }
}

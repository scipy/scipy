/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "qpsolver/scaling.hpp"

#include <algorithm>
#include <map>

static double largestpoweroftwo(double value) {
  double l = log2(value);
  HighsInt il = (HighsInt)l;
  return powf(1.0, il);
}

static void scale_rows(Runtime& rt) {
  if (!rt.settings.rowscaling) {
    return;
  }

  std::map<HighsInt, double> maxabscoefperrow;
  for (HighsInt row = 0; row < rt.scaled.num_con; row++) {
    maxabscoefperrow[row] = 0.0;
  }

  for (HighsInt var = 0; var < rt.scaled.num_var; var++) {
    for (HighsInt index = rt.scaled.A.mat.start[var];
         index < rt.scaled.A.mat.start[var + 1]; index++) {
      if (fabs(rt.scaled.A.mat.value[index]) >
          maxabscoefperrow[rt.scaled.A.mat.index[index]]) {
        maxabscoefperrow[rt.scaled.A.mat.index[index]] =
            fabs(rt.scaled.A.mat.value[index]);
      }
    }
  }

  for (HighsInt var = 0; var < rt.scaled.num_var; var++) {
    for (HighsInt index = rt.scaled.A.mat.start[var];
         index < rt.scaled.A.mat.start[var + 1]; index++) {
      double factor =
          largestpoweroftwo(maxabscoefperrow[rt.scaled.A.mat.index[index]]);
      rt.scaled.A.mat.value[index] /= factor;
    }
  }

  for (HighsInt row = 0; row < rt.scaled.num_con; row++) {
    double factor = largestpoweroftwo(maxabscoefperrow[row]);
    if (rt.scaled.con_lo[row] > -kHighsInf) {
      rt.scaled.con_lo[row] /= factor;
    }
    if (rt.scaled.con_up[row] < kHighsInf) {
      rt.scaled.con_up[row] /= factor;
    }
  }
}

static void scale_cols(Runtime& rt) {
  if (!rt.settings.varscaling) {
    return;
  }

  std::map<HighsInt, double> maxabscoefpervar;
  for (HighsInt var = 0; var < rt.scaled.num_var; var++) {
    maxabscoefpervar[var] = 0.0;
  }

  for (HighsInt var = 0; var < rt.scaled.num_var; var++) {
    // gather information about variable
    for (HighsInt index = rt.scaled.A.mat.start[var];
         index < rt.scaled.A.mat.start[var + 1]; index++) {
      if (fabs(rt.scaled.A.mat.value[index]) > maxabscoefpervar[var]) {
        maxabscoefpervar[var] = fabs(rt.scaled.A.mat.value[index]);
      }
    }

    for (HighsInt index = rt.scaled.Q.mat.start[var];
         index < rt.scaled.Q.mat.start[var + 1]; index++) {
      if (rt.scaled.Q.mat.index[index] == var) {
        maxabscoefpervar[var] = fmax(maxabscoefpervar[var],
                                     sqrt(fabs(rt.scaled.Q.mat.value[index])));
      }
    }
  }

  std::map<HighsInt, double> factorpervar;

  for (HighsInt var = 0; var < rt.scaled.num_var; var++) {
    factorpervar[var] = fmin(128.0, largestpoweroftwo(maxabscoefpervar[var]));
  }

  for (HighsInt var = 0; var < rt.scaled.num_var; var++) {
    // scale variable: constraint matrix
    double factor = factorpervar[var];
    for (HighsInt index = rt.scaled.A.mat.start[var];
         index < rt.scaled.A.mat.start[var + 1]; index++) {
      rt.scaled.A.mat.value[index] /= factor;
    }

    // scale variable: linear objective
    rt.scaled.c.value[var] /= factor;
  }

  // scale variable: hessian matrix
  for (HighsInt var1 = 0; var1 < rt.scaled.num_var; var1++) {
    double factor1 = factorpervar[var1];
    for (HighsInt index = rt.scaled.Q.mat.start[var1];
         index < rt.scaled.Q.mat.start[var1 + 1]; index++) {
      HighsInt var2 = rt.scaled.Q.mat.index[index];
      double factor2 = factorpervar[var2];
      rt.scaled.Q.mat.value[index] /= (factor1 * factor2);
    }
  }
}

void scale(Runtime& rt) {
  rt.scaled = rt.instance;
  scale_rows(rt);
  scale_cols(rt);
  scale_rows(rt);
}

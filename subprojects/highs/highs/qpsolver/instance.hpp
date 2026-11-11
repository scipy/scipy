/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_INSTANCE_HPP__
#define __SRC_LIB_INSTANCE_HPP__

#include <vector>

#include "matrix.hpp"
#include "qpvector.hpp"

struct SumNum {
  double sum = 0.0;
  HighsInt num = 0;
};

struct Instance {
  HighsInt sense = 1;  // Minimization
  HighsInt num_var = 0;
  HighsInt num_con = 0;
  double offset = 0;
  QpVector c;
  Matrix Q;
  std::vector<double> con_lo;
  std::vector<double> con_up;
  Matrix A;
  std::vector<double> var_lo;
  std::vector<double> var_up;

  Instance(HighsInt nv = 0, HighsInt nc = 0)
      : num_var(nv),
        num_con(nc),
        c(QpVector(nv)),
        Q(Matrix(nv, nv)),
        A(Matrix(nc, nv)) {}

  double objval(const QpVector& x) {
    return c * x + 0.5 * (Q.vec_mat(x) * x) + offset;
  }

  bool isEquality(const HighsInt iVar) {
    return iVar < this->num_con ? this->con_lo[iVar] == this->con_up[iVar]
                                : this->var_lo[iVar - this->num_con] ==
                                      this->var_up[iVar - this->num_con];
  }
  SumNum sumnumprimalinfeasibilities(const QpVector& x,
                                     const QpVector& rowactivity) {
    SumNum res;
    for (HighsInt row = 0; row < num_con; row++) {
      if (rowactivity.value[row] < con_lo[row]) {
        res.sum += (con_lo[row] - rowactivity.value[row]);
        res.num++;
      } else if (rowactivity.value[row] > con_up[row]) {
        res.sum += (rowactivity.value[row] - con_up[row]);
        res.num++;
      }
    }
    for (HighsInt var = 0; var < num_var; var++) {
      if (x.value[var] < var_lo[var]) {
        res.sum += (var_lo[var] - x.value[var]);
        res.num++;
      } else if (x.value[var] > var_up[var]) {
        res.sum += (x.value[var] - var_up[var]);
        res.num++;
      }
    }
    return res;
  }
};

#endif

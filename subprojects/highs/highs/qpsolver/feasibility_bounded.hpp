/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_FEASIBILITYBOUNDED_HPP__
#define __SRC_LIB_FEASIBILITYBOUNDED_HPP__

#include "Highs.h"
#include "qpsolver/a_asm.hpp"

static void computeStartingPointBounded(Instance& instance, Settings& settings,
                                        Statistics& stats,
                                        QpModelStatus& modelstatus,
                                        QpHotstartInformation& result,
                                        HighsTimer& timer) {
  // compute initial feasible point for problems with bounds only (no general
  // linear constraints)
  const bool debug_printing = false;
  // Solve  Qx + c = 0 --> x = -Q^-1c
  HighsInt dim = instance.num_var;
  QpVector res = -instance.c;
  MatrixBase& hessian = instance.Q.mat;
  assert(res.dim == dim);
  if (hessian.isDiagonal()) {
    // Diagonal Hessian
    for (HighsInt iRow = 0; iRow < dim; iRow++) {
      double value = hessian.diagonal(iRow);
      if (value <= 0) {
        modelstatus = QpModelStatus::kNonConvex;
        return;
      }
      res.value[iRow] /= value;
    }
  } else {
    // General Hessian: compute Cholesky factorization of Q
    std::vector<double> L;
    L.resize(dim * dim);

    auto printL = [&]() {
      if (!debug_printing) return;
      for (HighsInt iRow = 0; iRow < dim; iRow++) {
        for (HighsInt iCol = 0; iCol <= iRow; iCol++)
          printf(" %11.4g", L[iRow * dim + iCol]);
        printf("\n");
      }
    };

    // First copy the lower triangle of Q into L
    L.assign(dim * dim, 0);
    HighsInt num_entries;
    std::vector<double> value(dim);
    std::vector<HighsInt> index(dim);
    for (HighsInt iCol = 0; iCol < dim; iCol++) {
      hessian.getColumn(iCol, num_entries, index.data(), value.data());
      for (HighsInt iEl = 0; iEl < num_entries; iEl++) {
        HighsInt iRow = index[iEl];
        // Take the entries above or on the diagonal in column iCol of
        // Q as the entries before or on the diagonal in row iCol of
        // (row-wise) L
        if (iRow <= iCol) L[iCol * dim + iRow] = value[iEl];
      }
    }
    // Now compute Cholesky factorization of L
    for (HighsInt iRow = 0; iRow < dim; iRow++) {
      printL();
      for (HighsInt iCol = 0; iCol <= iRow; iCol++) {
        double sum = 0;
        for (HighsInt k = 0; k < iCol; k++)
          sum += L[iRow * dim + k] * L[iCol * dim + k];
        if (iCol < iRow) {
          double value = (L[iRow * dim + iCol] - sum) / L[iCol * dim + iCol];
          if (debug_printing)
            printf(
                "Computed L[%1d, %1d] = (%11.4g - %11.4g) / %11.4g = %11.4g\n",
                int(iRow), int(iCol), L[iRow * dim + iCol], sum,
                L[iCol * dim + iCol], value);
          L[iRow * dim + iCol] = value;
        } else {
          double value = L[iCol * dim + iCol] - sum;
          if (debug_printing)
            printf(
                "Computed L[%1d, %1d] = sqrt(%11.4g - %11.4g) = sqrt(%11.4g) = "
                "%11.4g\n",
                int(iCol), int(iCol), L[iCol * dim + iCol], sum, value,
                std::sqrt(value));
          if (value <= 0) {
            modelstatus = QpModelStatus::kNonConvex;
            return;
          }
          L[iCol * dim + iCol] = std::sqrt(value);
        }
      }
    }
    printL();
    // Solve for Qx = -c
    // Solve Ly = -c
    for (HighsInt iRow = 0; iRow < dim; iRow++) {
      double sum = 0.0;
      for (HighsInt iCol = 0; iCol < iRow; iCol++)
        sum += res.value[iCol] * L[iRow * dim + iCol];
      res.value[iRow] = (res.value[iRow] - sum) / L[iRow * dim + iRow];
    }
    // Solve L^Tx = y
    for (HighsInt iRow = dim - 1; iRow >= 0; iRow--) {
      res.value[iRow] /= L[iRow * dim + iRow];
      for (HighsInt iCol = 0; iCol < iRow; iCol++)
        res.value[iCol] -= res.value[iRow] * L[iRow * dim + iCol];
    }
  }

  // project solution to bounds and collect active bounds
  QpVector x0(instance.num_var);
  QpVector ra(instance.num_con);
  std::vector<HighsInt> initialactive;
  std::vector<HighsInt> initialinactive;
  std::vector<BasisStatus> atlower;

  for (int i = 0; i < instance.num_var; i++) {
    // Check for unboundedness so that qp-unbounded and test-qod pass:
    // assumes diagonal Hessian and
    // settings.hessian_regularization_value > 0
    if (settings.hessian_regularization_value > 0) {
      if (res.value[i] > 0.5 / settings.hessian_regularization_value &&
          instance.var_up[i] == kHighsInf && instance.c.value[i] < 0.0) {
        modelstatus = QpModelStatus::kUnbounded;
        return;
      } else if (res.value[i] < -0.5 / settings.hessian_regularization_value &&
                 instance.var_lo[i] == -kHighsInf &&
                 instance.c.value[i] > 0.0) {
        modelstatus = QpModelStatus::kUnbounded;
        return;
      }
    }
    if (res.value[i] <
        instance.var_lo[i] - settings.primal_feasibility_tolerance) {
      res.value[i] = instance.var_lo[i];
      initialactive.push_back(i + instance.num_con);
      atlower.push_back(BasisStatus::kActiveAtLower);
    } else if (res.value[i] >
               instance.var_up[i] + settings.primal_feasibility_tolerance) {
      res.value[i] = instance.var_up[i];
      initialactive.push_back(i + instance.num_con);
      atlower.push_back(BasisStatus::kActiveAtUpper);
    } else {
      initialinactive.push_back(i + instance.num_con);
    }
    if (fabs(res.value[i]) > 1e-4) {
      x0.value[i] = res.value[i];
      x0.index[x0.num_nz++] = i;
    }
  }

  // if no bounds are active, solution lies in the interior -> optimal
  if (initialactive.size() == 0) {
    modelstatus = QpModelStatus::kOptimal;
  }

  assert((HighsInt)(initialactive.size() + initialinactive.size()) ==
         instance.num_var);

  result.status = atlower;
  result.active = initialactive;
  result.inactive = initialinactive;
  result.primal = x0;
  result.rowact = ra;
}

#endif

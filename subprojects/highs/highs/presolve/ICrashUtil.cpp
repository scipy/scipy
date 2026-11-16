/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/ICrashUtil.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "presolve/ICrashUtil.h"

#include <algorithm>
#include <cmath>

#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsSolution.h"
#include "util/HighsUtils.h"

void convertToMinimization(HighsLp& lp) {
  if (lp.sense_ != ObjSense::kMinimize) {
    for (int col = 0; col < lp.num_col_; col++)
      lp.col_cost_[col] = -lp.col_cost_[col];
  }
}

bool isEqualityProblem(const HighsLp& lp) {
  for (int row = 0; row < lp.num_row_; row++)
    if (lp.row_lower_[row] != lp.row_upper_[row]) return false;

  return true;
}

double vectorProduct(const std::vector<double>& v1,
                     const std::vector<double>& v2) {
  assert(v1.size() == v2.size());
  double sum = 0;
  for (int i = 0; i < (int)v1.size(); i++) sum += v1[i] * v2[i];
  return sum;
}

void muptiplyByTranspose(const HighsLp& lp, const std::vector<double>& v,
                         std::vector<double>& result) {
  assert((int)result.size() == lp.num_col_);
  assert((int)v.size() == lp.num_row_);

  result.assign(lp.num_col_, 0);
  for (int col = 0; col < lp.num_col_; col++) {
    for (int k = lp.a_matrix_.start_[col]; k < lp.a_matrix_.start_[col + 1];
         k++) {
      const int row = lp.a_matrix_.index_[k];
      result.at(col) += lp.a_matrix_.value_[k] * v[row];
    }
  }
}

void printMinorIterationDetails(const double iteration, const double col,
                                const double old_value, const double update,
                                const double ctx, const std::vector<double>& r,
                                const double quadratic_objective,
                                const HighsLogOptions& options) {
  double rnorm = getNorm2(r);
  std::stringstream ss;
  ss << "iter " << iteration;
  ss << ", col " << col;
  ss << ", update " << update;
  ss << ", old_value " << old_value;
  ss << ", new_value " << old_value + update;
  ss << ", ctx " << ctx;
  ss << ", r " << rnorm;
  ss << ", quadratic_objective " << quadratic_objective;
  ss << std::endl;

  // std::cout << ss.str();
  highsLogUser(options, HighsLogType::kInfo, ss.str().c_str());
}

bool initialize(const HighsLp& lp, HighsSolution& solution,
                std::vector<double>& lambda) {
  // Clear and resize primal column solution.
  solution.clear();
  solution.col_value.resize(lp.num_col_);

  for (int col = 0; col < lp.num_col_; col++) {
    if (lp.col_lower_[col] <= 0 && lp.col_upper_[col] >= 0)
      solution.col_value[col] = 0;
    else if (lp.col_lower_[col] > 0)
      solution.col_value[col] = lp.col_lower_[col];
    else if (lp.col_upper_[col] < 0)
      solution.col_value[col] = lp.col_upper_[col];
    else {
      printf("ICrash error: setting initial value for column %d\n", col);
      return false;
    }
  }

  lambda.resize(lp.num_row_);
  lambda.assign(lp.num_row_, 0);

  return true;
}

double minimizeComponentQP(const int col, const double mu, const HighsLp& lp,
                           double& objective, std::vector<double>& residual,
                           HighsSolution& sol) {
  // Minimize quadratic for column col.

  // Formulas for a and b when minimizing for x_j
  // a = (1/(2*mu)) * sum_i a_ij^2
  // b = -(1/(2*mu)) sum_i (2 * a_ij * (sum_{k!=j} a_ik * x_k - b_i))
  // b / 2 = -(1/(2*mu)) sum_i (2 * a_ij
  double a = 0.0;
  double b = 0.0;

  for (int k = lp.a_matrix_.start_[col]; k < lp.a_matrix_.start_[col + 1];
       k++) {
    int row = lp.a_matrix_.index_[k];
    a += lp.a_matrix_.value_[k] * lp.a_matrix_.value_[k];
    // matlab but with b = b / 2
    double bracket =
        -residual[row] - lp.a_matrix_.value_[k] * sol.col_value[col];
    // clp minimizing for delta_x
    // double bracket_clp = - residual_[row];
    b += lp.a_matrix_.value_[k] * bracket;
  }

  a = (0.5 / mu) * a;
  b = (0.5 / mu) * b + 0.5 * lp.col_cost_[col];

  double theta = -b / a;
  double delta_x = 0;

  // matlab
  double new_x;
  if (theta > 0)
    new_x = std::min(theta, lp.col_upper_[col]);
  else
    new_x = std::max(theta, lp.col_lower_[col]);
  delta_x = new_x - sol.col_value[col];

  // clp minimizing for delta_x
  // if (theta > 0)
  //   delta_x = std::min(theta, lp_.col_upper_[col] - col_value_[col]);
  // else
  //   delta_x = std::max(theta, lp_.col_lower_[col] - col_value_[col]);

  sol.col_value[col] += delta_x;

  // std::cout << "col " << col << ": " << delta_x << std::endl;

  // Update objective, row_value, residual after each component update.
  objective += lp.col_cost_[col] * delta_x;
  for (int k = lp.a_matrix_.start_[col]; k < lp.a_matrix_.start_[col + 1];
       k++) {
    int row = lp.a_matrix_.index_[k];
    sol.row_value[row] += lp.a_matrix_.value_[k] * delta_x;
    residual[row] = std::fabs(lp.row_upper_[row] - sol.row_value[row]);
  }

  return delta_x;
}

double minimizeComponentIca(const int col, const double mu,
                            const std::vector<double>& lambda,
                            const HighsLp& lp, double& objective,
                            std::vector<double>& residual, HighsSolution& sol) {
  // Minimize quadratic for column col.

  // Formulas for a and b when minimizing for x_j
  // a = (1/(2*mu)) * sum_i a_ij^2
  // b = -(1/(2*mu)) sum_i (2 * a_ij * (sum_{k!=j} a_ik * x_k - b_i)) + c_j (\)
  //     + sum_i a_ij * lambda_i
  // b / 2 = -(1/(2*mu)) sum_i (2 * a_ij
  double a = 0.0;
  double b = 0.0;

  for (int k = lp.a_matrix_.start_[col]; k < lp.a_matrix_.start_[col + 1];
       k++) {
    int row = lp.a_matrix_.index_[k];
    a += lp.a_matrix_.value_[k] * lp.a_matrix_.value_[k];
    // matlab but with b = b / 2
    double bracket =
        -residual[row] - lp.a_matrix_.value_[k] * sol.col_value[col];
    bracket += lambda[row];
    // clp minimizing for delta_x
    // double bracket_clp = - residual_[row];
    b += lp.a_matrix_.value_[k] * bracket;
  }

  a = (0.5 / mu) * a;
  b = (0.5 / mu) * b + 0.5 * lp.col_cost_[col];

  double theta = -b / a;
  double delta_x = 0;

  // matlab
  double new_x;
  if (theta > 0)
    new_x = std::min(theta, lp.col_upper_[col]);
  else
    new_x = std::max(theta, lp.col_lower_[col]);
  delta_x = new_x - sol.col_value[col];

  // clp minimizing for delta_x
  // if (theta > 0)
  //   delta_x = std::min(theta, lp_.col_upper_[col] - col_value_[col]);
  // else
  //   delta_x = std::max(theta, lp_.col_lower_[col] - col_value_[col]);

  sol.col_value[col] += delta_x;

  // std::cout << "col " << col << ": " << delta_x << std::endl;

  // Update objective, row_value, residual after each component update.
  objective += lp.col_cost_[col] * delta_x;
  for (int k = lp.a_matrix_.start_[col]; k < lp.a_matrix_.start_[col + 1];
       k++) {
    int row = lp.a_matrix_.index_[k];
    residual[row] -= lp.a_matrix_.value_[k] * delta_x;
    sol.row_value[row] += lp.a_matrix_.value_[k] * delta_x;
    // residual[row] = fabs(lp.row_upper_[row] - sol.row_value[row]); # ~~~
  }

  return delta_x;
}

void updateResidual(bool piecewise, const HighsLp& lp, const HighsSolution& sol,
                    std::vector<double>& residual) {
  residual.clear();
  residual.assign(lp.num_row_, 0);

  if (!piecewise) {
    assert(isEqualityProblem(lp));
    for (int row = 0; row < lp.num_row_; row++)
      residual[row] = std::fabs(lp.row_upper_[row] - sol.row_value[row]);
  } else {
    // piecewise
    for (int row = 0; row < lp.num_row_; row++) {
      double value = 0;
      if (sol.row_value[row] <= lp.row_lower_[row])
        value = lp.row_lower_[row] - sol.row_value[row];
      else if (sol.row_value[row] >= lp.row_upper_[row])
        value = sol.row_value[row] - lp.row_upper_[row];

      residual[row] = value;
    }
  }
}

void updateResidualFast(const HighsLp& lp, const HighsSolution& sol,
                        std::vector<double>& residual) {
  assert(isEqualityProblem(lp));
  for (int row = 0; row < lp.num_row_; row++) {
    residual[row] = std::fabs(lp.row_upper_[row] - sol.row_value[row]);
  }
}

// Allows negative residuals
void updateResidualIca(const HighsLp& lp, const HighsSolution& sol,
                       std::vector<double>& residual) {
  assert(isEqualityProblem(lp));
  for (int row = 0; row < lp.num_row_; row++)
    residual[row] = lp.row_upper_[row] - sol.row_value[row];
}

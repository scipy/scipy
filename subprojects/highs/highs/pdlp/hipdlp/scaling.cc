/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/scaling.cc
 * @brief
 */
#include "scaling.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "io/HighsIO.h"  // For pdlpLogging
#include "linalg.hpp"

void Scaling::initialize(const HighsLp& lp) {
  col_scale_.assign(lp.num_col_, 1.0);
  row_scale_.assign(lp.num_row_, 1.0);
  is_scaled_ = false;

  // use linalg to compute norms
  norm_cost_ = linalg::computeCostNorm(lp, 2.0);
  norm_rhs_ = linalg::computeRhsNorm(lp, 2.0);
}

void Scaling::scaleProblem() {
  is_scaled_ = false;

  if (params_->use_ruiz_scaling) {
    highsLogDev(params_->log_options_, HighsLogType::kInfo,
                "Applying Ruiz scaling...\n");
    applyRuizScaling();
    is_scaled_ = true;
  }

  if (params_->use_pc_scaling) {
    highsLogDev(params_->log_options_, HighsLogType::kInfo,
                "Applying Pock-Chambolle scaling...\n");
    applyPockChambolleScaling();
    is_scaled_ = true;
  }

  if (params_->use_l2_scaling) {
    highsLogDev(params_->log_options_, HighsLogType::kInfo,
                "Applying L2 norm scaling...\n");
    applyL2Scaling();
    is_scaled_ = true;
  }
}

void Scaling::applyRuizScaling() {
  std::vector<double> current_col_scaling(lp_->num_col_);
  std::vector<double> current_row_scaling(lp_->num_row_);

  for (HighsInt iter = 0; iter < params_->ruiz_iterations; ++iter) {
    // Reset current scaling factors
    std::fill(current_col_scaling.begin(), current_col_scaling.end(), 0.0);
    std::fill(current_row_scaling.begin(), current_row_scaling.end(), 0.0);

    // Compute column norms (norm of each column)
    for (HighsInt col = 0; col < lp_->num_col_; ++col) {
      std::vector<double> col_values;
      for (HighsInt el = lp_->a_matrix_.start_[col];
           el < lp_->a_matrix_.start_[col + 1]; ++el) {
        col_values.push_back(lp_->a_matrix_.value_[el]);
      }

      if (!col_values.empty()) {
        current_col_scaling[col] = std::sqrt(computeNorm(
            col_values.data(), col_values.size(), params_->ruiz_norm));
      }

      if (current_col_scaling[col] == 0.0) {
        current_col_scaling[col] = 1.0;
      }
    }

    // Compute row norms (infinity norm for rows)
    const bool ruiz_norm_ok = params_->ruiz_norm == INFINITY;
    if (!ruiz_norm_ok) {
      highsLogUser(params_->log_options_, HighsLogType::kError,
                   "Currently only support infinity norm for Ruiz scaling\n");
      assert(ruiz_norm_ok);
    }
    if (params_->ruiz_norm == INFINITY) {
      // For infinity norm, find max absolute value in each row
      for (HighsInt col = 0; col < lp_->num_col_; ++col) {
        for (HighsInt el = lp_->a_matrix_.start_[col];
             el < lp_->a_matrix_.start_[col + 1]; ++el) {
          HighsInt row = lp_->a_matrix_.index_[el];
          double abs_val = std::abs(lp_->a_matrix_.value_[el]);
          current_row_scaling[row] =
              std::max(current_row_scaling[row], abs_val);
        }
      }

      for (HighsInt row = 0; row < lp_->num_row_; ++row) {
        if (current_row_scaling[row] == 0.0) {
          current_row_scaling[row] = 1.0;
        } else {
          current_row_scaling[row] = std::sqrt(current_row_scaling[row]);
        }
      }
    }

    // Apply the scaling
    applyScaling(current_col_scaling, current_row_scaling);

    // Update cumulative scaling factors
    for (HighsInt i = 0; i < lp_->num_col_; ++i) {
      col_scale_[i] *= current_col_scaling[i];
    }
    for (HighsInt i = 0; i < lp_->num_row_; ++i) {
      row_scale_[i] *= current_row_scaling[i];
    }
  }
}

void Scaling::applyPockChambolleScaling() {
  const bool pc_alpha_ok = 0.0 <= params_->pc_alpha && params_->pc_alpha <= 2.0;
  if (!pc_alpha_ok) {
    highsLogUser(params_->log_options_, HighsLogType::kError,
                 "PC alpha must be in [0, 2]\n");
    assert(pc_alpha_ok);
  }

  std::vector<double> current_col_scaling(lp_->num_col_, 0.0);
  std::vector<double> current_row_scaling(lp_->num_row_, 0.0);

  // Compute column scaling: (sum |A_ij|^alpha)^(1/alpha)
  for (HighsInt col = 0; col < lp_->num_col_; ++col) {
    for (HighsInt el = lp_->a_matrix_.start_[col];
         el < lp_->a_matrix_.start_[col + 1]; ++el) {
      current_col_scaling[col] +=
          std::pow(std::abs(lp_->a_matrix_.value_[el]), params_->pc_alpha);
    }

    if (current_col_scaling[col] > 0.0) {
      current_col_scaling[col] = std::sqrt(
          std::pow(current_col_scaling[col], 1.0 / params_->pc_alpha));
    } else {
      current_col_scaling[col] = 1.0;
    }
  }

  // Compute row scaling: (sum |A_ij|^(2-alpha))^(1/(2-alpha))
  for (HighsInt col = 0; col < lp_->num_col_; ++col) {
    for (HighsInt el = lp_->a_matrix_.start_[col];
         el < lp_->a_matrix_.start_[col + 1]; ++el) {
      HighsInt row = lp_->a_matrix_.index_[el];
      current_row_scaling[row] += std::pow(std::abs(lp_->a_matrix_.value_[el]),
                                           2.0 - params_->pc_alpha);
    }
  }

  for (HighsInt row = 0; row < lp_->num_row_; ++row) {
    if (current_row_scaling[row] > 0.0) {
      current_row_scaling[row] = std::sqrt(
          std::pow(current_row_scaling[row], 1.0 / (2.0 - params_->pc_alpha)));
    } else {
      current_row_scaling[row] = 1.0;
    }
  }

  // Apply the scaling
  applyScaling(current_col_scaling, current_row_scaling);

  // Update cumulative scaling factors
  for (HighsInt i = 0; i < lp_->num_col_; ++i) {
    col_scale_[i] *= current_col_scaling[i];
  }
  for (HighsInt i = 0; i < lp_->num_row_; ++i) {
    row_scale_[i] *= current_row_scaling[i];
  }
}

void Scaling::applyL2Scaling() {
  std::vector<double> current_col_scaling(lp_->num_col_, 1.0);
  std::vector<double> current_row_scaling(lp_->num_row_, 0.0);

  // Compute L2 norm of each column
  for (HighsInt col = 0; col < lp_->num_col_; ++col) {
    double sum_sq = 0.0;
    for (HighsInt el = lp_->a_matrix_.start_[col];
         el < lp_->a_matrix_.start_[col + 1]; ++el) {
      sum_sq += lp_->a_matrix_.value_[el] * lp_->a_matrix_.value_[el];
    }

    if (sum_sq > 0.0) {
      current_col_scaling[col] = std::sqrt(std::sqrt(sum_sq));
    } else {
      current_col_scaling[col] = 1.0;
    }
  }

  // Compute L2 norm of each row
  for (HighsInt col = 0; col < lp_->num_col_; ++col) {
    for (HighsInt el = lp_->a_matrix_.start_[col];
         el < lp_->a_matrix_.start_[col + 1]; ++el) {
      HighsInt row = lp_->a_matrix_.index_[el];
      current_row_scaling[row] +=
          lp_->a_matrix_.value_[el] * lp_->a_matrix_.value_[el];
    }
  }

  for (HighsInt row = 0; row < lp_->num_row_; ++row) {
    if (current_row_scaling[row] > 0.0) {
      current_row_scaling[row] = std::sqrt(std::sqrt(current_row_scaling[row]));
    } else {
      current_row_scaling[row] = 1.0;
    }
  }

  // Apply the scaling
  applyScaling(current_col_scaling, current_row_scaling);

  // Update cumulative scaling factors
  for (HighsInt i = 0; i < lp_->num_col_; ++i) {
    col_scale_[i] *= current_col_scaling[i];
  }
  for (HighsInt i = 0; i < lp_->num_row_; ++i) {
    row_scale_[i] *= current_row_scaling[i];
  }
}

void Scaling::applyScaling(const std::vector<double>& col_scaling,
                           const std::vector<double>& row_scaling) {
  // Scale cost vector: c_scaled = c / col_scaling
  for (HighsInt i = 0; i < lp_->num_col_; ++i) {
    lp_->col_cost_[i] /= col_scaling[i];
  }

  // Scale column bounds: l_scaled = l * col_scaling, u_scaled = u * col_scaling
  for (HighsInt i = 0; i < lp_->num_col_; ++i) {
    if (lp_->col_lower_[i] > -kHighsInf) {
      lp_->col_lower_[i] *= col_scaling[i];
    }
    if (lp_->col_upper_[i] < kHighsInf) {
      lp_->col_upper_[i] *= col_scaling[i];
    }
  }

  // Scale row bounds: b_scaled = b / row_scaling
  for (HighsInt i = 0; i < lp_->num_row_; ++i) {
    if (lp_->row_lower_[i] > -kHighsInf) {
      lp_->row_lower_[i] /= row_scaling[i];
    }
    if (lp_->row_upper_[i] < kHighsInf) {
      lp_->row_upper_[i] /= row_scaling[i];
    }
  }

  // Scale matrix: A_scaled = diag(1/row_scaling) * A * diag(1/col_scaling)
  for (HighsInt col = 0; col < lp_->num_col_; ++col) {
    for (HighsInt el = lp_->a_matrix_.start_[col];
         el < lp_->a_matrix_.start_[col + 1]; ++el) {
      HighsInt row = lp_->a_matrix_.index_[el];
      lp_->a_matrix_.value_[el] /= (row_scaling[row] * col_scaling[col]);
    }
  }

  // TO DO: Scale offset
}

void Scaling::unscaleSolution(std::vector<double>& x,
                              std::vector<double>& y) const {
  if (!is_scaled_) return;

  // Unscale primal variables: x_original = x_scaled / col_scale
  for (size_t i = 0; i < x.size(); ++i) {
    x[i] /= col_scale_[i];
  }

  // Unscale dual variables: y_original = y_scaled / row_scale
  for (size_t i = 0; i < y.size(); ++i) {
    y[i] /= row_scale_[i];
  }
}

double Scaling::computeNorm(const double* values, HighsInt size,
                            double norm_type) const {
  if (norm_type == INFINITY) {
    double max_val = 0.0;
    for (HighsInt i = 0; i < size; ++i) {
      max_val = std::max(max_val, std::abs(values[i]));
    }
    return max_val;
  } else if (norm_type == 2.0) {
    double sum_sq = 0.0;
    for (HighsInt i = 0; i < size; ++i) {
      sum_sq += values[i] * values[i];
    }
    return std::sqrt(sum_sq);
  } else {
    // General p-norm
    double sum = 0.0;
    for (HighsInt i = 0; i < size; ++i) {
      sum += std::pow(std::abs(values[i]), norm_type);
    }
    return std::pow(sum, 1.0 / norm_type);
  }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/linalg.cc
 * @brief
 */
#include "linalg.hpp"

#include "Highs.h"

namespace linalg {

double projectBox(double x, double l, double u) {
  return std::max(l, std::min(x, u));
}

double projectNonNegative(double x) { return std::max(0.0, x); }

void projectBounds(const HighsLp& lp, std::vector<double>& x) {
  for (HighsInt i = 0; i < lp.num_col_; ++i) {
    // Project to upper bound
    if (x[i] > lp.col_upper_[i]) {
      x[i] = lp.col_upper_[i];
    }
    // Project to lower bound
    if (x[i] < lp.col_lower_[i]) {
      x[i] = lp.col_lower_[i];
    }
  }
}

void ax(const HighsLp& lp, const std::vector<double>& x,
        std::vector<double>& result) {
  std::fill(result.begin(), result.end(), 0.0);
  // Assumes column-wise matrix format
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    for (HighsInt i = lp.a_matrix_.start_[col];
         i < lp.a_matrix_.start_[col + 1]; ++i) {
      const HighsInt row = lp.a_matrix_.index_[i];
      result[row] += lp.a_matrix_.value_[i] * x[col];
    }
  }
}

void aTy(const HighsLp& lp, const std::vector<double>& y,
         std::vector<double>& result) {
  std::fill(result.begin(), result.end(), 0.0);
  // Assumes column-wise matrix format. For each column `col` of A,
  // this loop calculates Dot(column_col, y) and adds it to result[col].
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    for (HighsInt i = lp.a_matrix_.start_[col];
         i < lp.a_matrix_.start_[col + 1]; ++i) {
      const HighsInt row = lp.a_matrix_.index_[i];
      result[col] += lp.a_matrix_.value_[i] * y[row];
    }
  }
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument(
        "Vectors must be of the same size for Dot product.");
  }
  double result = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    result += a[i] * b[i];
  }
  return result;
}

// =================================================================
// Norm Functions
// =================================================================

double norm2(const std::vector<double>& vec) {
  return std::sqrt(dot(vec, vec));
}

double vectorNormSquared(const std::vector<double>& vec) {
  return dot(vec, vec);
}

double vectorNorm(const double* values, size_t size, double p) {
  if (std::isinf(p)) {  // Infinity norm
    double max_val = 0.0;
    for (size_t i = 0; i < size; ++i) {
      max_val = std::max(max_val, std::abs(values[i]));
    }
    return max_val;
  }
  if (p == 1.0) {  // L1 norm
    double sum = 0.0;
    for (size_t i = 0; i < size; ++i) {
      sum += std::abs(values[i]);
    }
    return sum;
  }
  // General p-norm (including L2)
  double sum = 0.0;
  for (size_t i = 0; i < size; ++i) {
    sum += std::pow(std::abs(values[i]), p);
  }
  return std::pow(sum, 1.0 / p);
}

double vectorNorm(const std::vector<double>& vec, double p) {
  // For L2 norm, call the optimized Norm2 to avoid pow()
  if (p == 2.0) {
    return norm2(vec);
  }
  return vectorNorm(vec.data(), vec.size(), p);
}

// =================================================================
// Vector Operations
// =================================================================

void scale(std::vector<double>& vec, double factor) {
  for (double& val : vec) {
    val *= factor;
  }
}

void normalize(std::vector<double>& vec) {
  double norm = norm2(vec);
  if (norm > 1e-9) {  // Use a small tolerance
    scale(vec, 1.0 / norm);
  }
}

double diffTwoNorm(const std::vector<double>& v1,
                   const std::vector<double>& v2) {
  if (v1.size() != v2.size()) {
    throw std::invalid_argument("Vectors must have the same dimension.");
  }
  double norm_sq = 0.0;
  for (size_t i = 0; i < v1.size(); ++i) {
    double diff = v1[i] - v2[i];
    norm_sq += diff * diff;
  }
  return std::sqrt(norm_sq);
}

// =================================================================
// LP-related Norm Computations
// =================================================================

double computeCostNorm(const HighsLp& lp, double p) {
  return vectorNorm(lp.col_cost_, p);
}

double computeRhsNorm(const HighsLp& lp, double p) {
  if (std::isinf(p)) {
    double max_val = 0.0;
    for (double val : lp.row_lower_) {
      if (std::isfinite(val)) max_val = std::max(max_val, std::abs(val));
    }
    return max_val;
  }

  double sum = 0.0;
  for (double val : lp.row_lower_) {
    if (std::isfinite(val)) {
      if (p == 1.0)
        sum += std::abs(val);
      else
        sum += std::pow(std::abs(val), p);
    }
  }

  if (p == 2.0) return std::sqrt(sum);
  if (p == 1.0) return sum;
  return std::pow(sum, 1.0 / p);
}

std::vector<double> computeColumnNorms(const HighsLp& lp, double p) {
  std::vector<double> col_norms(lp.num_col_, 0.0);
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    HighsInt start = lp.a_matrix_.start_[col];
    HighsInt end = lp.a_matrix_.start_[col + 1];
    if (start < end) {
      col_norms[col] = vectorNorm(&lp.a_matrix_.value_[start], end - start, p);
    }
  }
  return col_norms;
}

std::vector<double> computeRowNorms(const HighsLp& lp, double p) {
  std::vector<double> row_norms(lp.num_row_, 0.0);
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    for (HighsInt i = lp.a_matrix_.start_[col];
         i < lp.a_matrix_.start_[col + 1]; ++i) {
      const HighsInt row = lp.a_matrix_.index_[i];
      const double abs_val = std::abs(lp.a_matrix_.value_[i]);
      if (std::isinf(p)) {
        row_norms[row] = std::max(row_norms[row], abs_val);
      } else if (p == 1.0) {
        row_norms[row] += abs_val;
      } else {
        row_norms[row] += std::pow(abs_val, p);
      }
    }
  }

  if (p != 1.0 && !std::isinf(p)) {
    const double p_inv = 1.0 / p;
    for (HighsInt row = 0; row < lp.num_row_; ++row) {
      row_norms[row] = std::pow(row_norms[row], p_inv);
    }
  }
  return row_norms;
}

std::vector<double> vectorSubtract(const std::vector<double>& a,
                                   const std::vector<double>& b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument(
        "Vectors must be of the same size for subtraction.");
  }
  std::vector<double> result(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}

}  // namespace linalg
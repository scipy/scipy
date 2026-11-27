/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSolverSelect.cpp
 * @brief Choose an LP algorithm from features of the (presolved) model
 */
#include "util/HighsSolverSelect.h"

#include <algorithm>
#include <cmath>

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "util/HighsSolverSelectModels.h"

namespace {

// Tolerance for coefficients being "almost identical": see
// relativeNumberOfAlmostIdentical
const double kAlmostIdenticalTol = 1e-9;
// A row is "dense" when its nonzero count exceeds this multiple of the
// mean row length, and is also at least kDenseRowMinCount long
const double kDenseRowFactor = 10.0;
const double kDenseRowMinCount = 10.0;

// Fraction of `values` that lie in a cluster of size >= 2, where
// neighbouring sorted values join a cluster when they differ by at most
// tol * max(1, |value|). `values` is sorted in place.
double relativeNumberOfAlmostIdentical(std::vector<double>& values,
                                       const double tol) {
  const size_t n = values.size();
  if (n < 2) return 0.0;
  std::sort(values.begin(), values.end());
  size_t num_in_cluster = 0;
  size_t run = 1;
  for (size_t i = 1; i <= n; i++) {
    const bool close = i < n && (values[i] - values[i - 1]) <=
                                    tol * std::max(1.0, std::fabs(values[i]));
    if (close) {
      run++;
    } else {
      if (run >= 2) num_in_cluster += run;
      run = 1;
    }
  }
  return static_cast<double>(num_in_cluster) / static_cast<double>(n);
}

double safeRatio(const double numerator, const double denominator) {
  if (denominator <= 0.0) return 0.0;
  return numerator / denominator;
}

}  // namespace

HighsLpFeatures computeLpFeatures(const HighsLp& lp) {
  HighsLpFeatures f;
  const HighsInt num_col = lp.num_col_;
  const HighsInt num_row = lp.num_row_;
  const HighsSparseMatrix& matrix = lp.a_matrix_;
  const HighsInt num_nz = matrix.numNz();

  f.num_row = num_row;
  f.num_col = num_col;
  f.num_nz = num_nz;

  const double d_num_col = static_cast<double>(num_col);
  const double d_num_row = static_cast<double>(num_row);
  const double d_num_nz = static_cast<double>(num_nz);

  // --------------------------------------------------------------------------
  // Column bounds and objective coefficients.
  // --------------------------------------------------------------------------
  HighsInt num_cols_without_upper = 0;
  HighsInt num_cols_without_lower = 0;
  HighsInt num_free_cols = 0;
  HighsInt num_boxed_cols = 0;
  HighsInt num_singly_bounded_cols = 0;
  HighsInt num_fixed_cols = 0;

  double max_abs_obj = 0.0;
  double min_abs_obj = kHighsInf;
  std::vector<double> obj_values;
  obj_values.reserve(num_col);

  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    const bool has_lower = lp.col_lower_[iCol] > -kHighsInf;
    const bool has_upper = lp.col_upper_[iCol] < kHighsInf;
    if (!has_upper) num_cols_without_upper++;
    if (!has_lower) num_cols_without_lower++;
    if (!has_lower && !has_upper) {
      num_free_cols++;
    } else if (has_lower && has_upper) {
      if (lp.col_lower_[iCol] == lp.col_upper_[iCol])
        num_fixed_cols++;
      else
        num_boxed_cols++;
    } else {
      num_singly_bounded_cols++;
    }

    const double cost = lp.col_cost_[iCol];
    obj_values.push_back(cost);
    const double abs_cost = std::fabs(cost);
    if (abs_cost > 0.0) {
      max_abs_obj = std::max(max_abs_obj, abs_cost);
      min_abs_obj = std::min(min_abs_obj, abs_cost);
    }
  }
  if (min_abs_obj == kHighsInf) min_abs_obj = 0.0;

  // --------------------------------------------------------------------------
  // Row activity bounds ("right-hand sides").
  // --------------------------------------------------------------------------
  HighsInt num_equalities = 0;
  HighsInt num_inequalities = 0;
  HighsInt num_ranged_rows = 0;
  HighsInt num_free_rows = 0;

  double max_abs_rhs = 0.0;
  double min_abs_rhs = kHighsInf;
  std::vector<double> rhs_values;  // one representative per row
  rhs_values.reserve(num_row);

  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    const double lower = lp.row_lower_[iRow];
    const double upper = lp.row_upper_[iRow];
    const bool has_lower = lower > -kHighsInf;
    const bool has_upper = upper < kHighsInf;
    if (!has_lower && !has_upper) {
      num_free_rows++;
    } else if (has_lower && has_upper) {
      if (lower == upper)
        num_equalities++;
      else
        num_ranged_rows++;
    } else {
      num_inequalities++;
    }

    // Representative rhs: finite upper if present, else finite lower.
    if (has_upper)
      rhs_values.push_back(upper);
    else if (has_lower)
      rhs_values.push_back(lower);

    for (int k = 0; k < 2; k++) {
      const double bound = k == 0 ? lower : upper;
      const double abs_bound = std::fabs(bound);
      if (abs_bound > 0.0 && abs_bound < kHighsInf) {
        max_abs_rhs = std::max(max_abs_rhs, abs_bound);
        min_abs_rhs = std::min(min_abs_rhs, abs_bound);
      }
    }
  }
  if (min_abs_rhs == kHighsInf) min_abs_rhs = 0.0;

  // --------------------------------------------------------------------------
  // Single pass over the constraint matrix: per-row / per-column nonzero
  // counts and coefficient magnitude statistics.
  // --------------------------------------------------------------------------
  std::vector<HighsInt> row_count(num_row, 0);
  std::vector<HighsInt> col_count(num_col, 0);
  double max_abs_matrix = 0.0;
  double min_abs_matrix = kHighsInf;
  std::vector<double> abs_matrix_values;
  abs_matrix_values.reserve(num_nz);

  const bool colwise = matrix.isColwise();
  const HighsInt num_major = colwise ? num_col : num_row;
  const bool have_matrix =
      num_nz > 0 &&
      static_cast<HighsInt>(matrix.start_.size()) >= num_major + 1;
  if (have_matrix) {
    for (HighsInt major = 0; major < num_major; major++) {
      const HighsInt from = matrix.start_[major];
      const HighsInt to = matrix.start_[major + 1];
      if (colwise)
        col_count[major] = to - from;
      else
        row_count[major] = to - from;
      for (HighsInt el = from; el < to; el++) {
        const HighsInt minor = matrix.index_[el];
        if (colwise)
          row_count[minor]++;
        else
          col_count[minor]++;
        const double abs_value = std::fabs(matrix.value_[el]);
        abs_matrix_values.push_back(abs_value);
        if (abs_value > 0.0) {
          max_abs_matrix = std::max(max_abs_matrix, abs_value);
          min_abs_matrix = std::min(min_abs_matrix, abs_value);
        }
      }
    }
  }
  if (min_abs_matrix == kHighsInf) min_abs_matrix = 0.0;

  HighsInt max_row_count = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    max_row_count = std::max(max_row_count, row_count[iRow]);
  HighsInt max_col_count = 0;
  for (HighsInt iCol = 0; iCol < num_col; iCol++)
    max_col_count = std::max(max_col_count, col_count[iCol]);

  const double mean_row_count = safeRatio(d_num_nz, d_num_row);
  const double dense_row_threshold =
      std::max(kDenseRowFactor * mean_row_count, kDenseRowMinCount);
  HighsInt num_dense_rows = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    if (static_cast<double>(row_count[iRow]) > dense_row_threshold)
      num_dense_rows++;

  // --------------------------------------------------------------------------
  // Assemble the feature record.
  // --------------------------------------------------------------------------
  // How constrained is the problem?
  f.relative_num_equalities = safeRatio(num_equalities, d_num_row);
  f.relative_num_inequalities = safeRatio(num_inequalities, d_num_row);
  f.relative_num_ranged_rows = safeRatio(num_ranged_rows, d_num_row);
  f.relative_num_free_rows = safeRatio(num_free_rows, d_num_row);
  f.relative_num_cols_without_upper =
      safeRatio(num_cols_without_upper, d_num_col);
  f.relative_num_cols_without_lower =
      safeRatio(num_cols_without_lower, d_num_col);
  f.relative_num_free_cols = safeRatio(num_free_cols, d_num_col);
  f.relative_num_boxed_cols = safeRatio(num_boxed_cols, d_num_col);
  f.relative_num_singly_bounded_cols =
      safeRatio(num_singly_bounded_cols, d_num_col);
  f.relative_num_fixed_cols = safeRatio(num_fixed_cols, d_num_col);

  // How degenerate is the problem?
  f.relative_num_almost_identical_nonzeros =
      relativeNumberOfAlmostIdentical(abs_matrix_values, kAlmostIdenticalTol);
  f.relative_num_identical_obj_values =
      relativeNumberOfAlmostIdentical(obj_values, 0);
  f.relative_num_identical_rhs_values =
      relativeNumberOfAlmostIdentical(rhs_values, 0);

  // How stable is the problem?
  f.max_abs_matrix_coefficient = max_abs_matrix;
  f.min_abs_matrix_coefficient = min_abs_matrix;
  f.max_abs_obj_coefficient = max_abs_obj;
  f.min_abs_obj_coefficient = min_abs_obj;
  f.max_abs_rhs_value = max_abs_rhs;
  f.min_abs_rhs_value = min_abs_rhs;
  f.relative_max_matrix_coefficient = safeRatio(max_abs_matrix, min_abs_matrix);
  f.relative_max_obj_coefficient = safeRatio(max_abs_obj, min_abs_obj);
  f.relative_max_rhs_value = safeRatio(max_abs_rhs, min_abs_rhs);

  // How dense is the problem?
  f.matrix_nonzero_density = safeRatio(d_num_nz, d_num_row * d_num_col);
  f.avg_nonzeros_per_col = safeRatio(d_num_nz, d_num_col);
  f.avg_nonzeros_per_row = safeRatio(d_num_nz, d_num_row);
  f.relative_num_dense_rows = safeRatio(num_dense_rows, d_num_row);
  f.max_nonzeros_in_row = static_cast<double>(max_row_count);
  f.max_nonzeros_in_col = static_cast<double>(max_col_count);

  return f;
}

namespace {

// Name and member of each feature, in the order that the fitted models
// index
const std::pair<const char*, double HighsLpFeatures::*> kLpFeatures[] = {
    {"num_row", &HighsLpFeatures::num_row},
    {"num_col", &HighsLpFeatures::num_col},
    {"num_nz", &HighsLpFeatures::num_nz},
    {"num_integer_col", &HighsLpFeatures::num_integer_col},

    {"relative_num_equalities", &HighsLpFeatures::relative_num_equalities},
    {"relative_num_inequalities", &HighsLpFeatures::relative_num_inequalities},
    {"relative_num_ranged_rows", &HighsLpFeatures::relative_num_ranged_rows},
    {"relative_num_free_rows", &HighsLpFeatures::relative_num_free_rows},
    {"relative_num_cols_without_upper",
     &HighsLpFeatures::relative_num_cols_without_upper},
    {"relative_num_cols_without_lower",
     &HighsLpFeatures::relative_num_cols_without_lower},
    {"relative_num_free_cols", &HighsLpFeatures::relative_num_free_cols},
    {"relative_num_boxed_cols", &HighsLpFeatures::relative_num_boxed_cols},
    {"relative_num_singly_bounded_cols",
     &HighsLpFeatures::relative_num_singly_bounded_cols},
    {"relative_num_fixed_cols", &HighsLpFeatures::relative_num_fixed_cols},

    {"relative_num_almost_identical_nonzeros",
     &HighsLpFeatures::relative_num_almost_identical_nonzeros},
    {"relative_num_identical_obj_values",
     &HighsLpFeatures::relative_num_identical_obj_values},
    {"relative_num_identical_rhs_values",
     &HighsLpFeatures::relative_num_identical_rhs_values},

    {"relative_max_matrix_coefficient",
     &HighsLpFeatures::relative_max_matrix_coefficient},
    {"relative_max_obj_coefficient",
     &HighsLpFeatures::relative_max_obj_coefficient},
    {"relative_max_rhs_value", &HighsLpFeatures::relative_max_rhs_value},
    {"max_abs_matrix_coefficient",
     &HighsLpFeatures::max_abs_matrix_coefficient},
    {"min_abs_matrix_coefficient",
     &HighsLpFeatures::min_abs_matrix_coefficient},
    {"max_abs_obj_coefficient", &HighsLpFeatures::max_abs_obj_coefficient},
    {"min_abs_obj_coefficient", &HighsLpFeatures::min_abs_obj_coefficient},
    {"max_abs_rhs_value", &HighsLpFeatures::max_abs_rhs_value},
    {"min_abs_rhs_value", &HighsLpFeatures::min_abs_rhs_value},

    {"matrix_nonzero_density", &HighsLpFeatures::matrix_nonzero_density},
    {"avg_nonzeros_per_col", &HighsLpFeatures::avg_nonzeros_per_col},
    {"avg_nonzeros_per_row", &HighsLpFeatures::avg_nonzeros_per_row},
    {"relative_num_dense_rows", &HighsLpFeatures::relative_num_dense_rows},
    {"max_nonzeros_in_row", &HighsLpFeatures::max_nonzeros_in_row},
    {"max_nonzeros_in_col", &HighsLpFeatures::max_nonzeros_in_col},
};

}  // namespace

std::vector<std::pair<std::string, double>> highsLpFeatureVector(
    const HighsLpFeatures& f) {
  std::vector<std::pair<std::string, double>> feature;
  for (const auto& name_member : kLpFeatures)
    feature.emplace_back(name_member.first, f.*name_member.second);
  return feature;
}

std::vector<std::string> highsLpFeatureNames() {
  std::vector<std::string> names;
  for (const auto& name_member : kLpFeatures)
    names.push_back(name_member.first);
  return names;
}

static HighsSolverSelect selectSolverByHeuristic(const HighsLpFeatures& f) {
  // NOTE: provisional heuristic
  if (f.num_nz == 0 || f.num_row == 0 || f.num_col == 0)
    return HighsSolverSelect::kDualSimplex;

  const bool large = f.num_nz > 500000;
  const bool very_large = f.num_nz > 5000000;
  const bool dense = f.matrix_nonzero_density > 1e-2 ||
                     f.avg_nonzeros_per_row > 25.0 ||
                     f.relative_num_dense_rows > 1e-2;
  const bool ill_conditioned =
      f.relative_max_matrix_coefficient > 1e8 || f.relative_max_rhs_value > 1e8;

  // Simplex is the safe default for small, sparse, well-scaled models: it
  // thrives on sparsity and warm-starts cheaply.
  if (!large && !dense) return HighsSolverSelect::kDualSimplex;

  // Large and/or dense: an interior-point method amortises the factorisation
  // cost better.  Keep ill-conditioned models on IPX, whose linear algebra and
  // crossover are the more battle-tested; send the very largest well-scaled
  // models to HiPO.
  if (very_large && !ill_conditioned) return HighsSolverSelect::kHipo;
  return HighsSolverSelect::kIpx;
}

// Each solver's trees predict its regret, and the smallest is selected
static HighsSolverSelect evaluateModel(const HighsSolverSelectModel& model,
                                       const HighsLpFeatures& f) {
  if (f.num_nz < model.min_num_nz) return HighsSolverSelect::kDualSimplex;
  HighsSolverSelect best_solver = HighsSolverSelect::kDualSimplex;
  double best_regret = kHighsInf;
  for (HighsInt s = 0; s < model.num_solver; s++) {
    double sum = 0;
    for (HighsInt t = 0; t < model.num_tree; t++) {
      HighsInt node = model.root[s * model.num_tree + t];
      while (model.left[node] >= 0) {
        // The trees were fitted by scikit-learn, which compares in single
        // precision
        const float value =
            static_cast<float>(f.*kLpFeatures[model.feature[node]].second);
        node = value <= model.threshold[node] ? model.left[node]
                                              : model.right[node];
      }
      sum += model.value[node];
    }
    const double regret = model.initial[s] + model.learning_rate * sum;
    if (regret < best_regret) {
      best_regret = regret;
      best_solver = model.solver[s];
    }
  }
  return best_solver;
}

HighsSolverSelect selectSolver(const HighsLp& lp,
                               const SolverSelectStrategy strategy) {
  switch (strategy) {
    case SolverSelectStrategy::kSimplexOnly:
      return HighsSolverSelect::kDualSimplex;
    case SolverSelectStrategy::kColsRowsNnz: {
      // Only the size features are used, so don't compute the rest
      HighsLpFeatures f;
      f.num_row = lp.num_row_;
      f.num_col = lp.num_col_;
      f.num_nz = lp.a_matrix_.numNz();
      return evaluateModel(kS1Model, f);
    }
    case SolverSelectStrategy::kRankedFeatures:
      // The model chooses dual simplex below min_num_nz, so don't pay for
      // the features
      if (lp.a_matrix_.numNz() < kS2Model.min_num_nz)
        return HighsSolverSelect::kDualSimplex;
      return evaluateModel(kS2Model, computeLpFeatures(lp));
    case SolverSelectStrategy::kHeuristic:
      return selectSolverByHeuristic(computeLpFeatures(lp));
    case SolverSelectStrategy::kDualSimplex:
      return HighsSolverSelect::kDualSimplex;
    case SolverSelectStrategy::kIpx:
      return HighsSolverSelect::kIpx;
    case SolverSelectStrategy::kHipo:
      return HighsSolverSelect::kHipo;
    case SolverSelectStrategy::kIpxIfLarge:
    case SolverSelectStrategy::kHipoIfLarge:
      if (lp.a_matrix_.numNz() < kS2Model.min_num_nz)
        return HighsSolverSelect::kDualSimplex;
      return strategy == SolverSelectStrategy::kIpxIfLarge
                 ? HighsSolverSelect::kIpx
                 : HighsSolverSelect::kHipo;
  }
  return HighsSolverSelect::kDualSimplex;
}

#include "PreProcess.h"

#include "Iterate.h"
#include "Model.h"
#include "model/HighsHessianUtils.h"

namespace hipo {

void PreprocessorPoint::assertConsistency(Int n, Int m) const {
  assert(static_cast<Int>(x.size()) == n);
  assert(static_cast<Int>(xl.size()) == n);
  assert(static_cast<Int>(xu.size()) == n);
  assert(static_cast<Int>(slack.size()) == m);
  assert(static_cast<Int>(y.size()) == m);
  assert(static_cast<Int>(zl.size()) == n);
  assert(static_cast<Int>(zu.size()) == n);
}

void PreprocessEmptyRows::apply(Model& model) {
  Int& n = model.n_;
  Int& m = model.m_;
  HighsSparseMatrix& A = model.A_;
  std::vector<double>& b = model.b_;
  std::vector<char>& constraints = model.constraints_;

  n_pre = n;
  m_pre = m;
  b_pre = b;

  // find empty rows
  std::vector<Int> entries_per_row(m, 0);
  for (Int col = 0; col < n; ++col) {
    for (Int el = A.start_[col]; el < A.start_[col + 1]; ++el) {
      const Int row = A.index_[el];
      ++entries_per_row[row];
    }
  }

  empty_rows = 0;
  for (Int i : entries_per_row)
    if (i == 0) ++empty_rows;

  if (empty_rows > 0) {
    rows_shift.assign(m, 0);
    for (Int i = 0; i < m; ++i) {
      if (entries_per_row[i] == 0) {
        // count how many empty rows there are before a given row
        for (Int j = i + 1; j < m; ++j) ++rows_shift[j];
        rows_shift[i] = -1;
      }
    }

    // shift each row index by the number of empty rows before it
    for (Int col = 0; col < n; ++col) {
      for (Int el = A.start_[col]; el < A.start_[col + 1]; ++el) {
        const Int row = A.index_[el];
        A.index_[el] -= rows_shift[row];
      }
    }
    A.num_row_ -= empty_rows;

    // shift entries in b and constraints
    for (Int i = 0; i < m; ++i) {
      // ignore entries to be removed
      if (rows_shift[i] == -1) continue;

      Int shifted_pos = i - rows_shift[i];
      b[shifted_pos] = b[i];
      constraints[shifted_pos] = constraints[i];
    }
    b.resize(A.num_row_);
    constraints.resize(A.num_row_);

    m = A.num_row_;
  }

  n_post = n;
  m_post = m;
}

void PreprocessEmptyRows::undo(PreprocessorPoint& point, const Model& model,
                               const Iterate& it) const {
  point.assertConsistency(n_post, m_post);
  if (empty_rows > 0) {
    // Add Lagrange multiplier for empty rows that were removed
    // Add slack for constraints that were removed

    std::vector<double> new_y(m_pre, 0.0);
    std::vector<double> new_slack(m_pre, 0.0);

    // position to read from y and slack
    Int pos = 0;

    for (Int i = 0; i < m_pre; ++i) {
      // empty row
      if (rows_shift[i] == -1) {
        new_slack[i] = b_pre[i];
        continue;
      }

      // re-align value of y and slack, considering empty rows
      new_y[pos + rows_shift[i]] = point.y[pos];
      new_slack[pos + rows_shift[i]] = point.slack[pos];
      ++pos;
    }

    point.y = std::move(new_y);
    point.slack = std::move(new_slack);
  }
  point.assertConsistency(n_pre, m_pre);
}

void PreprocessEmptyRows::print(std::stringstream& stream) const {
  if (empty_rows > 0) stream << "Removed " << empty_rows << " empty rows\n";
}

void PreprocessFixedVars::apply(Model& model) {
  Int& n = model.n_;
  Int& m = model.m_;
  HighsSparseMatrix& A = model.A_;
  std::vector<double>& b = model.b_;
  std::vector<double>& c = model.c_;
  std::vector<double>& lower = model.lower_;
  std::vector<double>& upper = model.upper_;
  HighsHessian& Q = model.Q_;
  double& offset = model.offset_;

  n_pre = n;
  m_pre = m;

  // See "Preprocessing for quadratic programming", Gould, Toint, Math Program
  fixed_vars = 0;
  for (Int j = 0; j < n; ++j) {
    if (lower[j] == upper[j]) ++fixed_vars;
    data[j].c = c[j];
  }

  // cannot remove all variables
  if (fixed_vars == n) fixed_vars = 0;

  if (fixed_vars > 0) {
    fixed_at.assign(n, kHighsInf);
    std::vector<Int> index_to_remove{};
    for (Int j = 0; j < n; ++j) {
      if (lower[j] == upper[j]) {
        fixed_at[j] = lower[j];
        index_to_remove.push_back(j);
        const double xcol = fixed_at[j];

        FixedVarsData& dataj = data[j];

        offset += dataj.c * xcol;
        if (model.qp()) offset += 0.5 * Q.diag(j) * xcol * xcol;

        for (Int el = A.start_[j]; el < A.start_[j + 1]; ++el) {
          const Int row = A.index_[el];
          const double val = A.value_[el];
          b[row] -= val * xcol;

          dataj.indA.push_back(row);
          dataj.valA.push_back(val);
        }

        if (model.qp()) {
          for (Int colQ = 0; colQ < j; ++colQ) {
            for (Int el = Q.start_[colQ]; el < Q.start_[colQ + 1]; ++el) {
              const Int rowQ = Q.index_[el];
              if (rowQ == j) {
                c[colQ] += Q.value_[el] * xcol;

                dataj.indQ.push_back(colQ);
                dataj.valQ.push_back(Q.value_[el]);
              }
            }
          }
          for (Int el = Q.start_[j]; el < Q.start_[j + 1]; ++el) {
            const Int rowQ = Q.index_[el];
            c[rowQ] += Q.value_[el] * xcol;

            dataj.indQ.push_back(rowQ);
            dataj.valQ.push_back(Q.value_[el]);
          }
        }
      }
    }

    HighsIndexCollection index_collection;
    create(index_collection, index_to_remove.size(), index_to_remove.data(), n);
    A.deleteCols(index_collection);
    if (model.qp()) Q.deleteCols(index_collection);

    Int next = 0;
    Int copy_to = 0;
    for (Int i = 0; i < n; ++i) {
      if (next < static_cast<Int>(index_to_remove.size()) &&
          i == index_to_remove[next]) {
        ++next;
        continue;
      } else {
        c[copy_to] = c[i];
        lower[copy_to] = lower[i];
        upper[copy_to] = upper[i];
        copy_to++;
      }
    }

    n -= fixed_vars;
    assert(A.num_col_ == n);
    assert(!model.qp() || Q.dim_ == n);
    c.resize(n);
    lower.resize(n);
    upper.resize(n);
  }

  n_post = n;
  m_post = m;
}

void PreprocessFixedVars::undo(PreprocessorPoint& point, const Model& model,
                               const Iterate& it) const {
  point.assertConsistency(n_post, m_post);

  if (fixed_vars > 0) {
    // Add primal and dual variables for fixed variables

    std::vector<double> new_x(n_pre, 0.0);
    std::vector<double> new_xl(n_pre, 0.0);
    std::vector<double> new_xu(n_pre, 0.0);
    std::vector<double> new_zl(n_pre, 0.0);
    std::vector<double> new_zu(n_pre, 0.0);

    Int pos{};
    for (Int j = 0; j < n_pre; ++j) {
      if (std::isfinite(fixed_at[j])) {
        new_x[j] = fixed_at[j];
        new_xl[j] = 0.0;
        new_xu[j] = 0.0;
      } else {
        new_x[j] = point.x[pos];
        new_xl[j] = point.xl[pos];
        new_xu[j] = point.xu[pos];
        new_zl[j] = point.zl[pos];
        new_zu[j] = point.zu[pos];
        ++pos;
      }
    }

    for (Int j = 0; j < n_pre; ++j) {
      if (std::isfinite(fixed_at[j])) {
        // compute dual variables so that they are dual feasible
        // z = c - A^T * y + Q * x
        // need to do this after all x have been computed, due to Q*x term
        const auto& dataj = data.at(j);
        double z = dataj.c;
        for (Int i = 0; i < static_cast<Int>(dataj.indA.size()); ++i) {
          const Int row = dataj.indA[i];
          const double val = dataj.valA[i];
          z -= val * point.y[row];
        }
        if (model.qp()) {
          for (Int i = 0; i < static_cast<Int>(dataj.indQ.size()); ++i) {
            const Int row = dataj.indQ[i];
            const double val = dataj.valQ[i];
            z += val * new_x[row];
          }
        }

        if (z >= 0.0) {
          new_zl[j] = z;
          new_zu[j] = 0.0;
        } else {
          new_zl[j] = 0.0;
          new_zu[j] = -z;
        }
      }
    }

    point.x = std::move(new_x);
    point.xl = std::move(new_xl);
    point.xu = std::move(new_xu);
    point.zl = std::move(new_zl);
    point.zu = std::move(new_zu);
  }

  point.assertConsistency(n_pre, m_pre);
}

void PreprocessFixedVars::print(std::stringstream& stream) const {
  if (fixed_vars > 0)
    stream << "Removed " << fixed_vars << " fixed variables\n";
}

static double roundToPowerOf2(double d) {
  int exp;  // int, not Int

  // d = x * 2^exp, d\in[0.5,1)
  double x = std::frexp(d, &exp);

  // d is already exact power of 2
  if (x == 0.5) return d;

  // return 1*2^exp
  return std::ldexp(1.0, exp);
}

void PreprocessScaling::apply(Model& model) {
  Int& n = model.n_;
  Int& m = model.m_;
  HighsSparseMatrix& A = model.A_;
  std::vector<double>& b = model.b_;
  std::vector<double>& c = model.c_;
  std::vector<double>& lower = model.lower_;
  std::vector<double>& upper = model.upper_;
  HighsHessian& Q = model.Q_;

  std::vector<double>& colscale = model.colscale_;
  std::vector<double>& rowscale = model.rowscale_;

  n_pre = n;
  m_pre = m;
  n_post = n;
  m_post = m;

  // *********************************************************************
  // Compute scaling
  // *********************************************************************
  // Transformation:
  // A -> R * A * C
  // b -> R * b
  // c -> C * c
  // x -> C^-1 * x
  // y -> R^-1 * y
  // z -> C * z
  // Q -> C * Q * C
  // where R is row scaling, C is col scaling.

  colscale.assign(n, 1.0);
  rowscale.assign(m, 1.0);

  auto colScaling = [&]() {
    // infinity norm of columns
    std::vector<double> norm_cols(n);
    for (Int col = 0; col < n; ++col) {
      for (Int el = A.start_[col]; el < A.start_[col + 1]; ++el) {
        const Int row = A.index_[el];
        double value = A.value_[el];
        value *= colscale[col];
        value *= rowscale[row];
        norm_cols[col] = std::max(norm_cols[col], std::abs(value));
      }
    }

    // apply col scaling
    for (Int i = 0; i < n; ++i) {
      double coeff = 1.0 / std::sqrt(norm_cols[i]);

      // use geomean of norm scaling and scaling that makes diagonal of Q equal
      // to 1
      if (model.qp()) {
        const double Qii = Q.diag(i) * colscale[i] * colscale[i];
        coeff = std::sqrt(coeff * 1.0 / std::sqrt(std::abs(Qii)));
      }

      // if bounds interval is too small, use geomean of norm scaling
      // and maximum scaling allowed by bounds
      if (std::isfinite(lower[i]) && std::isfinite(upper[i])) {
        const double l = lower[i] / colscale[i];
        const double u = upper[i] / colscale[i];
        const double diff = std::abs(u - l);

        if (diff / coeff < kScalingSmallBoundDiff)
          coeff = std::sqrt(coeff * diff / kScalingSmallBoundDiff);
      }

      if (!std::isinf(coeff) && !std::isnan(coeff)) colscale[i] *= coeff;

      colscale[i] = std::max(colscale[i], kScalingSmallCoeff);
      colscale[i] = std::min(colscale[i], kScalingLargeCoeff);
    }
  };
  auto rowScaling = [&]() {
    // infinity norm of rows
    std::vector<double> norm_rows(m);
    for (Int col = 0; col < n; ++col) {
      for (Int el = A.start_[col]; el < A.start_[col + 1]; ++el) {
        const Int row = A.index_[el];
        double value = A.value_[el];
        value *= colscale[col];
        value *= rowscale[row];
        norm_rows[row] = std::max(norm_rows[row], std::abs(value));
      }
    }

    // apply row scaling
    for (Int i = 0; i < m; ++i) {
      if (norm_rows[i] > 0.0) rowscale[i] *= 1.0 / std::sqrt(norm_rows[i]);
      rowscale[i] = std::max(rowscale[i], kScalingSmallCoeff);
      rowscale[i] = std::min(rowscale[i], kScalingLargeCoeff);
    }
  };

  const Int num_passes = 10;
  for (Int pass = 0; pass < num_passes; ++pass) {
    rowScaling();
    colScaling();
  }

  bool scaling_failed = isInfVector(colscale) || isNanVector(colscale) ||
                        isInfVector(rowscale) || isNanVector(rowscale);
  if (scaling_failed) {
    colscale.clear();
    rowscale.clear();
    return;
  }

  scaled = true;

  // *********************************************************************
  // Adjust scaling
  // *********************************************************************

  for (Int i = 0; i < n; ++i) {
    colscale[i] = std::max(colscale[i], kScalingSmallCoeff);
    colscale[i] = std::min(colscale[i], kScalingLargeCoeff);
    colscale[i] = roundToPowerOf2(colscale[i]);
  }
  for (Int i = 0; i < m; ++i) {
    rowscale[i] = std::max(rowscale[i], kScalingSmallCoeff);
    rowscale[i] = std::min(rowscale[i], kScalingLargeCoeff);
    rowscale[i] = roundToPowerOf2(rowscale[i]);
  }

  // *********************************************************************
  // Apply scaling
  // *********************************************************************

  // Column has been scaled up by colscale_[col], so cost is scaled up and
  // bounds are scaled down
  for (Int col = 0; col < n; ++col) {
    c[col] *= colscale[col];
    lower[col] /= colscale[col];
    upper[col] /= colscale[col];
  }

  // Row has been scaled up by rowscale_[row], so b is scaled up
  for (Int row = 0; row < m; ++row) b[row] *= rowscale[row];

  // Each entry of the matrix is scaled by the corresponding row and col
  // factor
  for (Int col = 0; col < n; ++col) {
    for (Int el = A.start_[col]; el < A.start_[col + 1]; ++el) {
      Int row = A.index_[el];
      A.value_[el] *= rowscale[row];
      A.value_[el] *= colscale[col];
    }
  }

  if (model.qp()) {
    for (Int col = 0; col < Q.dim_; ++col) {
      for (Int el = Q.start_[col]; el < Q.start_[col + 1]; ++el) {
        Int row = Q.index_[el];
        Q.value_[el] *= colscale[row];
        Q.value_[el] *= colscale[col];
      }
    }
  }
}

void PreprocessScaling::undo(PreprocessorPoint& point, const Model& model,
                             const Iterate& it) const {
  point.assertConsistency(n_post, m_post);
  if (model.scaled()) {
    const std::vector<double>& colscale = model.colscale_;
    const std::vector<double>& rowscale = model.rowscale_;

    for (Int i = 0; i < n_post; ++i) {
      point.x[i] *= colscale[i];
      point.xl[i] *= colscale[i];
      point.xu[i] *= colscale[i];
      point.zl[i] /= colscale[i];
      point.zu[i] /= colscale[i];
    }
    for (Int i = 0; i < m_post; ++i) {
      point.y[i] *= rowscale[i];
      point.slack[i] /= rowscale[i];
    }

    // set variables that were ignored
    for (Int i = 0; i < n_post; ++i) {
      if (!model.hasLb(i)) {
        point.xl[i] = kHighsInf;
        point.zl[i] = 0.0;
      }
      if (!model.hasUb(i)) {
        point.xu[i] = kHighsInf;
        point.zu[i] = 0.0;
      }
    }
    for (Int i : model.freeVars()) {
      point.xl[i] = kHighsInf;
      point.zl[i] = 0.0;
      point.xu[i] = kHighsInf;
      point.zu[i] = 0.0;
    }
  }
  point.assertConsistency(n_pre, m_pre);
}

void PreprocessScaling::print(std::stringstream& stream) const {}

void PreprocessFormulation::apply(Model& model) {
  Int& n = model.n_;
  Int& m = model.m_;
  HighsSparseMatrix& A = model.A_;
  std::vector<double>& b = model.b_;
  std::vector<double>& c = model.c_;
  std::vector<double>& lower = model.lower_;
  std::vector<double>& upper = model.upper_;
  std::vector<char>& constraints = model.constraints_;
  HighsHessian& Q = model.Q_;

  n_pre = n;
  m_pre = m;

  Int Annz = A.numNz();

  for (Int i = 0; i < m; ++i) {
    if (constraints[i] != '=') {
      // inequality constraint, add slack variable

      ++n;

      // lower/upper bound for new slack
      if (constraints[i] == '>') {
        lower.push_back(-kHighsInf);
        upper.push_back(0.0);
      } else {
        lower.push_back(0.0);
        upper.push_back(kHighsInf);
      }

      // cost for new slack
      c.push_back(0.0);

      // add column of identity to A_
      std::vector<Int> temp_ind{i};
      std::vector<double> temp_val{1.0};
      A.addVec(1, temp_ind.data(), temp_val.data());

      // set scaling to 1
      if (model.scaled()) model.colscale_.push_back(1.0);
    }
  }

  if (model.qp()) completeHessian(n, Q);

  n_post = n;
  m_post = m;
}

void PreprocessFormulation::undo(PreprocessorPoint& point, const Model& model,
                                 const Iterate& it) const {
  it.assertConsistency(n_post, m_post);

  // Copy x, xl, xu, zl, zu without slacks
  point.x = std::vector<double>(it.x.begin(), it.x.begin() + n_pre);
  point.xl = std::vector<double>(it.xl.begin(), it.xl.begin() + n_pre);
  point.xu = std::vector<double>(it.xu.begin(), it.xu.begin() + n_pre);
  point.zl = std::vector<double>(it.zl.begin(), it.zl.begin() + n_pre);
  point.zu = std::vector<double>(it.zu.begin(), it.zu.begin() + n_pre);

  // force unused entries to have correct value
  for (Int i = 0; i < n_pre; ++i) {
    if (!model.hasLb(i)) {
      point.xl[i] = kHighsInf;
      point.zl[i] = 0.0;
    }
    if (!model.hasUb(i)) {
      point.xu[i] = kHighsInf;
      point.zu[i] = 0.0;
    }
  }
  for (Int i : model.freeVars()) {
    point.xl[i] = kHighsInf;
    point.zl[i] = 0.0;
    point.xu[i] = kHighsInf;
    point.zu[i] = 0.0;
  }

  // For the Lagrange multipliers, use slacks from zl and zu, to get correct
  // sign. NB: there is no explicit slack stored for equality constraints.
  point.y.resize(m_pre);
  Int slack_pos = 0;
  for (Int i = 0; i < m_pre; ++i) {
    switch (model.constraint(i)) {
      case '=':
        point.y[i] = it.y[i];
        break;
      case '>':
        point.y[i] = it.zu[n_pre + slack_pos];
        ++slack_pos;
        break;
      case '<':
        point.y[i] = -it.zl[n_pre + slack_pos];
        ++slack_pos;
        break;
    }
  }

  // For x-slacks, use slacks from xl and xu, to get correct sign.
  // NB: there is no explicit slack stored for equality constraints.
  point.slack.resize(m_pre);
  slack_pos = 0;
  for (Int i = 0; i < m_pre; ++i) {
    switch (model.constraint(i)) {
      case '=':
        point.slack[i] = 0.0;
        break;
      case '>':
        point.slack[i] = -it.xu[n_pre + slack_pos];
        ++slack_pos;
        break;
      case '<':
        point.slack[i] = it.xl[n_pre + slack_pos];
        ++slack_pos;
        break;
    }
  }
  point.assertConsistency(n_pre, m_pre);
}

void PreprocessFormulation::print(std::stringstream& stream) const {
  stream << "Added " << n_post - n_pre << " slacks\n";
}

void PreprocessFreeVars::apply(Model& model) {
  Int& n = model.n_;
  std::vector<double>& lower = model.lower_;
  std::vector<double>& upper = model.upper_;
  std::vector<Int>& free_variables = model.free_variables_;

  free_variables.clear();
  for (Int i = 0; i < n; ++i) {
    if (!std::isfinite(lower[i]) && !std::isfinite(upper[i]) &&
        lower[i] != upper[i]) {
      // free variable
      free_variables.push_back(i);
      ++free_vars_count;
      lower[i] = -kFreeVarsInitialBound;
      upper[i] = kFreeVarsInitialBound;
    }
  }
}

void PreprocessFreeVars::undo(PreprocessorPoint& point, const Model& model,
                              const Iterate& it) const {}

void PreprocessFreeVars::print(std::stringstream& stream) const {
  if (free_vars_count > 0)
    stream << "Found " << free_vars_count << " free variables\n";
}

#define APPLY_ACTION(T)                                      \
  stack.push_back(std::unique_ptr<PreprocessAction>(new T)); \
  stack.back()->apply(model);

void Preprocessor::apply(Model& model) {
  // Remove fixed variables before removing empty rows, because removing columns
  // may create empty rows

  APPLY_ACTION(PreprocessFixedVars);
  APPLY_ACTION(PreprocessEmptyRows);
  APPLY_ACTION(PreprocessScaling);
  APPLY_ACTION(PreprocessFormulation);
  APPLY_ACTION(PreprocessFreeVars);
}
void Preprocessor::undo(PreprocessorPoint& point, const Model& model,
                        const Iterate& it) const {
  for (auto iterator = stack.rbegin(); iterator != stack.rend(); ++iterator) {
    (*iterator)->undo(point, model, it);
  }
}

void Preprocessor::print(std::stringstream& log_stream) const {
  for (auto iterator = stack.begin(); iterator != stack.end(); ++iterator) {
    (*iterator)->print(log_stream);
  }
}

}  // namespace hipo

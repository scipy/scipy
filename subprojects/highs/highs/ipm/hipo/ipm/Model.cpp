#include "Model.h"

#include "Parameters.h"
#include "Status.h"
#include "ipm/IpxWrapper.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "model/HighsHessianUtils.h"

namespace hipo {

Int Model::init(const HighsLp& lp, const HighsHessian& Q) {
  fillInIpxData(lp, n_, m_, offset_, c_, lower_, upper_, A_.start_, A_.index_,
                A_.value_, b_, constraints_);
  rhs_orig_ = b_;
  constraints_orig_ = constraints_;
  Q_ = Q;
  if (qp()) completeHessian(n_, Q_);
  sense_ = lp.sense_;

  if (checkData()) return kErrorModel;

  lp_orig_ = &lp;
  n_orig_ = n_;
  m_orig_ = m_;
  A_.num_col_ = n_;
  A_.num_row_ = m_;

  preprocess();

  nzBounds();
  computeNorms();

  // double transpose to sort indices of each column
  A_.ensureRowwise();
  A_.ensureColwise();

  if (checkData()) return kErrorModel;

  ready_ = true;

  return 0;
}

void Model::nzBounds() {
  // compute lower and upper bounds for the number of nonzeros in normal
  // equations.
  std::vector<HighsBool> mark(m_, false);
  NE_nz_lb_ = A_.num_row_;
  NE_nz_ub_ = A_.num_row_;
  for (Int col = 0; col < A_.num_col_; ++col) {
    Int64 used = 0;
    Int64 unused = 0;
    for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      const Int row = A_.index_[el];
      if (mark[row])
        ++used;
      else {
        mark[row] = true;
        ++unused;
      }
    }

    NE_nz_ub_ += (used + unused) * (used + unused - 1) / 2;
    NE_nz_lb_ += unused * (unused - 1) / 2 + used * unused;
  }
  NE_nz_ub_ = std::min(NE_nz_ub_, (Int64)A_.num_row_ * (A_.num_row_ + 1) / 2);

  AS_nz_ = A_.numNz() + A_.num_row_ + (qp() ? Q_.numNz() : A_.num_col_);
}

Int Model::checkData() const {
  // Check if model provided by the user is ok.
  // Return kErrorModel if something is wrong.

  // Dimensions are valid
  if (n_ <= 0 || m_ < 0) return kErrorModel;

  // Vectors are of correct size
  if (static_cast<Int>(c_.size()) != n_ || static_cast<Int>(b_.size()) != m_ ||
      static_cast<Int>(lower_.size()) != n_ ||
      static_cast<Int>(upper_.size()) != n_ ||
      static_cast<Int>(constraints_.size()) != m_ ||
      static_cast<Int>(A_.start_.size()) != n_ + 1 ||
      static_cast<Int>(A_.index_.size()) != A_.start_.back() ||
      static_cast<Int>(A_.value_.size()) != A_.start_.back())
    return kErrorModel;

  // Hessian is ok, for QPs only
  if (qp() && (Q_.dim_ != n_ || Q_.format_ != HessianFormat::kTriangular))
    return kErrorModel;

  // Vectors are valid
  for (Int i = 0; i < n_; ++i)
    if (!std::isfinite(c_[i])) return kErrorModel;
  for (Int i = 0; i < m_; ++i)
    if (!std::isfinite(b_[i])) return kErrorModel;
  for (Int i = 0; i < n_; ++i) {
    if (!std::isfinite(lower_[i]) && lower_[i] != -INFINITY) return kErrorModel;
    if (!std::isfinite(upper_[i]) && upper_[i] != INFINITY) return kErrorModel;
    if (lower_[i] > upper_[i]) return kErrorModel;
  }
  for (Int i = 0; i < m_; ++i)
    if (constraints_[i] != '<' && constraints_[i] != '=' &&
        constraints_[i] != '>')
      return kErrorModel;

  // Matrix is valid
  for (Int i = 0; i < A_.start_[n_]; ++i)
    if (!std::isfinite(A_.value_[i])) return kErrorModel;

  return 0;
}

void Model::preprocess() { preprocessor_.apply(*this); }

void Model::postprocess(std::vector<double>& x, std::vector<double>& xl,
                        std::vector<double>& xu, std::vector<double>& slack,
                        std::vector<double>& y, std::vector<double>& zl,
                        std::vector<double>& zu, const Iterate& it) const {
  PreprocessorPoint point{x, xl, xu, slack, y, zl, zu};
  preprocessor_.undo(point, *this, it);
}

void Model::computeNorms() {
  norm_scaled_obj_ = infNorm(c_);

  norm_unscaled_obj_ = 0.0;
  for (Int i = 0; i < n_; ++i) {
    double val = std::abs(c_[i]);
    if (scaled()) val /= colScale(i);
    norm_unscaled_obj_ = std::max(norm_unscaled_obj_, val);
  }

  norm_scaled_rhs_ = infNorm(b_);
  for (double d : lower_)
    if (std::isfinite(d))
      norm_scaled_rhs_ = std::max(norm_scaled_rhs_, std::abs(d));
  for (double d : upper_)
    if (std::isfinite(d))
      norm_scaled_rhs_ = std::max(norm_scaled_rhs_, std::abs(d));

  norm_unscaled_rhs_ = 0.0;
  for (Int i = 0; i < m_; ++i) {
    double val = std::abs(b_[i]);
    if (scaled()) val /= rowScale(i);
    norm_unscaled_rhs_ = std::max(norm_unscaled_rhs_, val);
  }
  for (Int i = 0; i < n_; ++i) {
    if (std::isfinite(lower_[i])) {
      double val = std::abs(lower_[i]);
      if (scaled()) val *= colScale(i);
      norm_unscaled_rhs_ = std::max(norm_unscaled_rhs_, val);
    }
    if (std::isfinite(upper_[i])) {
      double val = std::abs(upper_[i]);
      if (scaled()) val *= colScale(i);
      norm_unscaled_rhs_ = std::max(norm_unscaled_rhs_, val);
    }
  }

  // norms of rows and cols of A
  one_norm_cols_.resize(n_);
  one_norm_rows_.resize(m_);
  inf_norm_cols_.resize(n_);
  inf_norm_rows_.resize(m_);
  for (Int col = 0; col < n_; ++col) {
    for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      Int row = A_.index_[el];
      double val = A_.value_[el];
      one_norm_cols_[col] += std::abs(val);
      one_norm_rows_[row] += std::abs(val);
      inf_norm_rows_[row] = std::max(inf_norm_rows_[row], std::abs(val));
      inf_norm_cols_[col] = std::max(inf_norm_cols_[col], std::abs(val));
    }
  }
}

void Model::print(const Logger& logger) const {
  std::stringstream log_stream;

  log_stream << textline("Rows:") << sci(m_, 0, 1) << '\n';
  log_stream << textline("Cols:") << sci(n_, 0, 1) << '\n';
  log_stream << textline("Nnz A:") << sci(A_.numNz(), 0, 1) << '\n';

  if (qp()) {
    log_stream << textline("Nnz Q:") << sci(Q_.numNz(), 0, 1);
    if (nonSeparableQp())
      log_stream << ", non-separable\n";
    else
      log_stream << ", separable\n";
  }

  // compute max and min entry of A in absolute value
  double Amin = kHighsInf;
  double Amax = 0.0;
  std::vector<double> norm_inf_cols(n_);
  std::vector<double> norm_inf_rows(m_);
  std::vector<double> norm_one_cols(n_);
  std::vector<double> norm_one_rows(m_);
  for (Int col = 0; col < n_; ++col) {
    for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      const Int row = A_.index_[el];
      const double value = A_.value_[el];
      if (value != 0.0) {
        Amin = std::min(Amin, std::abs(value));
        Amax = std::max(Amax, std::abs(value));
        norm_inf_cols[col] = std::max(norm_inf_cols[col], std::abs(value));
        norm_inf_rows[row] = std::max(norm_inf_rows[row], std::abs(value));
        norm_one_cols[col] += std::abs(value);
        norm_one_rows[row] += std::abs(value);
      }
    }
  }
  if (std::isinf(Amin)) Amin = 0.0;

  double norm_inf_col_min = kHighsInf;
  double norm_inf_col_max = 0.0;
  double norm_inf_row_min = kHighsInf;
  double norm_inf_row_max = 0.0;
  for (double d : norm_inf_cols) {
    norm_inf_col_min = std::min(norm_inf_col_min, d);
    norm_inf_col_max = std::max(norm_inf_col_max, d);
  }
  for (double d : norm_inf_rows) {
    norm_inf_row_min = std::min(norm_inf_row_min, d);
    norm_inf_row_max = std::max(norm_inf_row_max, d);
  }

  double norm_one_col_min = kHighsInf;
  double norm_one_col_max = 0.0;
  double norm_one_row_min = kHighsInf;
  double norm_one_row_max = 0.0;
  for (double d : norm_one_cols) {
    norm_one_col_min = std::min(norm_one_col_min, d);
    norm_one_col_max = std::max(norm_one_col_max, d);
  }
  for (double d : norm_one_rows) {
    norm_one_row_min = std::min(norm_one_row_min, d);
    norm_one_row_max = std::max(norm_one_row_max, d);
  }

  // compute max and min entry of c
  double cmin = kHighsInf;
  double cmax = 0.0;
  for (Int i = 0; i < n_; ++i) {
    if (c_[i] != 0.0) {
      cmin = std::min(cmin, std::abs(c_[i]));
      cmax = std::max(cmax, std::abs(c_[i]));
    }
  }
  if (std::isinf(cmin)) cmin = 0.0;

  // compute max and min entry of b
  double bmin = kHighsInf;
  double bmax = 0.0;
  for (Int i = 0; i < m_; ++i) {
    if (b_[i] != 0.0) {
      bmin = std::min(bmin, std::abs(b_[i]));
      bmax = std::max(bmax, std::abs(b_[i]));
    }
  }
  if (std::isinf(bmin)) bmin = 0.0;

  // compute max and min entry of Q in absolute value
  double Qmin = kHighsInf;
  double Qmax = 0.0;
  for (double val : Q_.value_) {
    if (val != 0.0) {
      Qmin = std::min(Qmin, std::abs(val));
      Qmax = std::max(Qmax, std::abs(val));
    }
  }
  if (std::isinf(Qmin)) Qmin = 0.0;

  // compute max and min for bounds intervals
  double boundmin = kHighsInf;
  double boundmax = 0.0;
  for (Int i = 0; i < n_; ++i) {
    if (std::isfinite(lower_[i]) && std::isfinite(upper_[i])) {
      const double diff = std::abs(upper_[i] - lower_[i]);
      boundmin = std::min(boundmin, diff);
      boundmax = std::max(boundmax, diff);
    }
  }
  if (std::isinf(boundmin)) boundmin = 0.0;

  // compute max and min scaling
  double scalemin = kHighsInf;
  double scalemax = 0.0;
  if (scaled()) {
    for (Int i = 0; i < n_; ++i) {
      scalemin = std::min(scalemin, colScale(i));
      scalemax = std::max(scalemax, colScale(i));
    }
    for (Int i = 0; i < m_; ++i) {
      scalemin = std::min(scalemin, rowScale(i));
      scalemax = std::max(scalemax, rowScale(i));
    }
  }
  if (std::isinf(scalemin)) scalemin = 0.0;

  log_stream << textline("Range of A:") << "[" << sci(Amin, 5, 1) << ", "
             << sci(Amax, 5, 1) << "]\n";

  if (logger.debug(1)) {
    log_stream << textline("Inf-norm rows:") << "["
               << sci(norm_inf_row_min, 5, 1) << ", "
               << sci(norm_inf_row_max, 5, 1) << "]\n";

    log_stream << textline("Inf-norm cols:") << "["
               << sci(norm_inf_col_min, 5, 1) << ", "
               << sci(norm_inf_col_max, 5, 1) << "]\n";

    log_stream << textline("One-norm rows:") << "["
               << sci(norm_one_row_min, 5, 1) << ", "
               << sci(norm_one_row_max, 5, 1) << "]\n";

    log_stream << textline("One-norm cols:") << "["
               << sci(norm_one_col_min, 5, 1) << ", "
               << sci(norm_one_col_max, 5, 1) << "]\n";
  }

  log_stream << textline("Range of b:") << "[" << sci(bmin, 5, 1) << ", "
             << sci(bmax, 5, 1) << "]\n";

  log_stream << textline("Range of c:") << "[" << sci(cmin, 5, 1) << ", "
             << sci(cmax, 5, 1) << "]\n";

  if (qp()) {
    log_stream << textline("Range of Q:") << "[" << sci(Qmin, 5, 1) << ", "
               << sci(Qmax, 5, 1) << "]\n";
  }

  log_stream << textline("Range of bound intervals:") << "["
             << sci(boundmin, 5, 1) << ", " << sci(boundmax, 5, 1) << "]\n";

  log_stream << textline("Scaling coefficients:") << "[" << sci(scalemin, 5, 1)
             << ", " << sci(scalemax, 5, 1) << "]\n";

  if (logger.debug(1)) {
    preprocessor_.print(log_stream);

    log_stream << "Expected nnz: ";
    log_stream << "AS " << sci(AS_nz_, 0, 1) << "; ";
    log_stream << "NE [" << sci(NE_nz_lb_, 0, 1) << ", " << sci(NE_nz_ub_, 0, 1)
               << "]\n";
  }

  logger.print(log_stream.str().c_str());
}

Int Model::loadIntoIpx(ipx::LpSolver& lps) const {
  Int ipx_m, ipx_n;
  std::vector<double> ipx_b, ipx_c, ipx_lower, ipx_upper, ipx_A_vals;
  std::vector<Int> ipx_A_ptr, ipx_A_rows;
  std::vector<char> ipx_constraints;
  double ipx_offset;

  if (!lp_orig_) return kErrorInvalidPointer;

  fillInIpxData(*lp_orig_, ipx_n, ipx_m, ipx_offset, ipx_c, ipx_lower,
                ipx_upper, ipx_A_ptr, ipx_A_rows, ipx_A_vals, ipx_b,
                ipx_constraints);

  Int load_status = lps.LoadModel(
      ipx_n, ipx_offset, ipx_c.data(), ipx_lower.data(), ipx_upper.data(),
      ipx_m, ipx_A_ptr.data(), ipx_A_rows.data(), ipx_A_vals.data(),
      ipx_b.data(), ipx_constraints.data());

  return load_status;
}

void Model::printDense() const {
  std::vector<std::vector<double>> Adense(m_, std::vector<double>(n_, 0.0));
  std::vector<std::vector<double>> Qdense(n_, std::vector<double>(n_, 0.0));

  for (Int col = 0; col < n_; ++col) {
    for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      const Int row = A_.index_[el];
      const double val = A_.value_[el];
      Adense[row][col] = val;
    }
  }
  for (Int col = 0; col < n_; ++col) {
    for (Int el = Q_.start_[col]; el < Q_.start_[col + 1]; ++el) {
      const Int row = Q_.index_[el];
      const double val = Q_.value_[el];
      Qdense[row][col] = val;
    }
  }

  printf("\nA\n");
  for (Int i = 0; i < m_; ++i) {
    for (Int j = 0; j < n_; ++j) printf("%6.2f ", Adense[i][j]);
    printf("\n");
  }
  printf("b\n");
  for (Int i = 0; i < m_; ++i) printf("%6.2f ", b_[i]);
  printf("\n");
  printf("c\n");
  for (Int i = 0; i < n_; ++i) printf("%6.2f ", c_[i]);
  printf("\n");
  printf("lb\n");
  for (Int i = 0; i < n_; ++i) printf("%6.2f ", lower_[i]);
  printf("\n");
  printf("ub\n");
  for (Int i = 0; i < n_; ++i) printf("%6.2f ", upper_[i]);
  printf("\n");
  printf("Q\n");
  for (Int i = 0; i < n_; ++i) {
    for (Int j = 0; j < n_; ++j) printf("%6.2f ", Qdense[i][j]);
    printf("\n");
  }
  printf("offset %6.2f\n", offset_);
}

void Model::adjustFreeVars(std::vector<double>& x, std::vector<double>& xl,
                           std::vector<double>& xu, const Logger& logger) {
  for (Int i : free_variables_) {
    if (x[i] < lower_[i] * kFreeVarsCloseRatio) {
      // getting close to lower bound
      const double new_lower = x[i] / kFreeVarsCloseRatio;
      logger.printDetailed(
          "Free var %5d is at %8.1e with lb %8.1e, lb changed to %8.1e\n", i,
          x[i], lower_[i], new_lower);
      lower_[i] = new_lower;
      xl[i] = x[i] - lower_[i];
    }
    if (x[i] > upper_[i] * kFreeVarsCloseRatio) {
      // getting close to upper bound
      const double new_upper = x[i] / kFreeVarsCloseRatio;
      logger.printDetailed(
          "Free var %5d is at %8.1e with ub %8.1e, ub changed to %8.1e\n", i,
          x[i], upper_[i], new_upper);
      upper_[i] = new_upper;
      xu[i] = upper_[i] - x[i];
    }
  }
}

void Model::saveBounds() {
  best_lower_.resize(free_variables_.size());
  best_upper_.resize(free_variables_.size());
  Int next = 0;
  for (Int i : free_variables_) {
    best_lower_[next] = lower_[i];
    best_upper_[next] = upper_[i];
    ++next;
  }
}
void Model::resetBounds() {
  Int next = 0;
  for (Int i : free_variables_) {
    lower_[i] = best_lower_[next];
    upper_[i] = best_upper_[next];
    ++next;
  }
}

}  // namespace hipo
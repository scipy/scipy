#include "KrylovMethodsIpm.h"

namespace hipo {

void IpmMatrix::reset(const HighsSparseMatrix& A,
                      const std::vector<double>& scaling, bool use_as) {
  A_ = &A;
  scaling_ = &scaling;
  use_as_ = use_as;
}

void IpmMatrix::apply(std::vector<double>& x) const {
  const Int n = A_->num_col_;
  const Int m = A_->num_row_;

  if (use_as_) {
    // compute [y_1; y_2] = [-Theta^-1 * x_1 + A^T * x_2; A * x_1]

    // split x
    std::vector<double> x_1(x.begin(), x.begin() + n);
    std::vector<double> x_2(x.begin() + n, x.end());

    std::vector<double> y_1(n);
    std::vector<double> y_2(m);

    // y_1 = A^T * x_2
    A_->alphaProductPlusY(1.0, x_2, y_1, true);

    // y_1 -= Theta^-1 * x_1
    for (Int i = 0; i < n; ++i) y_1[i] -= (*scaling_)[i] * x_1[i];

    // y_2 = A * x_1
    A_->alphaProductPlusY(1.0, x_1, y_2);

    // put result back into x
    x.clear();
    x.insert(x.begin(), y_1.begin(), y_1.end());
    x.insert(x.begin() + n, y_2.begin(), y_2.end());

  } else {
    // compute A * Theta * A^T * x

    std::vector<double> w(n);

    // w = A^T * x
    A_->alphaProductPlusY(1.0, x, w, true);

    // w = Theta * w
    for (Int i = 0; i < n; ++i) w[i] /= (*scaling_)[i];

    // x = A * w
    A_->alphaProductPlusY(1.0, w, x);
  }
}

IpmFactor::IpmFactor(const Numeric& N) : N_{N} {}
void IpmFactor::apply(std::vector<double>& x) const { N_.solve(x.data()); }

void NeDiagPrec::reset(const HighsSparseMatrix& A,
                       const std::vector<double>& scaling) {
  // prepare diagonal preconditioner
  diag.assign(A.num_row_, 0.0);

  // build diagonal of normal equations
  for (Int c = 0; c < A.num_col_; ++c) {
    for (Int el = A.start_[c]; el < A.start_[c + 1]; ++el) {
      Int r = A.index_[el];
      double v = A.value_[el];

      diag[r] += v * v / scaling[c];
    }
  }

  // compute inverse of diagonal entries
  for (Int i = 0; i < static_cast<Int>(diag.size()); ++i)
    diag[i] = 1.0 / diag[i];
}
void NeDiagPrec::apply(std::vector<double>& x) const {
  // apply diagonal preconditioner
  for (Int i = 0; i < static_cast<Int>(diag.size()); ++i) x[i] *= diag[i];
}

}  // namespace hipo
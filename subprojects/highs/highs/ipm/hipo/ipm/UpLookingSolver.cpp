#include "UpLookingSolver.h"

#include "HighsExternalApi.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

UpLookingSolver::UpLookingSolver(KktMatrix& kkt, Info& info, IpmData& data,
                                 const Regularisation& regul,
                                 const Model& model, const Control& control)
    : kkt_{kkt},
      ptr_{kkt_.isNE() ? kkt_.ptrNE : kkt_.ptrAS},
      rows_{kkt_.isNE() ? kkt_.rowsNE : kkt_.rowsAS},
      val_{kkt_.isNE() ? kkt_.valNE : kkt_.valAS},
      n_{static_cast<Int>(ptr_.size() - 1)},
      info_{info},
      data_{data},
      regul_{regul},
      model_{model},
      control_{control} {}

void UpLookingSolver::etreeAndCounts(const std::vector<Int>& ptr,
                                     const std::vector<Int>& rows,
                                     std::vector<Int64>& colcount) {
  // compute elimination tree and column counts
  parent_.resize(n_);
  colcount.resize(n_);
  std::vector<Int> work(n_);
  HighsExtras::qdldl::etree(n_, ptr.data(), rows.data(), parent_.data(),
                            colcount.data(), work.data());
}

Int UpLookingSolver::setup() {
  // KktMatrix stores lower triangle, unpermuted.
  // Transform it into upper triangle and permute it.
  std::vector<Int> ptrT(ptr_.size()), rowsT(rows_.size());
  transpose(ptr_, rows_, ptrT, rowsT);
  std::vector<double> empty_val;
  permuteSym(kkt_.iperm(), ptrT, rowsT, empty_val, false);

  std::vector<Int64> colcount;
  etreeAndCounts(ptrT, rowsT, colcount);

  Int64 nzL{};
  flops_ = 0.0;
  for (Int64 i : colcount) {
    nzL += i;
    flops_ += (i - 1) * (i - 1);
  }

  ptrL_.resize(n_ + 1);
  rowsL_.resize(nzL);
  valL_.resize(nzL);

  counts2Ptr(n_, ptrL_.data(), colcount.data());

  regularisation_.assign(n_, 0.0);

  signs_.assign(n_, 1.0);
  if (kkt_.isAS())
    for (Int i = 0; i < model_.A().num_col_; ++i) signs_[i] = -1.0;
  permuteVectorInverse(signs_, kkt_.iperm());

  return kOk;
}

Int UpLookingSolver::factor(const std::vector<Int>& ptr,
                            const std::vector<Int>& rows,
                            const std::vector<double>& val) {
  // Up-looking factorisation.
  // It computes L one row at a time.
  // If the first k rows have been computed, then we already have a
  // factorisation of M(1:k,1:k) = LDL^T (in Matlab notation). The next row of L
  // is [l^T 1], the next element of D is d, and the next row of M is [b^T a]:
  //
  // [ L   ] [D  ] [L^T l] = [ M  b]
  // [l^T 1] [  d] [    1]   [b^T a]
  //
  // Then
  // - l is found solving LDl = b
  // - d is found solving l^TDl + d = a
  //
  // A must be upper triangular.
  // L is computed as lower triangular. The diagonal of L is used to store D^-1.
  //
  // Returns kErrorFactorise if the time limit is reached; L is then invalid.

  std::vector<HighsBool> mark(n_, false);
  std::vector<Int> stack(n_);
  std::vector<Int> revpattern(n_);
  std::vector<double> x(n_, 0.0);
  std::vector<Int64> next = ptrL_;

  // number of entries of L scanned since the last check of the time limit
  Int64 work = 0;

  for (Int row = 0; row < n_; ++row) {
    if (work >= kUplookTimeCheckWork) {
      work = 0;
      if (control_.timeLimitReached()) return kErrorFactorise;
    }

    double d = HighsExtras::qdldl::factorRow(
        row, ptr.data(), rows.data(), val.data(), parent_.data(),
        ptrL_.data(), rowsL_.data(), valL_.data(), next.data(), mark.data(),
        stack.data(), revpattern.data(), x.data(), &work);

    // diagonal entry
    double old_pivot = d;
    if (signs_[row] > 0)
      d += regul_.dual;
    else
      d -= regul_.primal;

    if (d * signs_[row] < reg_threshold_) {
      d = reg_value_ * signs_[row];
    }

    regularisation_[row] = d - old_pivot;

    valL_[ptrL_[row]] = 1.0 / d;
    rowsL_[next[row]] = row;
    next[row]++;
  }

  return kOk;
}

void UpLookingSolver::solve(std::vector<double>& x) {
  // forward solve
  for (Int i = 0; i < n_; ++i) {
    for (Int el = ptrL_[i] + 1; el < ptrL_[i + 1]; ++el) {
      x[rowsL_[el]] -= valL_[el] * x[i];
    }
  }

  // diag solve
  for (Int i = 0; i < n_; ++i) {
    x[i] *= valL_[ptrL_[i]];
  }

  // backward solve
  for (Int i = n_ - 1; i >= 0; --i) {
    for (Int el = ptrL_[i] + 1; el < ptrL_[i + 1]; ++el) {
      x[i] -= valL_[el] * x[rowsL_[el]];
    }
  }
}

Int UpLookingSolver::factorAS(const std::vector<double>& scaling) {
  assert(!valid_);
  assert(kkt_.isAS());

  Clock clock;
  kkt_.buildASvalues(scaling);
  info_.times[kMatrixValuesTime] += clock.stop();

  // construct permuted upper triangle
  std::vector<Int> ptrT(ptr_.size()), rowsT(rows_.size());
  std::vector<double> valT(val_.size());
  transpose(ptr_, rows_, val_, ptrT, rowsT, valT);
  permuteSym(kkt_.iperm(), ptrT, rowsT, valT, false);

  clock.start();
  if (Int status = factor(ptrT, rowsT, valT)) return status;
  info_.times[kFactoriseTime] += clock.stop();
  info_.factor_number++;

  valid_ = true;
  return kOk;
}
Int UpLookingSolver::factorNE(const std::vector<double>& scaling) {
  assert(!valid_);
  assert(kkt_.isNE());

  Clock clock;
  kkt_.buildNEvalues(scaling);
  info_.times[kMatrixValuesTime] += clock.stop();

  // construct permuted upper triangle
  std::vector<Int> ptrT(ptr_.size()), rowsT(rows_.size());
  std::vector<double> valT(val_.size());
  transpose(ptr_, rows_, val_, ptrT, rowsT, valT);
  permuteSym(kkt_.iperm(), ptrT, rowsT, valT, false);

  clock.start();
  if (Int status = factor(ptrT, rowsT, valT)) return status;
  info_.times[kFactoriseTime] += clock.stop();
  info_.factor_number++;

  valid_ = true;
  return kOk;
}

Int UpLookingSolver::solveAS(const std::vector<double>& rhs_x,
                             const std::vector<double>& rhs_y,
                             std::vector<double>& lhs_x,
                             std::vector<double>& lhs_y) {
  assert(valid_);

  Int n = rhs_x.size();

  Clock clock;
  as_buffer_.resize(rhs_x.size() + rhs_y.size());
  std::copy(rhs_x.begin(), rhs_x.end(), as_buffer_.begin());
  std::copy(rhs_y.begin(), rhs_y.end(), as_buffer_.begin() + n);
  info_.times[kInsertSplitTime] += clock.stop();

  clock.start();
  permuteVectorInverse(as_buffer_, kkt_.iperm());
  solve(as_buffer_);
  permuteVector(as_buffer_, kkt_.iperm());
  info_.times[kSolveTime] += clock.stop();
  info_.solve_number++;
  data_.back().num_solves++;

  clock.start();
  lhs_x.assign(as_buffer_.begin(), as_buffer_.begin() + n);
  lhs_y.assign(as_buffer_.begin() + n, as_buffer_.end());
  info_.times[kInsertSplitTime] += clock.stop();

  return kOk;
}
Int UpLookingSolver::solveNE(const std::vector<double>& rhs,
                             std::vector<double>& lhs) {
  assert(valid_);

  lhs = rhs;
  Clock clock;
  permuteVectorInverse(lhs, kkt_.iperm());
  solve(lhs);
  permuteVector(lhs, kkt_.iperm());
  info_.times[kSolveTime] += clock.stop();
  info_.solve_number++;
  data_.back().num_solves++;

  return kOk;
}

void UpLookingSolver::getReg(std::vector<double>& reg) {
  reg = regularisation_;
  permuteVector(reg, kkt_.iperm());
}

}  // namespace hipo
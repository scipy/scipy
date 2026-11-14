#include "KktMatrix.h"

#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

KktMatrix::KktMatrix(const Model& m, const Regularisation& r, Info& i,
                     const Logger& l, const Options& o)
    : model{m}, regul{r}, info{i}, logger{l}, options{o} {}

Int KktMatrix::buildASstructure() {
  // Build lower triangular structure of the augmented system.
  // Build values of AS that will not change during the iterations.
  Clock clock;

  const HighsSparseMatrix& A = model.A();
  const HighsHessian& Q = model.Q();
  const Int n = A.num_col_;
  const Int m = A.num_row_;
  const Int nzA = A.numNz();
  const Int nzQ = Q.numNz();

  logger.printInfo("Building AS structure\n");

  const Int nzBlock11 = model.qp() ? nzQ : n;

  // AS matrix must fit into HighsInt
  if ((Int64)nzBlock11 + m + nzA >= kHighsIInf) return kErrorOverflow;

  ptrAS.resize(n + m + 1);
  rowsAS.resize(nzBlock11 + nzA + m);
  valAS.resize(nzBlock11 + nzA + m);

  Int next = 0;

  for (Int i = 0; i < n; ++i) {
    // diagonal element
    rowsAS[next] = i;
    next++;

    // column of Q
    if (model.qp()) {
      assert(Q.index_[Q.start_[i]] == i);
      for (Int el = Q.start_[i] + 1; el < Q.start_[i + 1]; ++el) {
        rowsAS[next] = Q.index_[el];
        valAS[next] = -Q.value_[el];  // values of AS that will not change
        ++next;
      }
    }

    // column of A
    for (Int el = A.start_[i]; el < A.start_[i + 1]; ++el) {
      rowsAS[next] = n + A.index_[el];
      valAS[next] = A.value_[el];  // values of AS that will not change
      ++next;
    }

    ptrAS[i + 1] = next;
  }

  // 2,2 block
  for (Int i = 0; i < m; ++i) {
    rowsAS[next] = n + i;
    ++next;
    ptrAS[n + i + 1] = ptrAS[n + i] + 1;
  }

  info.times[kMatrixStructureTime_AS] = clock.stop();

  AS_ready = true;
  return kOk;
}

Int KktMatrix::buildASvalues(const std::vector<double>& scaling) {
  // build AS values that change during iterations.

  assert(AS_ready);
  Clock clock;

  const HighsHessian& Q = model.Q();
  const Int n = model.A().num_col_;

  for (Int i = 0; i < n; ++i) {
    valAS[ptrAS[i]] = scaling.empty() ? -1.0 : -scaling[i];
    if (model.qp()) valAS[ptrAS[i]] -= model.sense() * model.Q().diag(i);
  }

  info.times[kMatrixValuesTime] += clock.stop();

  return kOk;
}

Int KktMatrix::buildNEstructure() {
  // Build lower triangular structure of AAt.
  // This approach uses a column-wise copy of A, a partial row-wise copy and a
  // vector of corresponding indices.
  // NB: A must have sorted columns for this to work
  Clock clock;

  const HighsSparseMatrix& A = model.A();
  const Int n = A.num_col_;
  const Int m = A.num_row_;
  const Int nzA = A.numNz();

  logger.printInfo("Building NE structure\n");

  // create partial row-wise representation without values, and array or
  // corresponding indices between cw and rw representation

  ptrA_rw.assign(m + 1, 0);
  idxA_rw.assign(nzA, 0);

  // pointers of row-start
  for (Int el = 0; el < nzA; ++el) ptrA_rw[A.index_[el] + 1]++;
  for (Int i = 0; i < m; ++i) ptrA_rw[i + 1] += ptrA_rw[i];

  std::vector<Int> temp = ptrA_rw;
  corr_A.assign(nzA, 0);

  // rw-indices and corresponding indices created together
  for (Int col = 0; col < n; ++col) {
    for (Int el = A.start_[col]; el < A.start_[col + 1]; ++el) {
      Int row = A.index_[el];

      corr_A[temp[row]] = el;
      idxA_rw[temp[row]] = col;
      temp[row]++;
    }
  }

  ptrNE.clear();
  rowsNE.clear();

  std::vector<Int> work(m);
  std::atomic<Int> current_nz{0};
  std::atomic<bool> overflow{false};

  auto process_row = [&](Int row, std::vector<char>& is_nz,
                         std::vector<Int>& temp_index,
                         std::vector<Int>& row_space) {
    // go along the entries of the row, and then down each column.
    // this builds the lower triangular part of the row-th column of AAt.

    if (overflow.load(std::memory_order_relaxed)) return;

    Int nz_in_col = 0;

    for (Int el = ptrA_rw[row]; el < ptrA_rw[row + 1]; ++el) {
      Int col = idxA_rw[el];
      Int corr = corr_A[el];

      // for each nonzero in the row, go down corresponding column,
      // starting from current position
      for (Int colEl = corr; colEl < A.start_[col + 1]; ++colEl) {
        Int row2 = A.index_[colEl];

        // row2 is guaranteed to be larger or equal than row
        // (provided that the columns of A are sorted)

        // save information that there is nonzero in position (row2,row).
        if (!is_nz[row2]) {
          is_nz[row2] = true;
          temp_index[nz_in_col] = row2;
          ++nz_in_col;
        }
      }
    }
    // intersection of row with rows below finished.

    // if the total number of nonzeros exceeds the maximum, return error.
    if (current_nz.load(std::memory_order_relaxed) + nz_in_col >=
        NE_nz_limit.load(std::memory_order_relaxed)) {
      overflow.store(true, std::memory_order_relaxed);
      return;
    }
    current_nz.fetch_add(nz_in_col, std::memory_order_relaxed);

    // keep track of column counts
    work[row] = nz_in_col;

    // now assign indices
    for (Int i = 0; i < nz_in_col; ++i) {
      Int index = temp_index[i];
      // push_back is better then reserve, because the final length is not
      // known
      row_space.push_back(index);
      is_nz[index] = false;
    }
  };

  if (options.getParallel(ParallelTechnique::kNEStruct)) {
    std::vector<std::vector<Int>> rowsNE_local(m);

    highs::parallel::for_each(
        0, m,
        [&](Int start, Int end) {
          std::vector<char> is_nz(m, false);
          std::vector<Int> temp_index(m);
          for (Int row = start; row < end; ++row) {
            process_row(row, is_nz, temp_index, rowsNE_local[row]);
          }
        },
        // choose grainsize so that the number of tasks in the for_each loop
        // that execute the function is roughly kParallelNEStructTasks
        std::ceil((double)m / kParallelNEStructTasks));

    if (overflow) {
      info.times[kMatrixStructureTime_NE] = clock.stop();
      return kErrorOverflow;
    }

    ptrNE.resize(m + 1, 0);
    counts2Ptr(m, ptrNE.data(), work.data());

    rowsNE.reserve(ptrNE.back());
    for (const auto& v : rowsNE_local)
      rowsNE.insert(rowsNE.end(), v.begin(), v.end());

  } else {
    rowsNE.reserve(model.nzNElb());
    std::vector<char> is_nz(m, false);
    std::vector<Int> temp_index(m);
    for (Int row = 0; row < m; ++row) {
      process_row(row, is_nz, temp_index, rowsNE);
    }

    if (overflow) {
      info.times[kMatrixStructureTime_NE] = clock.stop();
      return kErrorOverflow;
    }

    ptrNE.resize(m + 1, 0);
    counts2Ptr(m, ptrNE.data(), work.data());
  }

  info.times[kMatrixStructureTime_NE] = clock.stop();

  NE_ready = true;
  return kOk;
}

Int KktMatrix::buildNEvalues(const std::vector<double>& scaling) {
  // given the NE structure already computed, fill in the NE values

  assert(NE_ready);
  Clock clock;

  const HighsSparseMatrix& A = model.A();
  const HighsHessian& Q = model.Q();
  const Int m = A.num_row_;

  valNE.resize(rowsNE.size());

  auto process_row = [&](Int row, std::vector<double>& work) {
    // go along the entries of the row, and then down each column.
    // this builds the lower triangular part of the row-th column of AAt.

    for (Int el = ptrA_rw[row]; el < ptrA_rw[row + 1]; ++el) {
      Int col = idxA_rw[el];
      Int corr = corr_A[el];

      double denom = scaling.empty() ? 1.0 : scaling[col];
      denom += regul.primal;
      if (model.qp()) denom += model.sense() * Q.diag(col);

      const double mult = 1.0 / denom;
      const double row_value = mult * A.value_[corr];

      // for each nonzero in the row, go down corresponding column,
      // starting from current position
      for (Int colEl = corr; colEl < A.start_[col + 1]; ++colEl) {
        Int row2 = A.index_[colEl];

        // row2 is guaranteed to be larger or equal than row
        // (provided that the columns of A are sorted)

        // compute and accumulate value
        double value = row_value * A.value_[colEl];
        work[row2] += value;
      }
    }
    // intersection of row with rows below finished.

    // read from work, using indices of column "row" of AAt
    for (Int el = ptrNE[row]; el < ptrNE[row + 1]; ++el) {
      Int index = rowsNE[el];
      valNE[el] = work[index];
      work[index] = 0.0;
    }
  };

  if (options.getParallel(ParallelTechnique::kNEValues)) {
    highs::parallel::for_each(
        0, m,
        [&](Int start, Int end) {
          std::vector<double> work(m, 0.0);
          for (Int row = start; row < end; ++row) process_row(row, work);
        },
        // choose grainsize so that the number of tasks in the for_each loop
        // that execute the function is roughly kParallelNEValuesTasks
        std::ceil((double)m / kParallelNEValuesTasks));

  } else {
    std::vector<double> work(m, 0.0);
    for (Int row = 0; row < m; ++row) process_row(row, work);
  }

  info.times[kMatrixValuesTime] += clock.stop();

  return kOk;
}

void KktMatrix::freeASmemory() {
  // Give up memory used for AS.
  freeVector(ptrAS);
  freeVector(rowsAS);
  freeVector(valAS);
  AS_ready = false;
}

void KktMatrix::freeNEmemory() {
  // Give up memory used for NE.
  freeVector(ptrNE);
  freeVector(rowsNE);
  freeVector(valNE);
  freeVector(ptrA_rw);
  freeVector(idxA_rw);
  freeVector(corr_A);
  NE_ready = false;
}

Int KktMatrix::n() const {
  if (isNE()) return ptrNE.size() - 1;
  if (isAS()) return ptrAS.size() - 1;
  return -1;
}
Int KktMatrix::nz() const {
  if (isNE()) return rowsNE.size();
  if (isAS()) return rowsAS.size();
  return -1;
}
std::string KktMatrix::nla() const {
  if (isNE()) return kHipoNormalEqString;
  if (isAS()) return kHipoAugmentedString;
  return "empty";
}

bool KktMatrix::isAS() const { return AS_ready && !NE_ready; }
bool KktMatrix::isNE() const { return NE_ready && !AS_ready; }

const std::vector<Int>& KktMatrix::iperm() const { return S.iperm(); }

}  // namespace hipo
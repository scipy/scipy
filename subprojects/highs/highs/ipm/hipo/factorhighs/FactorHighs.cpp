#include "FactorHighs.h"

#include "Analyse.h"
#include "DataCollector.h"
#include "Factorise.h"
#include "HighsExternalApi.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

FHsolver::~FHsolver() {
  if (logger_) {
    data_.printTimes(*logger_);
    data_.printIter(*logger_);
  }
  if (local_logger_ && logger_) delete logger_;
}

void FHsolver::newIter() { data_.append(); }

void FHsolver::setRegularisation(double reg_p, double reg_d) {
  options_.reg_p = reg_p;
  options_.reg_d = reg_d;
}

Int FHsolver::analyse(Symbolic& S, Int n, Int nz, const Int* rows,
                      const Int* ptr, const Int* signs, const Int* perm) {
  Analyse an_obj(n, nz, rows, ptr, signs, options_, logger_, data_, perm);
  return an_obj.run(S);
}

Int FHsolver::factorise(const Symbolic& S, Int n, Int nz, const Int* rows,
                        const Int* ptr, const double* vals) {
  Factorise fact_obj(S, n, nz, rows, ptr, vals, options_, logger_, data_,
                     sn_columns_, &serial_stack_);
  return fact_obj.run(N_);
}

Int FHsolver::solve(double* x, Int k) const { return N_.solve(x, k); }
Int FHsolver::forwardSolve(double* x, Int k) const {
  return N_.forwardSolve(x, k);
}
Int FHsolver::diagSolve(double* x, Int k) const { return N_.diagSolve(x, k); }
Int FHsolver::backwardSolve(double* x, Int k) const {
  return N_.backwardSolve(x, k);
}

void FHsolver::getRegularisation(double* reg) { N_.getReg(reg); }

void FHsolver::setBlockSize(Int nb) {
  if (nb > 0) options_.nb = nb;
}

void FHsolver::setPivoting(bool pivoting) { options_.pivoting = pivoting; }

void FHsolver::setParallel(bool tree, bool node) {
  options_.parallel_tree = tree;
  options_.parallel_node = node;
}

void FHsolver::setParallelSolve(bool forward, bool backward, bool diag) {
  options_.parallel_forward = forward;
  options_.parallel_backward = backward;
  options_.parallel_diag = diag;
}

void FHsolver::setLogger(const Logger* logger, bool use_printf) {
  if (local_logger_ && logger_) delete logger_;
  local_logger_ = false;

  logger_ = logger;
  if (!logger_) {
    logger_ = new Logger(use_printf);
    if (logger_) local_logger_ = true;
  }

#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  if (logger_)
    logger_->printw(
        "Running in debug mode: COLLECTING EXPENSIVE FACTORISATION DATA\n");
#endif
#if HIPO_TIMING_LEVEL > 0
  if (logger_)
    logger_->printw(
        "Running in debug mode: COLLECTING EXPENSIVE TIMING DATA\n");
#endif
}

void FHsolver::inertia(Int& pos, Int& neg, Int& zero, double tol) const {
  N_.inertia(pos, neg, zero, tol);
}

void FHsolver::setOneIndexing(bool one_indexing) {
  options_.one_indexing = one_indexing;
}

void FHsolver::iperm(const Symbolic& S, Int* ip) const {
  const std::vector<Int>& iperm = S.iperm();
  const Int offset = options_.one_indexing ? 1 : 0;
  for (Int i = 0; i < static_cast<Int>(iperm.size()); ++i) {
    ip[i] = iperm[i] + offset;
  }
}

static void getFull(Int n, Int nz, const Int* rows, const Int* ptr,
                    bool one_indexing, std::vector<Int>& full_ptr_v,
                    std::vector<Int>& full_rows_v) {
  std::vector<Int> rows_tmp(nz);
  std::vector<Int> ptr_tmp(n + 1);
  for (Int i = 0; i < nz; ++i) rows_tmp[i] = rows[i];
  for (Int i = 0; i < n + 1; ++i) ptr_tmp[i] = ptr[i];

  if (one_indexing) {
    for (Int& i : rows_tmp) --i;
    for (Int& i : ptr_tmp) --i;
  }

  // metis, amd, rcm require a full copy of the matrix, without diagonal
  // entries
  fullFromLower(n, nz, ptr_tmp.data(), rows_tmp.data(), full_ptr_v,
                full_rows_v);
}

Int FHsolver::reorderMetis(Int n, Int nz, const Int* rows, const Int* ptr,
                           Int* perm, bool full_matrix_0, Int seed) const {
  const Int *full_ptr, *full_rows;
  std::vector<Int> full_ptr_v, full_rows_v;
  if (full_matrix_0) {
    full_ptr = ptr;
    full_rows = rows;
  } else {
    getFull(n, nz, rows, ptr, options_.one_indexing, full_ptr_v, full_rows_v);
    full_ptr = full_ptr_v.data();
    full_rows = full_rows_v.data();
  }

  idx_t options[METIS_NOPTIONS];
  HighsExtras::metis::set_default_options(options);
  options[METIS_OPTION_SEED] = seed;

  options[METIS_OPTION_DBGLVL] = 0;

  // no2hop improves the quality of ordering in general
  options[METIS_OPTION_NO2HOP] = 1;

  std::vector<Int> iperm(n);

  Int status = HighsExtras::metis::nodeND(&n, full_ptr, full_rows, nullptr,
                                          options, perm, iperm.data());

  if (options_.one_indexing)
    for (Int i = 0; i < n; ++i) perm[i]++;

  return status != METIS_OK;
}

Int FHsolver::reorderAmd(Int n, Int nz, const Int* rows, const Int* ptr,
                         Int* perm, bool full_matrix_0) const {
  const Int *full_ptr, *full_rows;
  std::vector<Int> full_ptr_v, full_rows_v;
  if (full_matrix_0) {
    full_ptr = ptr;
    full_rows = rows;
  } else {
    getFull(n, nz, rows, ptr, options_.one_indexing, full_ptr_v, full_rows_v);
    full_ptr = full_ptr_v.data();
    full_rows = full_rows_v.data();
  }

  double control[AMD_CONTROL];
  HighsExtras::amd::set_defaults(control);
  double info[AMD_INFO];

  Int status =
      HighsExtras::amd::order(n, full_ptr, full_rows, perm, control, info);

  if (options_.one_indexing)
    for (Int i = 0; i < n; ++i) perm[i]++;

  return status != AMD_OK;
}

Int FHsolver::reorderRcm(Int n, Int nz, const Int* rows, const Int* ptr,
                         Int* perm, bool full_matrix_0) const {
  const Int *full_ptr, *full_rows;
  std::vector<Int> full_ptr_v, full_rows_v;
  if (full_matrix_0) {
    full_ptr = ptr;
    full_rows = rows;
  } else {
    getFull(n, nz, rows, ptr, options_.one_indexing, full_ptr_v, full_rows_v);
    full_ptr = full_ptr_v.data();
    full_rows = full_rows_v.data();
  }

  Int status =
      HighsExtras::rcm::genrcm(n, full_ptr[n], full_ptr, full_rows, perm);

  if (options_.one_indexing)
    for (Int i = 0; i < n; ++i) perm[i]++;

  return status != 0;
}

}  // namespace hipo
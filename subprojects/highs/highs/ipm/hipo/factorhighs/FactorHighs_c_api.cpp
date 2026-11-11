#include "FactorHighs_c_api.h"

#include "FactorHighs.h"
#include "HighsExternalApi.h"
#include "parallel/HighsParallel.h"

HighsInt FactorHighs_initialise(HighsInt threads) {
  highs::parallel::initialize_scheduler(threads);
  bool hipo_available = HighsExternalApi::isAvailable<HighsExtras::hipo>();
  if (hipo_available) HighsExtras::blas::openblas_set_num_threads(1);
  return !hipo_available;
}

void FactorHighs_terminate(void) { HighsTaskExecutor::shutdown(true); }

void* FactorHighs_create(void) { return new hipo::FHsolver(); }
void FactorHighs_destroy(void* FH) { delete (hipo::FHsolver*)FH; }

void* FactorHighs_symbolic_create(void) { return new hipo::Symbolic(); }
void FactorHighs_symbolic_destroy(void* S) { delete (hipo::Symbolic*)S; }

HighsInt FactorHighs_analyse(void* FH, void* S, HighsInt n, HighsInt nz,
                             const HighsInt* rows, const HighsInt* ptr,
                             const HighsInt* signs, const HighsInt* perm) {
  return ((hipo::FHsolver*)FH)
      ->analyse(*(hipo::Symbolic*)S, n, nz, rows, ptr, signs, perm);
}

HighsInt FactorHighs_factorise(void* FH, const void* S, HighsInt n, HighsInt nz,
                               const HighsInt* rows, const HighsInt* ptr,
                               const double* vals) {
  return ((hipo::FHsolver*)FH)
      ->factorise(*(hipo::Symbolic*)S, n, nz, rows, ptr, vals);
}

HighsInt FactorHighs_solve(void* FH, double* x, HighsInt k) {
  return ((hipo::FHsolver*)FH)->solve(x, k);
}

HighsInt FactorHighs_forwardSolve(void* FH, double* x, HighsInt k) {
  return ((hipo::FHsolver*)FH)->forwardSolve(x, k);
}

HighsInt FactorHighs_diagSolve(void* FH, double* x, HighsInt k) {
  return ((hipo::FHsolver*)FH)->diagSolve(x, k);
}

HighsInt FactorHighs_backwardSolve(void* FH, double* x, HighsInt k) {
  return ((hipo::FHsolver*)FH)->backwardSolve(x, k);
}

void FactorHighs_setRegularisation(void* FH, double reg_p, double reg_d) {
  ((hipo::FHsolver*)FH)->setRegularisation(reg_p, reg_d);
}

void FactorHighs_getRegularisation(void* FH, double* reg) {
  ((hipo::FHsolver*)FH)->getRegularisation(reg);
}

void FactorHighs_newIter(void* FH) { ((hipo::FHsolver*)FH)->newIter(); }

void FactorHighs_setBlockSize(void* FH, HighsInt nb) {
  ((hipo::FHsolver*)FH)->setBlockSize(nb);
}

void FactorHighs_setPivoting(void* FH, HighsInt pivoting) {
  ((hipo::FHsolver*)FH)->setPivoting(pivoting);
}

void FactorHighs_setLogging(void* FH, HighsInt display) {
  ((hipo::FHsolver*)FH)->setLogger(nullptr, display);
}

void FactorHighs_inertia(void* FH, HighsInt* pos, HighsInt* neg, HighsInt* zero,
                         double tol) {
  ((hipo::FHsolver*)FH)->inertia(*pos, *neg, *zero, tol);
}

void FactorHighs_symbolic_print(void* FH, void* S, HighsInt verbose) {
  ((hipo::Symbolic*)S)->print(*(((hipo::FHsolver*)FH)->getLogger()), verbose);
}

void FactorHighs_iperm(void* FH, void* S, HighsInt* ip) {
  ((hipo::FHsolver*)FH)->iperm(*(hipo::Symbolic*)S, ip);
}

void FactorHighs_setOneIndexing(void* FH, HighsInt one_indexing) {
  ((hipo::FHsolver*)FH)->setOneIndexing(one_indexing);
}

HighsInt FactorHighs_reorderMetis(void* FH, HighsInt n, HighsInt nz,
                                  const HighsInt* rows, const HighsInt* ptr,
                                  HighsInt* perm) {
  return ((hipo::FHsolver*)FH)->reorderMetis(n, nz, rows, ptr, perm, 0, 0);
}
HighsInt FactorHighs_reorderAmd(void* FH, HighsInt n, HighsInt nz,
                                const HighsInt* rows, const HighsInt* ptr,
                                HighsInt* perm) {
  return ((hipo::FHsolver*)FH)->reorderAmd(n, nz, rows, ptr, perm, 0);
}
HighsInt FactorHighs_reorderRcm(void* FH, HighsInt n, HighsInt nz,
                                const HighsInt* rows, const HighsInt* ptr,
                                HighsInt* perm) {
  return ((hipo::FHsolver*)FH)->reorderRcm(n, nz, rows, ptr, perm, 0);
}
#include "Numeric.h"

#include "DataCollector.h"
#include "FactorHighsSettings.h"
#include "HybridSolveHandler.h"
#include "ReturnValues.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"
#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"

namespace hipo {

void Numeric::computeAnySwaps() {
  // compute which blocks of columns require swaps
  any_swaps_.resize(S_->sn());
  const Int nb = options_->nb;
  for (Int sn = 0; sn < S_->sn(); ++sn) {
    const Int sn_size = S_->snStart(sn + 1) - S_->snStart(sn);
    const Int n_blocks = (sn_size - 1) / nb + 1;
    any_swaps_[sn].assign(n_blocks, 0);

    for (Int j = 0; j < n_blocks; ++j) {
      const Int jb = std::min(nb, sn_size - nb * j);
      for (Int i = 0; i < jb; ++i) {
        if (swaps_[sn][nb * j + i] != i) {
          any_swaps_[sn][j] = 1;
          break;
        }
      }
    }
  }
}

Int Numeric::prepare() {
  if (!sn_columns_ || !S_ || !data_ || !options_) return kRetInvalidPointer;
  SH_.reset(new HybridSolveHandler(*S_, *sn_columns_, swaps_, any_swaps_,
                                   pivot_2x2_, *data_, *options_));
  if (!SH_) return kRetGeneric;
  if (options_->pivoting) computeAnySwaps();
  return kRetOk;
}

Int Numeric::solve(double* x) const {
  // Return the number of solves performed

  if (!SH_) return kRetGeneric;

  HIPO_CLOCK_CREATE;

  // permute rhs
  HIPO_CLOCK_START(2);
  permuteVectorInverse(x, S_->iperm());
  HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare);

  // solve
  HIPO_CLOCK_START(2);
  SH_->forwardSolve(x);
  SH_->diagSolve(x);
  SH_->backwardSolve(x);
  HIPO_CLOCK_STOP(2, *data_, kTimeSolveSolve);

  // unpermute solution
  HIPO_CLOCK_START(2);
  permuteVector(x, S_->iperm());
  HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare);

  HIPO_CLOCK_STOP(1, *data_, kTimeSolve);
  return kRetOk;
}

Int Numeric::forwardSolve(double* x) const {
  if (!SH_) return kRetGeneric;
  permuteVectorInverse(x, S_->iperm());
  SH_->forwardSolve(x);
  return kRetOk;
}
Int Numeric::diagSolve(double* x) const {
  if (!SH_) return kRetGeneric;
  SH_->diagSolve(x);
  return kRetOk;
}
Int Numeric::backwardSolve(double* x) const {
  if (!SH_) return kRetGeneric;
  SH_->backwardSolve(x);
  permuteVector(x, S_->iperm());
  return kRetOk;
}

#define SOLVE_MULTIPLE(f)                                        \
  if (k == 1)                                                    \
    return f(x);                                                 \
  else {                                                         \
    highs::parallel::TaskGroup tg;                               \
    const Int n = S_->size();                                    \
    std::atomic<bool> fail{false};                               \
    for (Int i = 0; i < k; ++i) {                                \
      tg.spawn([=, &fail]() {                                    \
        Int status = f(&x[i * n]);                               \
        if (status) fail.store(true, std::memory_order_relaxed); \
      });                                                        \
    }                                                            \
    tg.taskWait();                                               \
    return fail;                                                 \
  }

Int Numeric::solve(double* x, Int k) const { SOLVE_MULTIPLE(solve); }
Int Numeric::forwardSolve(double* x, Int k) const {
  SOLVE_MULTIPLE(forwardSolve);
}
Int Numeric::diagSolve(double* x, Int k) const { SOLVE_MULTIPLE(diagSolve); }
Int Numeric::backwardSolve(double* x, Int k) const {
  SOLVE_MULTIPLE(backwardSolve);
}

void Numeric::getReg(double* reg) {
  std::memcpy(reg, total_reg_.data(), total_reg_.size() * sizeof(double));
}

void Numeric::inertia(Int& pos, Int& neg, Int& zero, double tol) const {
  if (!SH_) return;
  SH_->inertia(pos, neg, zero, tol);
}

}  // namespace hipo
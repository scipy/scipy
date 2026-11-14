#include "DgemmParallel.h"

#include "CallAndTimeBlas.h"
#include "FactorHighsSettings.h"
#include "parallel/HighsParallel.h"

namespace hipo {

dgemmParalleliser::dgemmParalleliser(const double* P, const double* R,
                                     double* Q, Int col, Int jb,
                                     DataCollector& data)
    : P_{P}, R_{R}, Q_{Q}, col_{col}, jb_{jb}, data_{data} {}

void dgemmParalleliser::run(Int start, Int end, double beta) const {
  callAndTime_dgemm('T', 'N', col_, end - start, jb_, -1.0, P_, jb_,
                    &R_[start * jb_], jb_, beta, &Q_[start * col_], col_,
                    data_);
}

void dgemmParallel(const double* P, const double* R, double* Q, Int col, Int jb,
                   Int row, Int nb, double beta, DataCollector& data) {
  // if there is enough work to be done, parallelise
  if (col >= nb / 2 && jb >= nb / 2 && row >= kBlockParallelThreshold * nb) {
    dgemmParalleliser gemmP(P, R, Q, col, jb, data);
    dgemmParalleliser* pt = &gemmP;

    // I need to use an object to call gemm, otherwise the task is too large and
    // static_assert in the parallel deque fails.
    highs::parallel::for_each(
        0, row, [pt, beta](Int start, Int end) { pt->run(start, end, beta); },
        kBlockGrainSize * jb);
  } else {
    callAndTime_dgemm('T', 'N', col, row, jb, -1.0, P, jb, R, jb, beta, Q, col,
                      data);
  }
}

}  // namespace hipo
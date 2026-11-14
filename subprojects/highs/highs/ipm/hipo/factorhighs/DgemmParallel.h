#ifndef FACTORHIGHS_DGEMM_PARALLEL_H
#define FACTORHIGHS_DGEMM_PARALLEL_H

#include "DataCollector.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// parallelise dgemm for use within factorisation
// Performs Q <- Q - R P^T in hybrid format.
// Parallelised over the rows of R and Q.
class dgemmParalleliser {
  const double* P_;
  const double* R_;
  double* Q_;
  const Int col_;
  const Int jb_;
  DataCollector& data_;

 public:
  dgemmParalleliser(const double* P, const double* R, double* Q, Int col,
                    Int jb, DataCollector& data);

  void run(Int start, Int end, double beta) const;
};

void dgemmParallel(const double* P, const double* R, double* Q, Int col, Int jb,
                   Int row, Int nb, double beta, DataCollector& data);

}  // namespace hipo

#endif
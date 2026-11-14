#ifndef HIPO_KRYLOV_METHODS_H
#define HIPO_KRYLOV_METHODS_H

#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Abstract class for matrices inside of a Krylov method
class AbstractMatrix {
 public:
  virtual void apply(std::vector<double>& x) const = 0;
};

// =================
// Krylov solvers
// =================

Int Gmres(const AbstractMatrix* M, const AbstractMatrix* P,
          const std::vector<double>& b, std::vector<double>& x, double tol,
          Int maxit);

Int Cg(const AbstractMatrix* M, const AbstractMatrix* P,
       const std::vector<double>& b, std::vector<double>& x, double tol,
       Int maxit);

}  // namespace hipo

#endif
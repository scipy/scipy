#ifndef FACTORHIGHS_KRYLOV_METHODS_IPM_H
#define FACTORHIGHS_KRYLOV_METHODS_IPM_H

#include <vector>

#include "Numeric.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/KrylovMethods.h"
#include "util/HighsSparseMatrix.h"

namespace hipo {

// Class to perform matrix-vector products with ipm matrix
class IpmMatrix : public AbstractMatrix {
  const HighsSparseMatrix* A_ = nullptr;
  const std::vector<double>* scaling_ = nullptr;
  bool use_as_ = false;

 public:
  void reset(const HighsSparseMatrix& A, const std::vector<double>& scaling,
             bool use_as = false);
  void apply(std::vector<double>& x) const override;
};

// Class to perform solves with the ipm factorisation
class IpmFactor : public AbstractMatrix {
  const Numeric& N_;

 public:
  IpmFactor(const Numeric& N);
  void apply(std::vector<double>& x) const override;
};

// Class to apply diagonal preconditioner
class NeDiagPrec : public AbstractMatrix {
  std::vector<double> diag;

 public:
  void reset(const HighsSparseMatrix& A, const std::vector<double>& scaling);
  void apply(std::vector<double>& x) const override;
};

}  // namespace hipo

#endif
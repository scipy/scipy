// Copyright (c) 2018-2019 ERGO-Code. See license.txt for license.

#ifndef IPX_DIAGONAL_PRECOND_H_
#define IPX_DIAGONAL_PRECOND_H_

#include "linear_operator.h"
#include "model.h"
#include "sparse_matrix.h"

namespace ipx {

// DiagonalPrecond provides inverse operations with the diagonal matrix
//
//   diag(AI*W*AI').                                (1)
//
// Here AI is the m-by-(n+m) matrix defined by the model, and W is a diagonal
// (weight) matrix that is provided by the user.

class DiagonalPrecond : public LinearOperator {
public:
    // Constructor stores a reference to the model. No data is copied. The model
    // must be valid as long as the preconditioner is used.
    explicit DiagonalPrecond(const Model& model);

    // Factorizes the preconditioner. W must either hold n+m entries, or be
    // NULL, in which case the first n entries are assumed 1.0 and the last
    // m entries are assumed 0.0.
    void Factorize(const double* W, Info* info);

    // Returns computation time for calls to Apply() since last reset_time().
    double time() const;
    void reset_time();

private:
    void _Apply(const Vector& rhs, Vector& lhs, double* rhs_dot_lhs) override;

    const Model& model_;
    bool factorized_{false};    // preconditioner factorized?
    Vector diagonal_;           // diagonal of normal matrix
    double time_{0.0};
};

}  // namespace ipx

#endif  // IPX_DIAGONAL_PRECOND_H_

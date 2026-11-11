// Copyright (c) 2018-2019 ERGO-Code. See license.txt for license.

#ifndef IPX_KKT_SOLVER_DIAG_H_
#define IPX_KKT_SOLVER_DIAG_H_

#include "ipm/ipx/control.h"
#include "ipm/ipx/diagonal_precond.h"
#include "ipm/ipx/kkt_solver.h"
#include "ipm/ipx/model.h"
#include "ipm/ipx/normal_matrix.h"

namespace ipx {

// KKTSolverDiag implements a KKT solver that applies the Conjugate Residuals
// method with diagonal preconditioning to the normal equations. If the (1,1)
// block of the KKT matrix is not positive definite, regularization is applied.
//
// In the call to Factorize() @iterate is allowed to be NULL, in which case the
// (1,1) block of the KKT matrix is the identity matrix.

class KKTSolverDiag : public KKTSolver {
public:
    KKTSolverDiag(const Control& control, const Model& model);

    Int maxiter() const { return maxiter_; }
    void maxiter(Int new_maxiter) { maxiter_ = new_maxiter; }

private:
    void _Factorize(Iterate* iterate, Info* info) override;
    void _Solve(const Vector& a, const Vector& b, double tol,
                Vector& x, Vector& y, Info* info) override;
    Int _iter() const override { return iter_; };

    const Control& control_;
    const Model& model_;
    NormalMatrix normal_matrix_;
    DiagonalPrecond precond_;

    Vector W_;               // diagonal matrix in AI*W*AI'
    Vector resscale_;        // residual scaling factors for CR termination test
    bool factorized_{false}; // KKT matrix factorized?
    Int maxiter_{-1};
    Int iter_{0};               // # CR iterations since last Factorize()
};

}  // namespace ipx

#endif  // IPX_KKT_SOLVER_DIAG_H_

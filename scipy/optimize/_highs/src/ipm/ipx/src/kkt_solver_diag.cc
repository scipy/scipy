// Copyright (c) 2018-2019 ERGO-Code. See license.txt for license.

#include "kkt_solver_diag.h"
#include <cassert>
#include <cmath>
#include "conjugate_residuals.h"

namespace ipx {

KKTSolverDiag::KKTSolverDiag(const Control& control, const Model& model) :
    control_(control), model_(model), normal_matrix_(model), precond_(model) {
    Int m = model_.rows();
    Int n = model_.cols();
    W_.resize(m+n);
    resscale_.resize(m);
}

void KKTSolverDiag::_Factorize(Iterate* pt, Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    iter_ = 0;
    factorized_ = false;

    if (pt) {
        const Vector& xl = pt->xl();
        const Vector& xu = pt->xu();
        const Vector& zl = pt->zl();
        const Vector& zu = pt->zu();

        // Build matrix W for AI*W*AI'. For free variables set W[j] to
        // 1.0/regval, where regval is a regularization value. regval is chosen
        // as the minimum of the complementarity measure mu and the smallest
        // nonzero diagonal entry of the (1,1) block of the KKT matrix.
        double regval = pt->mu();
        for (Int j = 0; j < n+m; j++) {
            assert(xl[j] > 0.0);
            assert(xu[j] > 0.0);
            double g = zl[j]/xl[j] + zu[j]/xu[j];
            assert(std::isfinite(g));
            if (g != 0.0 && g < regval)
                regval = g;
            W_[j] = 1.0 / g;        // infinity if g is zero
        }
        for (Int j = 0; j < n+m; j++) {
            if (std::isinf(W_[j]))
                W_[j] = 1.0 / regval;
            assert(std::isfinite(W_[j]));
            assert(W_[j] > 0.0);
        }
    } else {
        W_ = 1.0;
    }

    // Residual scaling factors for termination test of CR method (see below).
    for (Int i = 0; i < m; i++)
        resscale_[i] = 1.0 / std::sqrt(W_[n+i]);

    // Build normal matrix and preconditioner.
    normal_matrix_.Prepare(&W_[0]);
    precond_.Factorize(&W_[0], info);
    if (info->errflag)
        return;

    factorized_ = true;
}

// Reduces the KKT system
//
//   [ W^{-1}  AI' ] (x) = (a) + (res)
//   [ AI       0  ] (y)   (b)   ( 0 )
//
// to normal equations
//
//   C * y := (AI*W*AI') * y = -b + AI*W*(a+res)
//
// and solves by the CR method. The solution to the KKT system is recovered so
// that the first n entries of res are zero. Therefore the residual in the
// normal equations is W[B]*res[B], where B is the slack basis. By multiplying
// by resscale, the CR method termination criterion tests the condition required
// from the KKT solver (see kkt_solver.h).
//
void KKTSolverDiag::_Solve(const Vector& a, const Vector& b, double tol,
                            Vector& x, Vector& y, Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const SparseMatrix& AI = model_.AI();
    assert(factorized_);

    // Compose right-hand side AI*W*a-b.
    Vector rhs = -b;
    for (Int j = 0; j < n+m; j++)
        ScatterColumn(AI, j, W_[j]*a[j], rhs);

    // Solve normal equations.
    y = 0.0;
    normal_matrix_.reset_time();
    precond_.reset_time();
    ConjugateResiduals cr(control_);
    cr.Solve(normal_matrix_, precond_, rhs, tol, &resscale_[0], maxiter_, y);
    info->errflag = cr.errflag();
    info->kktiter1 += cr.iter();
    info->time_cr1 += cr.time();
    info->time_cr1_AAt += normal_matrix_.time();
    info->time_cr1_pre += precond_.time();
    iter_ += cr.iter();

    // Recover solution to KKT system.
    for (Int i = 0; i < m; i++)
        x[n+i] = b[i];
    for (Int j = 0; j < n; j++) {
        double aty = DotColumn(AI, j, y);
        x[j] = W_[j] * (a[j]-aty);
        for (Int p = AI.begin(j); p < AI.end(j); p++) {
            Int i = AI.index(p);
            x[n+i] -= x[j] * AI.value(p);
        }
    }
}

}  // namespace ipx

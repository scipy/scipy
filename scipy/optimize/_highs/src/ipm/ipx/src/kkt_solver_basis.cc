// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "kkt_solver_basis.h"
#include <cassert>
#include <cmath>
#include "conjugate_residuals.h"
#include "maxvolume.h"
#include "starting_basis.h"

namespace ipx {

KKTSolverBasis::KKTSolverBasis(const Control& control, Basis& basis)
    : control_(control), model_(basis.model()), basis_(basis),
      splitted_normal_matrix_(model_) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    colscale_.resize(n+m);
}

void KKTSolverBasis::_Factorize(Iterate* iterate, Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    info->errflag = 0;
    factorized_ = false;
    iter_ = 0;
    basis_changes_ = 0;

    for (Int j = 0; j < n+m; j++)
        colscale_[j] = iterate->ScalingFactor(j);

    // Remove degenerate variables unless the primal objective is smaller than
    // the dual objective. In the latter case the model might be infeasible or
    // unbounded. In order not to affect an infeasibility test in the IPM (not
    // yet implemented) we do not remove variables if the model looks
    // infeasible.
    if (iterate->pobjective() >= iterate->dobjective()) {
        DropPrimal(iterate, info);
        if (info->errflag)
            return;
        DropDual(iterate, info);
        if (info->errflag)
            return;
    }

    // Run maxvolume ("Russian algorithm").
    Maxvolume maxvol(control_);
    if (control_.update_heuristic() == 0) {
        info->errflag = maxvol.RunSequential(&colscale_[0], basis_);
    } else {
        info->errflag = maxvol.RunHeuristic(&colscale_[0], basis_);
    }
    info->updates_ipm += maxvol.updates();
    info->time_maxvol += maxvol.time();
    basis_changes_ += maxvol.updates();
    if (info->errflag)
        return;

    // Refactorize and build preconditioned normal matrix.
    if (!basis_.FactorizationIsFresh()) {
        info->errflag = basis_.Factorize();
        if (info->errflag)
            return;
    }
    splitted_normal_matrix_.Prepare(basis_, &colscale_[0]);

    factorized_ = true;
}

// Reduces the KKT system to preconditioned normal equations while taking free
// variables (which must be basic) into account. See [1, Section 6.4].
//
// [1] L. Schork, "Basis Preconditioning in Interior Point Methods", PhD thesis
//     (2018)
//
void KKTSolverBasis::_Solve(const Vector& a, const Vector& b, double tol,
                             Vector& x, Vector& y, Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const SparseMatrix& AI = model_.AI();
    Vector rhs(m);      // unpermuted right-hand side
    Vector work(m);
    info->errflag = 0;
    assert(factorized_);

    // Compute work = inverse(B')*v, where v[p] = a[basis[p]] if variable
    // basis[p] is free, and v[p] = 0 otherwise.
    Int num_free = 0;
    work = 0.0;
    for (Int p = 0; p < m; p++) {
        Int j = basis_[p];
        if (basis_.StatusOf(j) == Basis::BASIC_FREE) {
            work[p] = a[j];
            num_free++;
        }
    }
    if (num_free > 0)
        basis_.SolveDense(work, work, 'T');

    // Compute rhs = inverse(B)*(N*D2[nonbasic]*(a[nonbasic]-N'*work)).
    rhs = 0.0;
    if (num_free > 0) {
        for (Int j = 0; j < n+m; j++) {
            if (basis_.StatusOf(j) == Basis::NONBASIC) {
                double d2 = colscale_[j] * colscale_[j];
                double alpha = a[j] - DotColumn(AI, j, work);
                alpha *= d2;
                assert(std::isfinite(alpha));
                ScatterColumn(AI, j, alpha, rhs);
            }
        }
    } else {
        for (Int j = 0; j < n+m; j++) {
            if (basis_.StatusOf(j) == Basis::NONBASIC) {
                double d2 = colscale_[j] * colscale_[j];
                double alpha = d2 * a[j];
                assert(std::isfinite(alpha));
                ScatterColumn(AI, j, alpha, rhs);
            }
        }
    }
    basis_.SolveDense(rhs, rhs, 'N');

    // Compute work = inverse(B)*b.
    basis_.SolveDense(b, work, 'N');

    // Build rhs[p] = (rhs[p]-work[p])/D[j] + D[j]*a[j], where j = basis[p]
    // is not a free variable, and rhs[p] = 0 otherwise.
    for (Int p = 0; p < m; p++) {
        Int j = basis_[p];
        if (basis_.StatusOf(j) == Basis::BASIC) {
            double d = colscale_[j];
            rhs[p] = (rhs[p]-work[p]) / d + a[j] * d;
            assert(std::isfinite(rhs[p]));
        } else {
            assert(basis_.StatusOf(j) == Basis::BASIC_FREE);
            rhs[p] = 0.0;
        }
    }

    // Build permuted rhs in work.
    const Int* colperm = splitted_normal_matrix_.colperm();
    for (Int k = 0; k < m; k++)
        work[k] = rhs[colperm[k]];

    // Solve normal equations.
    splitted_normal_matrix_.reset_time();
    Vector lhs = std::move(rhs); // don't need rhs any more
    lhs = 0.0;
    ConjugateResiduals cr(control_);
    cr.Solve(splitted_normal_matrix_, work, tol, nullptr, maxiter_, lhs);
    info->errflag = cr.errflag();
    info->kktiter2 += cr.iter();
    info->time_cr2 += cr.time();
    info->time_cr2_NNt += splitted_normal_matrix_.time_NNt();
    info->time_cr2_B += splitted_normal_matrix_.time_B();
    info->time_cr2_Bt += splitted_normal_matrix_.time_Bt();
    iter_ += cr.iter();

    // Permute back solution to normal equations.
    for (Int k = 0; k < m; k++)
        y[colperm[k]] = lhs[k];

    // Recover dual solution to KKT system.
    for (Int p = 0; p < m; p++) {
        Int j = basis_[p];
        if (basis_.StatusOf(j) == Basis::BASIC) {
            y[p] /= colscale_[j];
            assert(std::isfinite(y[p]));
        } else {
            assert(basis_.StatusOf(j) == Basis::BASIC_FREE);
            assert(y[p] == 0.0);
            y[p] = a[j];        // slot in solution to free basic variable
        }
    }
    basis_.SolveDense(y, y, 'T');

    // Compute x[nonbasic] and work = b - N*x[nonbasic].
    work = b;
    for (Int j = 0; j < n+m; j++) {
        double xj = 0.0;
        if (basis_.StatusOf(j) == Basis::NONBASIC) {
            xj = a[j] - DotColumn(AI, j, y);
            xj *= colscale_[j] * colscale_[j];
            assert(std::isfinite(xj));
            ScatterColumn(AI, j, -xj, work);
        }
        x[j] = xj;
    }

    // Compute x[basic].
    basis_.SolveDense(work, work, 'N');
    for (Int p = 0; p < m; p++)
        x[basis_[p]] = work[p];
}

void KKTSolverBasis::DropPrimal(Iterate* iterate, Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& xl = iterate->xl();
    const Vector& xu = iterate->xu();
    const Vector& zl = iterate->zl();
    const Vector& zu = iterate->zu();
    IndexedVector btran(m), row(n+m);
    const double drop_primal = control_.ipm_drop_primal();
    const double volume_tol = 2.0;
    info->errflag = 0;

    std::vector<Int> candidates;
    for (Int p = 0; p < m; p++) {
        Int jb = basis_[p];
        if (basis_.StatusOf(jb) != Basis::BASIC) // ignore free variables
            continue;
        assert(std::isfinite(xl[jb]) || std::isfinite(xu[jb]));
        assert(xl[jb] > 0.0);
        assert(xu[jb] > 0.0);
        double xj, zj;          // choose which bound is nearer
        if (xl[jb] <= xu[jb]) {
            xj = xl[jb];
            zj = zl[jb];
        } else {
            xj = xu[jb];
            zj = zu[jb];
        }
        if (xj < 0.01*zj && xj <= drop_primal)
            candidates.push_back(jb);
    }
    if (candidates.empty())
        return;

    // Maintain a copy of the inverse scaling factors of basic variables for
    // faster access.
    Vector invscale_basic(m);
    for (Int p = 0; p < m; p++) {
        Int j = basis_[p];
        invscale_basic[p] = 1.0 / colscale_[j];
        assert(std::isfinite(invscale_basic[p]));
        assert(invscale_basic[p] >= 0.0);
    }

    while (!candidates.empty()) {
        Int jb = candidates.back();
        Int p = basis_.PositionOf(jb);
        assert(p >= 0);
        // Pivot jb out of the basis if the volume increase sufficiently.
        const double s = invscale_basic[p];
        basis_.TableauRow(jb, btran, row, true);
        Int jmax = -1;
        double vmax = volume_tol;
        auto search_pivot = [&](Int j, double pivot) {
            pivot = std::abs(pivot);
            if (pivot > kPivotZeroTol) {
                double v = pivot * colscale_[j] * s;
                if (v > vmax) {
                    vmax = v;
                    jmax = j;
                }
            }
        };
        for_each_nonzero(row, search_pivot);
        if (jmax >= 0) {
            // Pivot jb out of the basis.
            double pivot = row[jmax];
            if (std::abs(pivot) < 1e-3)
                control_.Debug(3)
                    << " |pivot| = " << sci2(std::abs(pivot))
                    << " (primal basic variable close to bound)\n";
            assert(basis_.StatusOf(jmax) == Basis::NONBASIC);
            bool exchanged;
            info->errflag = basis_.ExchangeIfStable(jb, jmax, pivot, 1,
                                                    &exchanged);
            if (info->errflag)
                return;
            if (!exchanged)     // factorization was unstable, try again
                continue;
            invscale_basic[p] = 1.0 / colscale_[jmax];
            assert(std::isfinite(invscale_basic[p]));
            assert(invscale_basic[p] >= 0.0);
            info->updates_ipm++;
            basis_changes_++;
        } else {
            // Make variable jb "implied" at a bound.
            if (zl[jb]/xl[jb] > zu[jb]/xu[jb])
                iterate->make_implied_lb(jb);
            else
                iterate->make_implied_ub(jb);
            basis_.FreeBasicVariable(jb);
            invscale_basic[p] = 0.0;
            colscale_[jb] = INFINITY;
            info->primal_dropped++;
        }
        candidates.pop_back();
    }
}

void KKTSolverBasis::DropDual(Iterate* iterate, Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& xl = iterate->xl();
    const Vector& xu = iterate->xu();
    const Vector& zl = iterate->zl();
    const Vector& zu = iterate->zu();
    IndexedVector ftran(m);
    const double drop_dual = control_.ipm_drop_dual();
    const double volume_tol = 2.0;
    info->errflag = 0;

    std::vector<Int> candidates;
    for (Int jn = 0; jn < n+m; jn++) {
        if (basis_.StatusOf(jn) != Basis::NONBASIC)
            continue;
        assert(std::isfinite(xl[jn]) || std::isfinite(xu[jn]));
        assert(xl[jn] > 0.0);
        assert(xu[jn] > 0.0);
        assert(zl[jn] > 0.0 || zu[jn] > 0.0);
        double xj, zj;
        if (zl[jn] >= zu[jn]) { // choose larger dual variable
            xj = xl[jn];
            zj = zl[jn];
        } else {
            xj = xu[jn];
            zj = zu[jn];
        }
        if (zj < 0.01*xj && zj <= drop_dual)
            candidates.push_back(jn);
    }
    if (candidates.empty())
        return;

    // Maintain a copy of the inverse scaling factors of basic variables for
    // faster access.
    Vector invscale_basic(m);
    for (Int p = 0; p < m; p++) {
        Int j = basis_[p];
        invscale_basic[p] = 1.0 / colscale_[j];
        assert(std::isfinite(invscale_basic[p]));
        assert(invscale_basic[p] >= 0.0);
    }

    while (!candidates.empty()) {
        Int jn = candidates.back();
        // Pivot jn into the basis if volume increases sufficiently.
        const double s = colscale_[jn];
        basis_.SolveForUpdate(jn, ftran);
        Int pmax = -1;
        double vmax = volume_tol;
        auto search_pivot = [&](Int p, double pivot) {
            pivot = std::abs(pivot);
            if (pivot > kPivotZeroTol) {
                double v = pivot * invscale_basic[p] * s;
                if (v > vmax) {
                    vmax = v;
                    pmax = p;
                }
            }
        };
        for_each_nonzero(ftran, search_pivot);
        if (pmax >= 0) {
            double pivot = ftran[pmax];
            if (std::abs(pivot) < 1e-3)
                control_.Debug(3)
                    << " |pivot| = " << sci2(std::abs(pivot))
                    << " (dual nonbasic variable close to zero)\n";
            Int jb = basis_[pmax];
            // Pivot jn into the basis.
            assert(basis_.StatusOf(jb) == Basis::BASIC);
            bool exchanged;
            info->errflag = basis_.ExchangeIfStable(jb, jn, pivot, -1,
                                                    &exchanged);
            if (info->errflag)
                return;
            if (!exchanged)     // factorization was unstable, try again
                continue;
            invscale_basic[pmax] = 1.0 / colscale_[jn];
            assert(std::isfinite(invscale_basic[pmax]));
            assert(invscale_basic[pmax] >= 0.0);
            info->updates_ipm++;
            basis_changes_++;
        } else {
            // Make variable jn "fixed" at its current value.
            iterate->make_fixed(jn);
            basis_.FixNonbasicVariable(jn);
            colscale_[jn] = 0.0;
            info->dual_dropped++;
        }
        candidates.pop_back();
    }
}

}  // namespace ipx

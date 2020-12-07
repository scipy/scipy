// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "conjugate_residuals.h"
#include <algorithm>
#include <cmath>
#include "timer.h"
#include "utils.h"

namespace ipx {

ConjugateResiduals::ConjugateResiduals(const Control& control) :
    control_(control) {}

void ConjugateResiduals::Solve(LinearOperator& C, const Vector& rhs,
                              double tol, const double* resscale, Int maxiter,
                              Vector& lhs) {
    const Int m = rhs.size();
    Vector residual(m);  // rhs - C*lhs
    Vector step(m);      // update to lhs
    Vector Cresidual(m); // C * residual
    Vector Cstep(m);     // C * step
    double cdot = 0.0;   // dot product from C.Apply
    Timer timer;

    errflag_ = 0;
    iter_ = 0;
    time_ = 0.0;
    if (maxiter < 0)
        maxiter = m+100;

    // Initialize residual, step and Cstep.
    if (Infnorm(lhs) == 0.0) {
        residual = rhs;         // saves a matrix-vector op
    } else {
        C.Apply(lhs, residual, nullptr);
        residual = rhs-residual;
    }
    C.Apply(residual, Cresidual, &cdot);
    step = residual;
    Cstep = Cresidual;

    while (true) {
        // Termination check.
        double resnorm = 0.0;
        if (resscale)
            for (Int i = 0; i < m; i++)
                resnorm = std::max(resnorm, std::abs(resscale[i]*residual[i]));
        else
            resnorm = Infnorm(residual);
        if (resnorm <= tol)
            break;
        if (iter_ == maxiter) {
            control_.Debug(3)
                << " CR method not converged in " << maxiter << " iterations."
                << " residual = " << sci2(resnorm) << ','
                << " tolerance = " << sci2(tol) << '\n';
            errflag_ = IPX_ERROR_cr_iter_limit;
            break;
        }
        if (cdot <= 0.0) {
            errflag_ = IPX_ERROR_cr_matrix_not_posdef;
            break;
        }

        // Update lhs, residual and Cresidual.
        const double denom = Dot(Cstep,Cstep);
        const double alpha = cdot/denom;
        if (!std::isfinite(alpha)) {
            errflag_ = IPX_ERROR_cr_inf_or_nan;
            break;
        }
        lhs += alpha*step;
        residual -= alpha*Cstep;
        double cdotnew;
        C.Apply(residual, Cresidual, &cdotnew);

        // Update step and Cstep.
        const double beta = cdotnew/cdot;
        step = residual + beta*step;
        Cstep = Cresidual + beta*Cstep;
        cdot = cdotnew;

        iter_++;
        if ((errflag_ = control_.InterruptCheck()) != 0)
            break;
    }
    time_ = timer.Elapsed();
}

void ConjugateResiduals::Solve(LinearOperator& C, LinearOperator& P,
                              const Vector& rhs, double tol,
                              const double* resscale, Int maxiter, Vector& lhs){
    const Int m = rhs.size();
    Vector residual(m);   // rhs - C*lhs
    Vector sresidual(m);  // preconditioned residual
    Vector step(m);       // update to lhs
    Vector Csresidual(m); // C * sresidual
    Vector Cstep(m);      // C * step
    double cdot = 0.0;    // dot product from C.Apply
    Timer timer;

    // resnorm_precond_system = residual' * P * residual.
    // This quantity is minimized by the preconditioned CR method over a Krylov
    // subspace. Hence (in theory) must decrease strictly monotonically. If it
    // does not, then the method stagnated due to round-off errors. This
    // happened in a few cases with augmented diagonal preconditioning (i.e.
    // diagonal preconditioning with dense columns treated as low-rank update)
    // because operations with P were not sufficiently accurate.
    double resnorm_precond_system = 0.0;

    errflag_ = 0;
    iter_ = 0;
    time_ = 0.0;
    if (maxiter < 0)
        maxiter = m+100;

    // Initialize residual, sresidual, step and Cstep.
    if (Infnorm(lhs) == 0.0) {
        residual = rhs;         // saves a matrix-vector op
    } else {
        C.Apply(lhs, residual, nullptr);
        residual = rhs-residual;
    }
    P.Apply(residual, sresidual, &resnorm_precond_system);
    C.Apply(sresidual, Csresidual, &cdot);
    step = sresidual;
    Cstep = Csresidual;

    while (true) {
        // Termination check.
        double resnorm = 0.0;
        if (resscale)
            for (Int i = 0; i < m; i++)
                resnorm = std::max(resnorm, std::abs(resscale[i]*residual[i]));
        else
            resnorm = Infnorm(residual);
        if (resnorm <= tol)
            break;
        if (iter_ == maxiter) {
            control_.Debug(3)
                << " PCR method not converged in " << maxiter << " iterations."
                << " residual = " << sci2(resnorm) << ','
                << " tolerance = " << sci2(tol) << '\n';
            errflag_ = IPX_ERROR_cr_iter_limit;
            break;
        }
        if (cdot <= 0.0) {
            control_.Debug(3)
                << " matrix in PCR method not posdef. cdot = " << sci2(cdot)
                << ", infnorm(sresidual) = " << sci2(Infnorm(sresidual))
                << ", infnorm(residual) = "  << sci2(Infnorm(residual))
                << '\n';
            errflag_ = IPX_ERROR_cr_matrix_not_posdef;
            break;
        }

        // Update lhs, residual, sresidual and Csresidual.
        double cdotnew;
        {
            // Uses Csresidual as storage for preconditioned Cstep.
            Vector& precond_Cstep = Csresidual;
            double pdot;
            P.Apply(Cstep, precond_Cstep, &pdot);
            if (pdot <= 0.0) {
                errflag_ = IPX_ERROR_cr_precond_not_posdef;
                break;
            }
            const double alpha = cdot/pdot;
            if (!std::isfinite(alpha)) {
                errflag_ = IPX_ERROR_cr_inf_or_nan;
                break;
            }
            lhs += alpha*step;
            residual -= alpha*Cstep;
            sresidual -= alpha*precond_Cstep;
            C.Apply(sresidual, Csresidual, &cdotnew);
            // Now Csresidual is restored and alias goes out of scope.
        }

        // Update step and Cstep.
        const double beta = cdotnew/cdot;
        step = sresidual + beta*step;
        Cstep = Csresidual + beta*Cstep;
        cdot = cdotnew;

        iter_++;
        if (iter_%5 == 0) {
            // We recompute the preconditioned residual from its definition
            // all 5 iterations. Otherwise it happened that the preconditioned
            // residual approached zero but the true residual stagnated. This
            // was caused by the update to the preconditioned residual being not
            // sufficiently accurate.
            // As a second safeguard, we check that resnorm_precond_system
            // decreased during the last 5 iterations.
            double rsdot;
            P.Apply(residual, sresidual, &rsdot);
            if (rsdot >= resnorm_precond_system) {
                control_.Debug(3)
                    << " resnorm_precond_system old = "
                    << sci2(resnorm_precond_system) << '\n'
                    << " resnorm_precond_system new = "
                    << sci2(rsdot) << '\n';
                errflag_ = IPX_ERROR_cr_no_progress;
                break;
            }
            resnorm_precond_system = rsdot;
        }

        if ((errflag_ = control_.InterruptCheck()) != 0)
            break;
    }
    time_ = timer.Elapsed();
}

Int ConjugateResiduals::errflag() const { return errflag_; }
Int ConjugateResiduals::iter() const { return iter_; }
double ConjugateResiduals::time() const { return time_; }

}  // namespace ipx

// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "ipm.h"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <limits>
#include "timer.h"
#include "utils.h"

namespace ipx {

struct IPM::Step {
    Step(Int m, Int n) : x(n+m), xl(n+m), xu(n+m), y(m), zl(n+m), zu(n+m) {}
    Vector x, xl, xu, y, zl, zu;
    Step& operator+=(const Step& rhs) {
        x += rhs.x; xl += rhs.xl; xu += rhs.xu;
        y += rhs.y; zl += rhs.zl; zu += rhs.zu;
        return *this;
    }
};

IPM::IPM(const Control& control) : control_(control) {}

void IPM::StartingPoint(KKTSolver* kkt, Iterate* iterate, Info* info) {
    kkt_ = kkt;
    iterate_ = iterate;
    info_ = info;
    PrintHeader();
    ComputeStartingPoint();
    if (info->errflag == 0)
        PrintOutput();
    // Set status_ipm.
    if (info->errflag == IPX_ERROR_interrupt_time) {
        info->errflag = 0;
        info->status_ipm = IPX_STATUS_time_limit;
    } else if (info->errflag) {
        info->status_ipm = IPX_STATUS_failed;
    } else {
        info->status_ipm = IPX_STATUS_not_run;
    }
}

void IPM::Driver(KKTSolver* kkt, Iterate* iterate, Info* info) {
    const Model& model = iterate->model();
    const Int m = model.rows();
    const Int n = model.cols();
    Step step(m, n);
    kkt_ = kkt;
    iterate_ = iterate;
    info_ = info;
    num_bad_iter_ = 0;

    while (true) {
        if (iterate->term_crit_reached()) {
            info->status_ipm = IPX_STATUS_optimal;
            break;
        }
        if (info->iter >= maxiter_) {
            info->status_ipm = IPX_STATUS_iter_limit;
            break;
        }
        if (num_bad_iter_ >= 5) {
            info->status_ipm = IPX_STATUS_no_progress;
            break;
        }
        if ((info->errflag = control_.InterruptCheck()) != 0)
            break;
        kkt->Factorize(iterate, info);
        if (info->errflag)
            break;
        Predictor(step);
        if (info->errflag)
            break;
        AddCorrector(step);
        if (info->errflag)
            break;
        MakeStep(step);
        info->iter++;
        PrintOutput();
    }

    // Set status_ipm if errflag terminated IPM.
    if (info->errflag) {
        if (info->errflag == IPX_ERROR_interrupt_time) {
            info->errflag = 0;
            info->status_ipm = IPX_STATUS_time_limit;
        } else {
            info->status_ipm = IPX_STATUS_failed;
        }
    }
}

void IPM::ComputeStartingPoint() {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const SparseMatrix& AI = model.AI();
    const Vector& b = model.b();
    const Vector& c = model.c();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    Vector x(n+m), xl(n+m), xu(n+m), y(m), zl(n+m), zu(n+m);
    Vector rb(m);               // workspace

    // Factorize the KKT matrix with the identity matrix in the (1,1) block.
    kkt_->Factorize(nullptr, info_);
    if (info_->errflag)
        return;

    // Set x within its bounds and compute the minimum norm solution dx to
    // AI*dx = (b-AI*x). Then update x := x + dx to obtain a feasible point.
    rb = b;
    for (Int j = 0; j < n+m; j++) {
        double xj = 0.0;
        if (xj < lb[j])
            xj = lb[j];
        if (xj > ub[j])
            xj = ub[j];
        x[j] = xj;
        if (xj != 0.0)
            ScatterColumn(AI, j, -xj, rb);
    }
    double tol = 0.1 * Infnorm(rb);
    zl = 0.0;
    kkt_->Solve(zl, rb, tol, xl, y, info_);
    if (info_->errflag)
        return;
    x += xl;

    // Compute xl, xu and shift to become positive.
    double xinfeas = 0.0;
    for (Int j = 0; j < n+m; j++) {
        xl[j] = x[j]-lb[j];
        xinfeas = std::max(xinfeas, -xl[j]);
        xu[j] = ub[j]-x[j];
        xinfeas = std::max(xinfeas, -xu[j]);
    }
    double xshift1 = 1.0 + 1.5*xinfeas;
    xl += xshift1;
    xu += xshift1;

    const double cnorm = Twonorm(c);
    if (cnorm == 0.0) {
        // Special treatment for zero objective.
        for (Int j = 0; j < n+m; j++) {
            zl[j] = std::isfinite(lb[j]) ? 1.0 : 0.0;
            zu[j] = std::isfinite(ub[j]) ? 1.0 : 0.0;
        }
    } else {
        // Compute y as the least-squares solution to AI'*y=c.
        // Recompute zl = c-AI'*y because the KKT system is solved
        // approximately with a residual in the first block equation.
        rb = 0.0;
        double tol = 0.1 * Infnorm(c);
        kkt_->Solve(c, rb, tol, zl, y, info_);
        if (info_->errflag)
            return;
        zl = c;
        MultiplyAdd(AI, y, -1.0, zl, 'T');

        // When c lies in range(AI'), then the dual slack variables are (close
        // to) zero, and the initial point would be almost complementary but
        // ususally not primal feasible. To prevent this from happening, add
        // a fraction of the objective to zl and adjust y. In exact computation
        // this does not affect dual feasiblity.
        const double znorm = Twonorm(zl);
        const double rho = 0.05;
        if (znorm < rho*cnorm) {
            zl += rho * c;
            y *= (1.0-rho);
        }

        // Split dual slack solution into zl, zu and shift to become positive.
        double zinfeas = 0.0;
        for (Int j = 0; j < n+m; j++) {
            double zval = zl[j];
            zl[j] = 0.0;
            zu[j] = 0.0;
            if (std::isfinite(lb[j]) && std::isfinite(ub[j])) {
                zl[j] = 0.5*zval;
                zu[j] = -0.5*zval;
            }
            else if (std::isfinite(lb[j]))
                zl[j] = zval;
            else if (std::isfinite(ub[j]))
                zu[j] = -zval;
            zinfeas = std::max(zinfeas, -zl[j]);
            zinfeas = std::max(zinfeas, -zu[j]);
        }
        double zshift1 = 1.0 + 1.5*zinfeas;
        for (Int j = 0; j < n+m; j++) {
            if (std::isfinite(lb[j]))
                zl[j] += zshift1;
            if (std::isfinite(ub[j]))
                zu[j] += zshift1;
        }
    }

    // Level pairwise complementarity products.
    double xsum = 1.0;
    double zsum = 1.0;
    double mu = 1.0;
    for (Int j = 0; j < n+m; j++) {
        if (std::isfinite(lb[j])) {
            xsum += xl[j];
            zsum += zl[j];
            mu += xl[j]*zl[j];
        }
        if (std::isfinite(ub[j])) {
            xsum += xu[j];
            zsum += zu[j];
            mu += xu[j]*zu[j];
        }
    }
    double xshift2 = 0.5*mu/zsum;
    double zshift2 = 0.5*mu/xsum;
    xl += xshift2;
    xu += xshift2;
    for (Int j = 0; j < n+m; j++) {
        if (std::isfinite(lb[j]))
            zl[j] += zshift2;
        if (std::isfinite(ub[j]))
            zu[j] += zshift2;
    }
    iterate_->Initialize(x, xl, xu, y, zl, zu);
}

// Computes maximum alpha such that x + alpha*dx >= 0.
// The blocking index is returned in blocking_index if not NULL.
static double StepToBoundary(const Vector& x, const Vector& dx,
                             Int* blocking_index, double alpha = 1.0) {
    const Int n = x.size();
    const double damp = 1.0 - std::numeric_limits<double>::epsilon();
    assert(damp < 1.0);

    Int iblock = -1;
    for (Int i = 0; i < n; i++) {
        assert(x[i] >= 0.0);
        if (x[i]+alpha*dx[i] < 0.0) {
            alpha = -(x[i]*damp) / dx[i];
            assert(x[i]+alpha*dx[i] >= 0.0);
            iblock = i;
        }
    }
    assert(alpha >= 0.0);
    if (blocking_index)
        *blocking_index = iblock;
    return alpha;
}

void IPM::Predictor(Step& step) {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& xl = iterate_->xl();
    const Vector& xu = iterate_->xu();
    const Vector& zl = iterate_->zl();
    const Vector& zu = iterate_->zu();

    // sl = -xl.*zl
    Vector sl(n+m);
    for (Int j = 0; j < n+m; j++)
        if (iterate_->has_barrier_lb(j))
            sl[j] = -xl[j]*zl[j];
        else
            sl[j] = 0.0;
    assert(AllFinite(sl));

    // su = -xu.*zu
    Vector su(n+m);
    for (Int j = 0; j < n+m; j++)
        if (iterate_->has_barrier_ub(j))
            su[j] = -xu[j]*zu[j];
        else
            su[j] = 0.0;
    assert(AllFinite(su));

    SolveNewtonSystem(&iterate_->rb()[0], &iterate_->rc()[0],
                      &iterate_->rl()[0], &iterate_->ru()[0], &sl[0], &su[0],
                      step);
}

void IPM::AddCorrector(Step& step) {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& xl = iterate_->xl();
    const Vector& xu = iterate_->xu();
    const Vector& zl = iterate_->zl();
    const Vector& zu = iterate_->zu();
    const Vector& dxl = step.xl;
    const Vector& dxu = step.xu;
    const Vector& dzl = step.zl;
    const Vector& dzu = step.zu;
    const double mu = iterate_->mu();

    // Choose centering parameter.
    double step_xl = StepToBoundary(xl, dxl, nullptr);
    double step_xu = StepToBoundary(xu, dxu, nullptr);
    double step_zl = StepToBoundary(zl, dzl, nullptr);
    double step_zu = StepToBoundary(zu, dzu, nullptr);
    double maxp = std::min(step_xl, step_xu);
    double maxd = std::min(step_zl, step_zu);
    double muaff = 0.0;
    Int num_finite = 0;
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->has_barrier_lb(j)) {
            assert(std::isfinite(xl[j]));
            assert(xl[j] != 0.0);
            muaff += (xl[j]+maxp*dxl[j]) * (zl[j]+maxd*dzl[j]);
            num_finite++;
        }
        if (iterate_->has_barrier_ub(j)) {
            assert(std::isfinite(xu[j]));
            assert(xu[j] != 0.0);
            muaff += (xu[j]+maxp*dxu[j]) * (zu[j]+maxd*dzu[j]);
            num_finite++;
        }
    }
    assert(std::isfinite(muaff));
    muaff /= num_finite;
    double ratio = muaff / mu;
    double sigma = ratio * ratio * ratio;

    // sl = -xl.*zl + sigma*mu - dxl.*dzl
    Vector sl(n+m);
    for (Int j = 0; j < n+m; j++)
        if (iterate_->has_barrier_lb(j))
            sl[j] = -xl[j]*zl[j] + sigma*mu - dxl[j]*dzl[j];
        else
            sl[j] = 0.0;
    assert(AllFinite(sl));

    // su = -xu.*zu + sigma*mu - dxu.*dzu
    Vector su(n+m);
    for (Int j = 0; j < n+m; j++)
        if (iterate_->has_barrier_ub(j))
            su[j] = -xu[j]*zu[j] + sigma*mu - dxu[j]*dzu[j];
        else
            su[j] = 0.0;
    assert(AllFinite(su));

    SolveNewtonSystem(&iterate_->rb()[0], &iterate_->rc()[0],
                      &iterate_->rl()[0], &iterate_->ru()[0], &sl[0], &su[0],
                      step);
}

void IPM::StepSizes(const Step& step) {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& xl = iterate_->xl();
    const Vector& xu = iterate_->xu();
    const Vector& zl = iterate_->zl();
    const Vector& zu = iterate_->zu();
    const Vector& dxl = step.xl;
    const Vector& dxu = step.xu;
    const Vector& dzl = step.zl;
    const Vector& dzu = step.zu;
    const double mu = iterate_->mu();
    const double gammaf = 0.9;
    const double gammaa = 1.0 / (1.0-gammaf);

    Int block_xl, block_xu, block_zl, block_zu;
    double step_xl = StepToBoundary(xl, dxl, &block_xl);
    double step_xu = StepToBoundary(xu, dxu, &block_xu);
    double step_zl = StepToBoundary(zl, dzl, &block_zl);
    double step_zu = StepToBoundary(zu, dzu, &block_zu);
    double maxp = std::fmin(step_xl, step_xu);
    double maxd = std::fmin(step_zl, step_zu);
    double mufull = 0.0;
    Int num_finite = 0;
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->has_barrier_lb(j)) {
            assert(std::isfinite(xl[j]));
            assert(xl[j] != 0.0);
            mufull += (xl[j]+maxp*dxl[j]) * (zl[j]+maxd*dzl[j]);
            num_finite++;
        }
        if (iterate_->has_barrier_ub(j)) {
            assert(std::isfinite(xu[j]));
            assert(xu[j] != 0.0);
            mufull += (xu[j]+maxp*dxu[j]) * (zu[j]+maxd*dzu[j]);
            num_finite++;
        }
    }
    assert(std::isfinite(mufull));
    mufull /= num_finite;
    mufull /= gammaa;

    double alphap = 1.0;
    double alphad = 1.0;
    Int blockp = -1;
    Int blockd = -1;
    if (maxp < 1.0) {
        if (step_xl <= step_xu) {
            double buffer;
            blockp = block_xl;
            buffer = mufull / (zl[blockp] + maxd*dzl[blockp]);
            alphap = (xl[blockp]-buffer) / (-dxl[blockp]);
        } else {
            double buffer;
            blockp = block_xu;
            buffer = mufull / (zu[blockp] + maxd*dzu[blockp]);
            alphap = (xu[blockp]-buffer) / (-dxu[blockp]);
        }
        alphap = std::max(alphap, gammaf*maxp);
        alphap = std::min(alphap, 1.0);
        assert(blockp >= 0.0);
    }
    if (maxd < 1.0) {
        if (step_zl <= step_zu) {
            double buffer;
            blockd = block_zl;
            buffer = mufull / (xl[blockd] + maxp*dxl[blockd]);
            alphad = (zl[blockd]-buffer) / (-dzl[blockd]);
        } else {
            double buffer;
            blockd = block_zu;
            buffer = mufull / (xu[blockd] + maxp*dxu[blockd]);
            alphad = (zu[blockd]-buffer) / (-dzu[blockd]);
        }
        alphad = std::max(alphad, gammaf*maxd);
        alphad = std::min(alphad, 1.0);
        assert(blockd >= 0.0);
    }
    step_primal_ = std::min(alphap, 1.0-1e-6);
    step_dual_   = std::min(alphad, 1.0-1e-6);
    (void)(mu);
}

void IPM::MakeStep(const Step& step) {
    StepSizes(step);
    iterate_->Update(step_primal_, &step.x[0], &step.xl[0], &step.xu[0],
                     step_dual_,   &step.y[0], &step.zl[0], &step.zu[0]);
    if (std::min(step_primal_, step_dual_) < 0.05)
        num_bad_iter_++;
    else
        num_bad_iter_ = 0;
}

void IPM::SolveNewtonSystem(const double* rb, const double* rc,
                                    const double* rl, const double* ru,
                                    const double* sl, const double* su,
                                    Step& step) {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const SparseMatrix& AI = model.AI();
    const Vector& xl = iterate_->xl();
    const Vector& xu = iterate_->xu();
    const Vector& zl = iterate_->zl();
    const Vector& zu = iterate_->zu();
    Vector& dx = step.x;
    Vector& dxl = step.xl;
    Vector& dxu = step.xu;
    Vector& dy = step.y;
    Vector& dzl = step.zl;
    Vector& dzu = step.zu;

    // Build RHS for KKT system.
    Vector rhs1(n+m), rhs2(m);
    if (rc) {
        for (Int j = 0; j < n+m; j++)
            rhs1[j] = -rc[j];
    }
    for (Int j = 0; j < n+m; j++) {
        double rlj = rl ? rl[j] : 0.0;
        double ruj = ru ? ru[j] : 0.0;
        if (iterate_->has_barrier_lb(j))
            rhs1[j] += (sl[j] + zl[j]*rlj) / xl[j];
        if (iterate_->has_barrier_ub(j))
            rhs1[j] -= (su[j] - zu[j]*ruj) / xu[j];
        if (iterate_->StateOf(j) == Iterate::State::fixed)
            rhs1[j] = 0.0;
    }
    assert(AllFinite(rhs1));
    if (rb)
        std::copy(rb, rb+m, std::begin(rhs2));

    // Solve KKT system.
    double tol =  control_.kkt_tol() * std::sqrt(iterate_->mu());
    kkt_->Solve(rhs1, rhs2, tol, dx, dy, info_);
    if (info_->errflag)
        return;

    // Recover solution to Newton system.
    dy *= -1.0;
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->StateOf(j) == Iterate::State::fixed) {
            assert(dx[j] == 0.0);
            dxl[j] = 0.0;
            dzl[j] = 0.0;
        } else if (iterate_->StateOf(j) == Iterate::State::free) {
            assert(!rl || rl[j] == 0.0);
            dxl[j] = 0.0;
            dzl[j] = 0.0;
        } else {
            double rlj = rl ? rl[j] : 0.0;
            dxl[j] = dx[j] - rlj;
            dzl[j] = (sl[j] - zl[j]*dxl[j]) / xl[j];
        }
    }
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->StateOf(j) == Iterate::State::fixed) {
            assert(dx[j] == 0.0);
            dxu[j] = 0.0;
            dzu[j] = 0.0;
        } else if (iterate_->StateOf(j) == Iterate::State::free) {
            assert(!ru || ru[j] == 0.0);
            dxu[j] = 0.0;
            dzu[j] = 0.0;
        } else {
            double ruj = ru ? ru[j] : 0.0;
            dxu[j] = ruj - dx[j];
            dzu[j] = (su[j] - zu[j]*dxu[j]) / xu[j];
        }
    }
    assert(AllFinite(dxl));
    assert(AllFinite(dxu));
    assert(AllFinite(dzl));
    assert(AllFinite(dzu));

    // Shift residual to the last two block equations.
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->StateOf(j) == Iterate::State::barrier) {
            assert(std::isfinite(xl[j]) || std::isfinite(xu[j]));
            double atdy = DotColumn(AI, j, dy);
            double rcj = rc ? rc[j] : 0.0;
            if (std::isfinite(xl[j]) && std::isfinite(xu[j])) {
                if (zl[j]*xu[j] >= zu[j]*xl[j])
                    dzl[j] = rcj + dzu[j] - atdy;
                else
                    dzu[j] = -rcj + dzl[j] + atdy;
            } else if (std::isfinite(xl[j])) {
                dzl[j] = rcj + dzu[j] - atdy;
            } else {
                dzu[j] = -rcj + dzl[j] + atdy;
            }
        }
    }

    #ifndef NDEBUG
    // Check solution for free and fixed variables.
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->StateOf(j) == Iterate::State::fixed) {
            assert(dx[j] == 0.0);
            assert(dxl[j] == 0.0 && dxu[j] == 0.0);
            assert(dzl[j] == 0.0 && dzu[j] == 0.0);
        }
        if (iterate_->StateOf(j) == Iterate::State::free)
            assert(dzl[j] == 0.0 && dzu[j] == 0.0);
    }
    #endif
}

void IPM::PrintHeader() {
    control_.Log()
        << " "  << Format("Iter", 4)
        << "  " << Format("P.res", 8) << " " << Format("D.res", 8)
        << "  " << Format("P.obj", 15) << " " << Format("D.obj", 15)
        << "  " << Format("mu", 8)
        << "  " << Format("Time", 7);
    control_.Debug()
        << "  " << Format("stepsizes", 9)
        << "  " << Format("pivots", 7) << " " << Format("kktiter", 7)
        << "  " << Format("P.fixed", 7) << " " << Format("D.fixed", 7);
    control_.Debug(4) << "  " << Format("svdmin(B)", 9);
    control_.Debug(4) << "  " << Format("density", 8);
    control_.Log() << '\n';
}

void IPM::PrintOutput() {
    const bool ipm_optimal = iterate_->feasible() && iterate_->optimal();

    control_.Log()
        << " "  << Format(info_->iter, 3)
        << (ipm_optimal ? "*" : " ")
        << "  " << Scientific(iterate_->presidual(), 8, 2)
        << " "  << Scientific(iterate_->dresidual(), 8, 2)
        << "  " << Scientific(iterate_->pobjective_after_postproc(), 15, 8)
        << " "  << Scientific(iterate_->dobjective_after_postproc(), 15, 8)
        << "  " << Scientific(iterate_->mu(), 8, 2)
        << "  " << Fixed(control_.Elapsed(), 6, 0) << "s";
    control_.Debug()
        << "  " << Fixed(step_primal_, 4, 2) << " " << Fixed(step_dual_, 4, 2)
        << "  " << Format(kkt_->basis_changes(), 7)
        << " "  << Format(kkt_->iter(), 7);
    control_.Debug()
        << "  " << Format(info_->dual_dropped, 7)
        << " "  << Format(info_->primal_dropped, 7);

    const Basis* basis = kkt_->basis();
    if (basis) {
        if (control_.Debug(4)) {
            control_.Debug(4) << "  "
                              << Scientific(basis->MinSingularValue(), 9, 2);
            Timer timer;
            double density = basis->DensityInverse();
            info_->time_symb_invert += timer.Elapsed();
            control_.Debug(4) << "  " << Scientific(density, 8, 2);
        }
    } else {
        control_.Debug(4) << "  " << Format("-", 9);
        control_.Debug(4) << "  " << Format("-", 8);
    }
    control_.Log() << '\n';
}

}  // namespace ipx

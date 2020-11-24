// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "crossover.h"
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <valarray>
#include "time.h"
#include "utils.h"

namespace ipx {

Crossover::Crossover(const Control& control) : control_(control) {}

void Crossover::PushAll(Basis* basis, Vector& x, Vector& y, Vector& z,
                        const double* weights, Info* info) {
    const Model& model = basis->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    std::vector<Int> perm = Sortperm(n+m, weights, false);

    control_.Log()
        << Textline("Primal residual before push phase:")
        << sci2(PrimalResidual(model, x)) << '\n'
        << Textline("Dual residual before push phase:")
        << sci2(DualResidual(model, y, z)) << '\n';

    // Run dual push phase.
    std::vector<Int> dual_superbasics;
    for (Int p = 0; p < (Int) perm.size(); p++) {
        Int j = perm[p];
        if (basis->IsBasic(j) && z[j] != 0.0)
            dual_superbasics.push_back(j);
    }
    control_.Log()
        << Textline("Number of dual pushes required:")
        << dual_superbasics.size() << '\n';
    PushDual(basis, y, z, dual_superbasics, x, info);
    assert(DualInfeasibility(model, x, z) == 0.0);
    if (info->status_crossover != IPX_STATUS_optimal)
        return;

    // Run primal push phase. Because z[j]==0 for all basic variables, none of
    // the primal variables is fixed at its bound.
    std::vector<Int> primal_superbasics;
    for (Int p = perm.size()-1; p >= 0; p--) {
        Int j = perm[p];
        if (basis->IsNonbasic(j) && x[j] != lb[j] && x[j] != ub[j] &&
            !(std::isinf(lb[j]) && std::isinf(ub[j]) && x[j] == 0.0))
            primal_superbasics.push_back(j);
    }
    control_.Log()
        << Textline("Number of primal pushes required:")
        << primal_superbasics.size() << '\n';
    PushPrimal(basis, x, primal_superbasics, nullptr, info);
    assert(PrimalInfeasibility(model, x) == 0.0);
    if (info->status_crossover != IPX_STATUS_optimal)
        return;

    control_.Debug()
        << Textline("Primal residual after push phase:")
        << sci2(PrimalResidual(model, x)) << '\n'
        << Textline("Dual residual after push phase:")
        << sci2(DualResidual(model, y, z)) << '\n';
    info->status_crossover = IPX_STATUS_optimal;
}

void Crossover::PushPrimal(Basis* basis, Vector& x,
                           const std::vector<Int>& variables,
                           const bool* fixed_at_bound, Info* info) {
    Timer timer;
    const Model& model = basis->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    IndexedVector ftran(m);
    const double feastol = model.dualized() ?
        control_.dfeasibility_tol() : control_.pfeasibility_tol();
    primal_pushes_ = 0;
    primal_pivots_ = 0;

    // Check that variables are nonbasic and that x satisfies bound condition.
    for (Int j : variables) {
        if (!basis->IsNonbasic(j))
            throw std::logic_error("invalid variable in Crossover::PushPrimal");
    }
    for (Int j = 0; j < n+m; j++) {
        if (x[j] < lb[j] || x[j] > ub[j])
            throw std::logic_error(
                "bound condition violated in Crossover::PushPrimal");
        if (fixed_at_bound && fixed_at_bound[j] &&
            x[j] != lb[j] && x[j] != ub[j])
            throw std::logic_error(
                "bound condition violated in Crossover::PushPrimal");
    }

    // Maintain a copy of primal basic variables and their bounds for faster
    // ratio test. Fixed-at-bound variables are handled by setting their bounds
    // equal.
    Vector xbasic  = CopyBasic(x,  *basis);
    Vector lbbasic = CopyBasic(lb, *basis);
    Vector ubbasic = CopyBasic(ub, *basis);
    if (fixed_at_bound) {
        for (Int p = 0; p < m; p++) {
            Int j = (*basis)[p];
            if (fixed_at_bound[j])
                lbbasic[p] = ubbasic[p] = x[j];
        }
    }

    control_.ResetPrintInterval();
    Int next = 0;
    while (next < (Int) variables.size()) {
        if ((info->errflag = control_.InterruptCheck()) != 0)
            break;

        const Int jn = variables[next];
        if (x[jn] == lb[jn] || x[jn] == ub[jn] ||
            (x[jn] == 0.0 && std::isinf(lb[jn]) && std::isinf(ub[jn]))) {
            // nothing to do
            next++;
            continue;
        }
        // Choose bound to push to. If the variable has two finite bounds, move
        // to the nearer. If it has none, move to zero.
        double move_to = 0.0;
        if (std::isfinite(lb[jn]) && std::isfinite(ub[jn]))
            move_to = x[jn]-lb[jn] <= ub[jn]-x[jn] ? lb[jn] : ub[jn];
        else if (std::isfinite(lb[jn]))
            move_to = lb[jn];
        else if (std::isfinite(ub[jn]))
            move_to = ub[jn];

        // A full step is such that x[jn]-step is at its bound.
        double step = x[jn]-move_to;

        basis->SolveForUpdate(jn, ftran);
        bool block_at_lb;
        Int pblock = PrimalRatioTest(xbasic, ftran, lbbasic, ubbasic,
                                     step, feastol, &block_at_lb);
        Int jb = pblock >= 0 ? (*basis)[pblock] : -1;

        // If step was blocked, update basis and compute step size.
        if (pblock >= 0) {
            double pivot = ftran[pblock];
            assert(pivot != 0.0);
            if (std::abs(pivot) < 1e-4)
                control_.Debug(3)
                    << " |pivot| = " << sci2(std::abs(pivot)) << '\n';
            bool exchanged;
            info->errflag = basis->ExchangeIfStable(jb, jn, pivot, -1,
                                                    &exchanged);
            if (info->errflag) {
                control_.Debug()
                    << Textline("Minimum singular value of basis matrix:")
                    << sci2(basis->MinSingularValue()) << '\n';
                break;
            }
            if (!exchanged)     // factorization was unstable, try again
                continue;
            primal_pivots_++;
            // We must use lbbasic[pblock] and ubbasic[pblock] (and not lb[jb]
            // and ub[jb]) so that step is 0.0 if a fixed-at-bound variable
            // blocked.
            if (block_at_lb)
                step = (lbbasic[pblock]-xbasic[pblock]) / ftran[pblock];
            else
                step = (ubbasic[pblock]-xbasic[pblock]) / ftran[pblock];
        }
        // Update solution.
        if (step != 0.0) {
            auto update = [&](Int p, double pivot) {
                xbasic[p] += step * pivot;
                xbasic[p] = std::max(xbasic[p], lbbasic[p]);
                xbasic[p] = std::min(xbasic[p], ubbasic[p]);
            };
            for_each_nonzero(ftran, update);
            x[jn] -= step;
        }
        if (pblock >= 0) {
            // make clean
            x[jb] = block_at_lb ? lbbasic[pblock] : ubbasic[pblock];
            assert(std::isfinite(x[jb]));
            // Update copy of basic variables and bounds. Note: jn cannot be
            // a fixed-at-bound variable since it was pushed to a bound.
            xbasic[pblock] = x[jn];
            lbbasic[pblock] = lb[jn];
            ubbasic[pblock] = ub[jn];
        } else {
            x[jn] = move_to; // make clean
            assert(std::isfinite(x[jn]));
        }

        primal_pushes_++;
        next++;
        control_.IntervalLog()
            << " " << Format(static_cast<Int>(variables.size()-next), 8)
            << " primal pushes remaining"
            << " (" << Format(primal_pivots_, 7) << " pivots)\n";
    }
    for (Int p = 0; p < m; p++)
        x[(*basis)[p]] = xbasic[p];

    // Set status flag.
    if (info->errflag == IPX_ERROR_interrupt_time) {
        info->errflag = 0;
        info->status_crossover = IPX_STATUS_time_limit;
    } else if (info->errflag != 0) {
        info->status_crossover = IPX_STATUS_failed;
    } else {
        info->status_crossover = IPX_STATUS_optimal;
    }
    time_primal_ = timer.Elapsed();
}

void Crossover::PushPrimal(Basis* basis, Vector& x,
                           const std::vector<Int>& variables,
                           const Vector& z, Info* info) {
    std::valarray<bool> bound_restrict = z != 0.0;
    PushPrimal(basis, x, variables, &bound_restrict[0], info);
}

void Crossover::PushDual(Basis* basis, Vector& y, Vector& z,
                         const std::vector<Int>& variables,
                         const int sign_restrict[], Info* info) {
    Timer timer;
    const Model& model = basis->model();
    const Int m = model.rows();
    const Int n = model.cols();
    IndexedVector btran(m), row(n+m);
    const double feastol = model.dualized() ?
        control_.pfeasibility_tol() : control_.dfeasibility_tol();
    dual_pushes_ = 0;
    dual_pivots_ = 0;

    // Check that variables are basic and that z satisfies sign condition.
    for (Int j : variables) {
        if (!basis->IsBasic(j))
            throw std::logic_error("invalid variable in Crossover::PushDual");
    }
    for (Int j = 0; j < n+m; j++) {
        if (((sign_restrict[j] & 1) && z[j] < 0.0) ||
            ((sign_restrict[j] & 2) && z[j] > 0.0))
            throw std::logic_error(
                "sign condition violated in Crossover::PushDual");
    }

    control_.ResetPrintInterval();
    Int next = 0;
    while (next < (Int) variables.size()) {
        if ((info->errflag = control_.InterruptCheck()) != 0)
            break;

        const Int jb = variables[next];
        if (z[jb] == 0.0) {
            // nothing to do
            next++;
            continue;
        }
        // The update operation applied below is
        // y := y + step*btran, z := z - step*row, z[jb] := z[jb] - step,
        // where row is the tableau row for variable jb. In exact arithmetic
        // this leaves A'y+z unchanged.
        basis->TableauRow(jb, btran, row);
        double step = z[jb];
        Int jn = DualRatioTest(z, row, sign_restrict, step, feastol);

        // If step was blocked, update basis and compute step size.
        if (jn >= 0) {
            assert(basis->IsNonbasic(jn));
            double pivot = row[jn];
            assert(pivot);
            if (std::abs(pivot) < 1e-4)
                control_.Debug(3)
                    << " |pivot| = " << sci2(std::abs(pivot)) << '\n';
            bool exchanged;
            info->errflag = basis->ExchangeIfStable(jb, jn, pivot, 1,
                                                    &exchanged);
            if (info->errflag) {
                control_.Debug()
                    << Textline("Minimum singular value of basis matrix:")
                    << sci2(basis->MinSingularValue()) << '\n';
                break;
            }
            if (!exchanged)     // factorization was unstable, try again
                continue;
            dual_pivots_++;
            step = z[jn]/row[jn];
            // Update must move z[jb] toward zero.
            if (sign_restrict[jb] & 1)
                assert(step >= 0.0);
            if (sign_restrict[jb] & 2)
                assert(step <= 0.0);
        }
        // Update solution.
        if (step != 0.0) {
            auto update_y = [&](Int i, double x) {
                y[i] += step*x;
            };
            for_each_nonzero(btran, update_y);
            auto update_z = [&](Int j, double pivot) {
                z[j] -= step * pivot;
                if (sign_restrict[j] & 1)
                    z[j] = std::max(z[j], 0.0);
                if (sign_restrict[j] & 2)
                    z[j] = std::min(z[j], 0.0);
            };
            for_each_nonzero(row, update_z);
            z[jb] -= step;
        }
        if (jn >= 0)
            z[jn] = 0.0; // make clean
        else
            assert(z[jb] == 0.0);

        dual_pushes_++;
        next++;
        control_.IntervalLog()
            << " " << Format(static_cast<Int>(variables.size()-next), 8)
            << " dual pushes remaining"
            << " (" << Format(dual_pivots_, 7) << " pivots)\n";
    }

    // Set status flag.
    if (info->errflag == IPX_ERROR_interrupt_time) {
        info->errflag = 0;
        info->status_crossover = IPX_STATUS_time_limit;
    } else if (info->errflag != 0) {
        info->status_crossover = IPX_STATUS_failed;
    } else {
        info->status_crossover = IPX_STATUS_optimal;
    }
    time_dual_ = timer.Elapsed();
}

void Crossover::PushDual(Basis* basis, Vector& y, Vector& z,
                         const std::vector<Int>& variables,
                         const Vector& x, Info* info) {
    const Model& model = basis->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();

    std::vector<int> sign_restrict(n+m);
    for (Int j = 0; j < (Int) sign_restrict.size(); j++) {
        if (x[j] != ub[j]) sign_restrict[j] |= 1;
        if (x[j] != lb[j]) sign_restrict[j] |= 2;
    }
    PushDual(basis, y, z, variables, sign_restrict.data(), info);
}

Int Crossover::PrimalRatioTest(const Vector& xbasic, const IndexedVector& ftran,
                               const Vector& lbbasic, const Vector& ubbasic,
                               double step, double feastol, bool* block_at_lb) {
    Int pblock = -1;            // return value
    *block_at_lb = true;

    // First pass: determine maximum step size exploiting feasibility tol.
    auto update_step = [&](Int p, double pivot) {
        if (std::abs(pivot) > kPivotZeroTol) {
            // test block at lower bound
            if (xbasic[p] + step*pivot < lbbasic[p]-feastol) {
                step = (lbbasic[p]-xbasic[p]-feastol) / pivot;
                pblock = p;
                *block_at_lb = true;
            }
            // test block at upper bound
            if (xbasic[p] + step*pivot > ubbasic[p]+feastol) {
                step = (ubbasic[p]-xbasic[p]+feastol) / pivot;
                pblock = p;
                *block_at_lb = false;
            }
        }
    };
    for_each_nonzero(ftran, update_step);

    // If the step was not blocked, we are done.
    if (pblock < 0)
        return pblock;

    // Second pass: choose maximum pivot among all that block within step.
    pblock = -1;
    double max_pivot = kPivotZeroTol;
    auto update_max = [&](Int p, double pivot) {
        if (std::abs(pivot) > max_pivot) {
            // test block at lower bound
            if (step*pivot < 0.0) {
                double step_p = (lbbasic[p]-xbasic[p]) / pivot;
                if (std::abs(step_p) <= std::abs(step)) {
                    pblock = p;
                    *block_at_lb = true;
                    max_pivot = std::abs(pivot);
                }
            }
            // test block at upper bound
            if (step*pivot > 0.0) {
                double step_p = (ubbasic[p]-xbasic[p]) / pivot;
                if (std::abs(step_p) <= std::abs(step)) {
                    pblock = p;
                    *block_at_lb = false;
                    max_pivot = std::abs(pivot);
                }
            }
        }
    };
    for_each_nonzero(ftran, update_max);
    assert(pblock >= 0);
    return pblock;
}

Int Crossover::DualRatioTest(const Vector& z, const IndexedVector& row,
                             const int sign_restrict[], double step,
                             double feastol) {
    Int jblock = -1;            // return value

    // First pass: determine maximum step size exploiting feasibility tol.
    auto update_step = [&](Int j, double pivot) {
        if (std::abs(pivot) > kPivotZeroTol) {
            if ((sign_restrict[j] & 1) && z[j]-step*pivot < -feastol) {
                step = (z[j]+feastol) / pivot;
                jblock = j;
                assert(z[j] >= 0.0);
                assert(step*pivot > 0.0);
            }
            if ((sign_restrict[j] & 2) && z[j]-step*pivot > feastol) {
                step = (z[j]-feastol) / pivot;
                jblock = j;
                assert(z[j] <= 0.0);
                assert(step*pivot < 0.0);
            }
        }
    };
    for_each_nonzero(row, update_step);

    // If step was not block, we are done.
    if (jblock < 0)
        return jblock;

    // Second pass: choose maximum pivot among all that block within step.
    jblock = -1;
    double max_pivot = kPivotZeroTol;
    auto update_max = [&](Int j, double pivot) {
        if (std::abs(pivot) > max_pivot &&
            std::abs(z[j]/pivot) <= std::abs(step)) {
            if ((sign_restrict[j] & 1) && step*pivot > 0.0) {
                jblock = j;
                max_pivot = std::abs(pivot);
            }
            if ((sign_restrict[j] & 2) && step*pivot < 0.0) {
                jblock = j;
                max_pivot = std::abs(pivot);
            }
        }
    };
    for_each_nonzero(row, update_max);
    assert(jblock >= 0);
    return jblock;
}

}  // namespace ipx

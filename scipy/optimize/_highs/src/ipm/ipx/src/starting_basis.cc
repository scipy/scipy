// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "starting_basis.h"
#include <cassert>
#include <cmath>
#include <vector>
#include "timer.h"

namespace ipx {

// Asserts that variables have a state and basic status that is consistent with
// their bounds in the model.
static void AssertConsistency(const Iterate& iterate, const Basis& basis) {
    const Model& model = basis.model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();

    for (Int j = 0; j < n+m; j++) {
        if (lb[j] == ub[j]) {
            switch (iterate.StateOf(j)) {
            case Iterate::State::fixed:
                assert(basis.StatusOf(j) == Basis::NONBASIC_FIXED);
                break;
            case Iterate::State::free:
                assert(j >= n); // must be slack variable
                assert(basis.StatusOf(j) == Basis::BASIC_FREE);
                break;
            default:
                assert(0);
            }
        } else if (std::isinf(lb[j]) && std::isinf(ub[j])) {
            switch (iterate.StateOf(j)) {
            case Iterate::State::fixed:
                assert(basis.StatusOf(j) == Basis::NONBASIC_FIXED);
                break;
            case Iterate::State::free:
                assert(basis.StatusOf(j) == Basis::BASIC_FREE);
                break;
            default:
                assert(0);
            }
        } else {
            assert(iterate.StateOf(j) == Iterate::State::barrier);
            assert(basis.StatusOf(j) == Basis::NONBASIC ||
                   basis.StatusOf(j) == Basis::BASIC);
        }
    }
}

static void PostprocessDependencies(Iterate* iterate, Basis* p_basis,
                                    Info* info) {
    const Model& model = iterate->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const SparseMatrix& AI = model.AI();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    const Vector& x = iterate->x();
    const Vector& y = iterate->y();
    const Vector& zl = iterate->zl();
    const Vector& zu = iterate->zu();
    Basis& basis = *p_basis;
    std::vector<Int> dependent_rows, dependent_cols;
    Vector dx(n+m), dy(m);

    // Postprocess linearly dependent columns. If a free variable could not be
    // made basic, then the variable is redundant and can be fixed at an
    // arbitrary value. We fix it at zero and adapt the basic variables to keep
    // AI*x unchanged.
    if (info->dependent_cols > 0) {
        Vector dxbasic(m);
        for (Int j = 0; j < n; j++) {
            if (std::isinf(lb[j]) && std::isinf(ub[j]) && basis.IsNonbasic(j)) {
                assert(iterate->StateOf(j) == Iterate::State::free);
                assert(basis.StatusOf(j) == Basis::NONBASIC_FIXED);
                dx[j] = -x[j];
                ScatterColumn(AI, j, x[j], dxbasic);
                dependent_cols.push_back(j);
            }
        }
        basis.SolveDense(dxbasic, dxbasic, 'N');
        for (Int p = 0; p < m; p++)
            dx[basis[p]] = dxbasic[p];
    }

    // Postprocess linearly dependent rows. If the slack variable to an equality
    // constraint cannot be made nonbasic, then the constraint is redundant. We
    // move its dual variable y[i] to zero without altering A'*y and make the
    // constraint "implied" (i.e. the slack variable is treated as free by the
    // interior point solver).
    if (info->dependent_rows > 0) {
        for (Int p = 0; p < m; p++) {
            Int j = basis[p];
            // The starting basis is constructed such that structural variables
            // with weight zero will always be nonbasic.
            if (lb[j] == ub[j])
                assert(j >= n);
            if (j >= n && lb[j] == ub[j]) {
                assert(basis.StatusOf(j) == Basis::BASIC_FREE);
                assert(iterate->has_barrier_lb(j));
                assert(iterate->has_barrier_ub(j));
                dy[p] = -y[j-n];
                dependent_rows.push_back(j-n);
            }
        }
        basis.SolveDense(dy, dy, 'T');
        for (Int i : dependent_rows)
            dy[i] = -y[i];      // would be already in exact arithmetic
    }

    // Now adjust the iterate.
    iterate->Update(1.0, &dx[0], nullptr, nullptr, 1.0, &dy[0], nullptr,
                    nullptr);

    for (Int j : dependent_cols) {
        assert(x[j] == 0.0);
        iterate->make_fixed(j, 0.0);
    }
    for (Int i : dependent_rows) {
        assert(y[i] == 0.0);
        iterate->make_implied_eq(n+i); // sets zl[n+i] = zu[n+i] = 0
        assert(zl[n+i] == 0.0);
        assert(zu[n+i] == 0.0);
    }
    (void)(zl);
    (void)(zu);
}

void StartingBasis(Iterate* iterate, Basis* p_basis, Info* info) {
    const Model& model = iterate->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    Vector colscale(n+m);
    Basis& basis = *p_basis;
    info->errflag = 0;
    Timer timer;

    // Construct starting basis. The column weights for the crash procedure are
    // the interior point scaling factors from the current iterate, except that
    // fixed variables get weight zero.
    for (Int j = 0; j < n+m; j++) {
        colscale[j] = iterate->ScalingFactor(j);
        if (std::isinf(lb[j]) && std::isinf(ub[j]))
            assert(colscale[j] == INFINITY);
        else
            assert(std::isfinite(colscale[j]));
        if (lb[j] == ub[j])
            colscale[j] = 0.0;
    }
    basis.ConstructBasisFromWeights(&colscale[0], info);
    if (info->errflag)
        return;

    // Change status of free variables either to BASIC_FREE or NONBASIC_FIXED.
    // Change status of fixed variables either to NONBASIC_FIXED or BASIC_FREE.
    for (Int j = 0; j < n+m; j++) {
        if (colscale[j] == 0.0 || std::isinf(colscale[j])) {
            if (basis.IsBasic(j))
                basis.FreeBasicVariable(j);
            else
                basis.FixNonbasicVariable(j);
        }
    }

    // Variables with equal lower and upper bounds can be removed from the
    // interior point solve if their basic status is NONBASIC_FIXED. Setting
    // x[j] to its bound alters the residual in AI*x=b; this can be ignored
    // because the iterate has not been primal feasible if x[j] was
    // significantly away from its bound.
    for (Int j = 0; j < n+m; j++) {
        if (lb[j] == ub[j] && basis.StatusOf(j) == Basis::NONBASIC_FIXED)
            iterate->make_fixed(j, lb[j]);
    }

    // Linearly dependent rows and columns require adjusting the iterate.
    PostprocessDependencies(iterate, p_basis, info);

    AssertConsistency(*iterate, basis);
    info->time_starting_basis += timer.Elapsed();
}

}  // namespace ipx

// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include <algorithm>
#include "iterate.h"
#include <cassert>
#include <cmath>
#include "utils.h"

namespace ipx {

// Iterate::kBarrierMin is odr-used because std::max() takes references as
// arguments. Hence we require a namespace scope definition.
constexpr double Iterate::kBarrierMin;

Iterate::Iterate(const Model& model) : model_(model) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    x_.resize(n+m);
    xl_.resize(n+m);
    xu_.resize(n+m);
    y_.resize(m);
    zl_.resize(n+m);
    zu_.resize(n+m);
    rb_.resize(m);
    rl_.resize(n+m);
    ru_.resize(n+m);
    rc_.resize(n+m);
    variable_state_.resize(n+m);

    for (Int j = 0; j < n+m; j++) {
        if (std::isfinite(lb[j]) && std::isfinite(ub[j])) {
            variable_state_[j] = StateDetail::BARRIER_BOXED;
            xl_[j] = 1.0;
            xu_[j] = 1.0;
            zl_[j] = 1.0;
            zu_[j] = 1.0;
        } else if (std::isfinite(lb[j])) {
            variable_state_[j] = StateDetail::BARRIER_LB;
            xl_[j] = 1.0;
            xu_[j] = INFINITY;
            zl_[j] = 1.0;
            zu_[j] = 0.0;
        } else if (std::isfinite(ub[j])) {
            variable_state_[j] = StateDetail::BARRIER_UB;
            xl_[j] = INFINITY;
            xu_[j] = 1.0;
            zl_[j] = 0.0;
            zu_[j] = 1.0;
        } else {
            variable_state_[j] = StateDetail::FREE;
            xl_[j] = INFINITY;
            xu_[j] = INFINITY;
            zl_[j] = 0.0;
            zu_[j] = 0.0;
        }
    }
    assert_consistency();
}

void Iterate::Initialize(const Vector& x, const Vector& xl, const Vector& xu,
                         const Vector& y, const Vector& zl, const Vector& zu) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    assert((int)x.size() == n+m);
    assert((int)xl.size() == n+m);
    assert((int)xu.size() == n+m);
    assert((int)y.size() == m);
    assert((int)zl.size() == n+m);
    assert((int)zu.size() == n+m);
    x_ = x; xl_ = xl; xu_ = xu; y_ = y; zl_ = zl; zu_ = zu;

    // Set variable statuses.
    for (Int j = 0; j < n+m; j++) {
        if (lb[j] == ub[j]) {
            variable_state_[j] = StateDetail::BARRIER_BOXED;
        } else if (std::isfinite(lb[j]) && std::isfinite(ub[j])) {
            variable_state_[j] = StateDetail::BARRIER_BOXED;
        } else if (std::isfinite(lb[j])) {
            variable_state_[j] = StateDetail::BARRIER_LB;
        } else if (std::isfinite(ub[j])) {
            variable_state_[j] = StateDetail::BARRIER_UB;
        } else {
            variable_state_[j] = StateDetail::FREE;
        }
    }
    assert_consistency();
    evaluated_ = false;
    postprocessed_ = false;
}

void Iterate::Update(double sp, const double* dx, const double* dxl,
                     const double* dxu, double sd, const double* dy,
                     const double* dzl, const double* dzu) {
    const Int m = model_.rows();
    const Int n = model_.cols();

    if (dx) {
        for (Int j = 0; j < n+m; j++)
            if (StateOf(j) != State::fixed)
                x_[j] += sp*dx[j];
    }
    if (dxl) {
        for (Int j = 0; j < n+m; j++)
            if (has_barrier_lb(j)) {
                xl_[j] += sp*dxl[j];
                xl_[j] = std::max(xl_[j], kBarrierMin);
            }
    }
    if (dxu) {
        for (Int j = 0; j < n+m; j++)
            if (has_barrier_ub(j)) {
                xu_[j] += sp*dxu[j];
                xu_[j] = std::max(xu_[j], kBarrierMin);
            }
    }
    if (dy) {
        for (Int i = 0; i < m; i++)
            y_[i] += sd*dy[i];
    }
    if (dzl) {
        for (Int j = 0; j < n+m; j++)
            if (has_barrier_lb(j)) {
                zl_[j] += sd*dzl[j];
                zl_[j] = std::max(zl_[j], kBarrierMin);
            }
    }
    if (dzu) {
        for (Int j = 0; j < n+m; j++)
            if (has_barrier_ub(j)) {
                zu_[j] += sd*dzu[j];
                zu_[j] = std::max(zu_[j], kBarrierMin);
            }
    }
    assert_consistency();
    evaluated_ = false;
}

const Vector& Iterate::rb() const { Evaluate(); return rb_; }
const Vector& Iterate::rl() const { Evaluate(); return rl_; }
const Vector& Iterate::ru() const { Evaluate(); return ru_; }
const Vector& Iterate::rc() const { Evaluate(); return rc_; }

void Iterate::make_fixed(Int j) {
    xl_[j] = 0.0;
    xu_[j] = 0.0;
    zl_[j] = 0.0;
    zu_[j] = 0.0;
    variable_state_[j] = StateDetail::FIXED;
    evaluated_ = false;
}

void Iterate::make_fixed(Int j, double value) {
    x_[j] = value;
    make_fixed(j);
}

void Iterate::make_implied_lb(Int j) {
    xl_[j] = INFINITY;
    xu_[j] = INFINITY;
    variable_state_[j] = StateDetail::IMPLIED_LB;
    evaluated_ = false;
}

void Iterate::make_implied_ub(Int j) {
    xl_[j] = INFINITY;
    xu_[j] = INFINITY;
    variable_state_[j] = StateDetail::IMPLIED_UB;
    evaluated_ = false;
}

void Iterate::make_implied_eq(Int j) {
    xl_[j] = INFINITY;
    xu_[j] = INFINITY;
    zl_[j] = 0.0;
    zu_[j] = 0.0;
    variable_state_[j] = StateDetail::IMPLIED_EQ;
    evaluated_ = false;
}

double Iterate::ScalingFactor(Int j) const {
    switch (StateOf(j)) {
    case State::fixed:
        return 0.0;
    case State::free:
        return INFINITY;
    default:
        assert(xl_[j] > 0.0);
        assert(xu_[j] > 0.0);
        double g = zl_[j]/xl_[j] + zu_[j]/xu_[j];
        double d = 1.0 / std::sqrt(g);
        assert(std::isfinite(d));
        assert(d > 0.0);
        return d;
    }
}

double Iterate::pobjective() const { Evaluate(); return pobjective_; }
double Iterate::dobjective() const { Evaluate(); return dobjective_; }

double Iterate::pobjective_after_postproc() const {
    Evaluate();
    return pobjective_ + offset_;
}

double Iterate::dobjective_after_postproc() const {
    Evaluate();
    return dobjective_ + offset_;
}

double Iterate::presidual() const { Evaluate(); return presidual_; }
double Iterate::dresidual() const { Evaluate(); return dresidual_; }

double Iterate::complementarity() const { Evaluate(); return complementarity_; }
double Iterate::mu() const { Evaluate(); return mu_; }
double Iterate::mu_min() const { Evaluate(); return mu_min_; }
double Iterate::mu_max() const { Evaluate(); return mu_max_; }

bool Iterate::feasible() const {
    Evaluate();
    return
        presidual_ <= feasibility_tol_ * (1.0+model_.norm_bounds()) &&
        dresidual_ <= feasibility_tol_ * (1.0+model_.norm_c());
}

bool Iterate::optimal() const {
    Evaluate();
    double pobj = pobjective_after_postproc();
    double dobj = dobjective_after_postproc();
    double obj = 0.5 * (pobj + dobj);
    double gap = pobj - dobj;
    return std::abs(gap) <= optimality_tol_ * (1.0+std::abs(obj)); 
}

bool Iterate::term_crit_reached() const {
    if (feasible() && optimal()) {
        if (crossover_start_ <= 0.0)
            return true;
        double pres, dres;
        ResidualsFromDropping(&pres, &dres);
        if (pres <= crossover_start_ * (1.0+model_.norm_bounds()) &&
            dres <= crossover_start_ * (1.0+model_.norm_c()))
            return true;
    }
    return false;
}

void Iterate::Postprocess() {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& c = model_.c();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    const SparseMatrix& AI = model_.AI();

    // For fixed variables compute xl[j] and xu[j] from x[j]. If the lower and
    // upper bound are equal, set zl[j] or zu[j] such that the variable is dual
    // feasibile. Otherwise leave them zero.
    for (Int j = 0; j < n+m; j++) {
        if (StateOf(j) == State::fixed) {
            xl_[j] = x_[j] - lb[j];
            xu_[j] = ub[j] - x_[j];
            assert(zl_[j] == 0.0);
            assert(zu_[j] == 0.0);
            if (lb[j] == ub[j]) {
                double z = c[j] - DotColumn(AI, j, y_);
                if (z >= 0.0)
                    zl_[j] = z;
                else
                    zu_[j] = -z;
            }
        }
    }
    // For implied variables set x[j] to the bound at which it was implied and
    // compute zl[j] or zu[j]. If the variable was implied at both bounds,
    // choose between zl and zu depending on sign.
    for (Int j = 0; j < n+m; j++) {
        if (is_implied(j)) {
            double z = c[j] - DotColumn(AI, j, y_);
            switch (variable_state_[j]) {
            case StateDetail::IMPLIED_EQ:
                assert(lb[j] == ub[j]);
                if (z >= 0.0) {
                    zl_[j] = z;
                    zu_[j] = 0.0;
                } else {
                    zl_[j] = 0.0;
                    zu_[j] = -z;
                }
                x_[j] = lb[j];
                break;
            case StateDetail::IMPLIED_LB:
                zl_[j] = z;
                zu_[j] = 0.0;
                x_[j] = lb[j];
                break;
            case StateDetail::IMPLIED_UB:
                zl_[j] = 0.0;
                zu_[j] = -z;
                x_[j] = ub[j];
                break;
            default:
                assert(0);
            }
            xl_[j] = x_[j] - lb[j];
            xu_[j] = ub[j] - x_[j];
        }
    }
    postprocessed_ = true;
    evaluated_ = false;
}

void Iterate::EvaluatePostsolved(Info* info) const {
    model_.EvaluateInteriorSolution(x_, xl_, xu_, y_, zl_, zu_, info);
}

void Iterate::DropToComplementarity(Vector& x, Vector& y, Vector& z) const {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    assert(postprocessed_);

    assert((int)x.size() == n+m);
    assert((int)y.size() == m);
    assert((int)z.size() == n+m);
    y = y_;
    for (Int j = 0; j < n+m; j++) {
        const double xlj = xl_[j];
        const double xuj = xu_[j];
        const double zlj = zl_[j];
        const double zuj = zu_[j];
        assert(xlj >= 0.0);
        assert(xuj >= 0.0);
        assert(zlj >= 0.0);
        assert(zuj >= 0.0);
        double xj = x_[j];
        xj = std::max(xj, lb[j]);
        xj = std::min(xj, ub[j]);

        if (lb[j] == ub[j]) {
            // fixed variable
            x[j] = lb[j];
            z[j] = zlj-zuj;
        } else if (std::isfinite(lb[j]) && std::isfinite(ub[j])) {
            // boxed variable
            if (zlj*xuj >= zuj*xlj) {
                // either active at lower bound or inactive
                if (zlj >= xlj) {
                    x[j] = lb[j];
                    z[j] = std::max(0.0, zlj-zuj);
                } else {
                    x[j] = xj;
                    z[j] = 0.0;
                }
            } else {
                // either active at upper bound or inactive
                if (zuj >= xuj) {
                    x[j] = ub[j];
                    z[j] = std::min(0.0, zlj-zuj);
                } else {
                    x[j] = xj;
                    z[j] = 0.0;
                }
            }
        } else if (std::isfinite(lb[j])) {
            // lower bound only
            if (zlj >= xlj) {
                x[j] = lb[j];
                z[j] = std::max(0.0, zlj-zuj);
            } else {
                x[j] = xj;
                z[j] = 0.0;
            }
        } else if (std::isfinite(ub[j])) {
            // upper bound only
            if (zuj >= xuj) {
                x[j] = ub[j];
                z[j] = std::min(0.0, zlj-zuj);
            } else {
                x[j] = xj;
                z[j] = 0.0;
            }
        } else {
            // free variable
            x[j] = xj;
            z[j] = 0.0;
        }
    }
}

void Iterate::ResidualsFromDropping(double* pres, double* dres) const {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const SparseMatrix& AI = model_.AI();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    double presmax = 0.0;
    double dresmax = 0.0;

    for (Int j = 0; j < n+m; j++) {
        double xdrop = 0.0;     // xnew = xold - xdrop
        double zdrop = 0.0;
        switch (variable_state_[j]) {
        case StateDetail::BARRIER_LB:
            if (zl_[j] >= xl_[j])
                xdrop = x_[j]-lb[j];   // active at lower bound
            else
                zdrop = zl_[j]-zu_[j]; // inactive
            break;
        case StateDetail::BARRIER_UB:
            if (zu_[j] >= xu_[j])
                xdrop = x_[j]-ub[j];   // active at upper bound
            else
                zdrop = zl_[j]-zu_[j]; // inactive
            break;
        case StateDetail::BARRIER_BOXED:
            if (zl_[j]/xl_[j] >= zu_[j]/xu_[j]) {
                if (zl_[j] >= xl_[j])
                    xdrop = x_[j]-lb[j];   // active at lower bound
                else
                    zdrop = zl_[j]-zu_[j]; // inactive
            } else {
                if (zu_[j] >= xu_[j])
                    xdrop = x_[j]-ub[j];   // active at upper bound
                else
                    zdrop = zl_[j]-zu_[j]; // inactive
            }
            break;
        default:
            break;
        }
        double amax = 0.0;
        for (Int p = AI.begin(j); p < AI.begin(j+1); p++)
            amax = std::max(amax, std::abs(AI.value(p)));
        presmax = std::max(presmax, std::abs(xdrop) * amax);
        dresmax = std::max(dresmax, std::abs(zdrop));
    }
    if (pres)
        *pres = presmax;
    if (dres)
        *dres = dresmax;
}

void Iterate::assert_consistency() {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    assert(AllFinite(x_));
    assert(AllFinite(y_));
    assert(AllFinite(zl_));
    assert(AllFinite(zu_));

    for (Int j = 0; j < n+m; j++) {
        switch (variable_state_[j]) {
        case StateDetail::BARRIER_LB:
            assert(std::isfinite(lb[j]));
            assert(std::isinf(ub[j]));
            assert(std::isfinite(xl_[j]) && xl_[j] > 0.0);
            assert(std::isinf(xu_[j]));
            assert(zl_[j] > 0.0);
            assert(zu_[j] == 0.0);
            break;
        case StateDetail::BARRIER_UB:
            assert(std::isinf(lb[j]));
            assert(std::isfinite(ub[j]));
            assert(std::isinf(xl_[j]));
            assert(std::isfinite(xu_[j]) && xu_[j] > 0.0);
            assert(zl_[j] == 0.0);
            assert(zu_[j] > 0.0);
            break;
        case StateDetail::BARRIER_BOXED:
            assert(std::isfinite(lb[j]));
            assert(std::isfinite(ub[j]));
            assert(std::isfinite(xl_[j]) && xl_[j] > 0.0);
            assert(std::isfinite(xu_[j]) && xu_[j] > 0.0);
            assert(zl_[j] > 0.0);
            assert(zu_[j] > 0.0);
            break;
        case StateDetail::FREE:
            assert(std::isinf(lb[j]));
            assert(std::isinf(ub[j]));
            assert(std::isinf(xl_[j]));
            assert(std::isinf(xu_[j]));
            assert(zl_[j] == 0.0);
            assert(zu_[j] == 0.0);
            break;
        case StateDetail::FIXED:
            assert(xl_[j] == 0.0);
            assert(xu_[j] == 0.0);
            assert(zl_[j] == 0.0);
            assert(zu_[j] == 0.0);
            break;
        case StateDetail::IMPLIED_LB:
            assert(std::isfinite(lb[j]));
            assert(std::isinf(xl_[j]));
            assert(std::isinf(xu_[j]));
            assert(zl_[j] >= 0.0);
            assert(zu_[j] >= 0.0);
            break;
        case StateDetail::IMPLIED_UB:
            assert(std::isfinite(ub[j]));
            assert(std::isinf(xl_[j]));
            assert(std::isinf(xu_[j]));
            assert(zl_[j] >= 0.0);
            assert(zu_[j] >= 0.0);
            break;
        case StateDetail::IMPLIED_EQ:
            assert(lb[j] == ub[j]);
            assert(std::isinf(xl_[j]));
            assert(std::isinf(xu_[j]));
            assert(zl_[j] == 0.0);
            assert(zu_[j] == 0.0);
            break;
        default:
            assert(false);
        }
    }
    (void)(lb);
    (void)(ub);
}
 
void Iterate::Evaluate() const {
    if (!evaluated_) {
        ComputeResiduals();
        ComputeObjectives();
        ComputeComplementarity();
        evaluated_ = true;
    }
}

void Iterate::ComputeResiduals() const {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    const SparseMatrix& AI = model_.AI();

    // Primal residual: rb = b-AI*x.
    rb_ = model_.b();
    MultiplyAdd(AI, x_, -1.0, rb_, 'N');

    // Dual residual: rc = c-AI'y-zl+zu. If the iteate has not been
    // postprocessed, then the dual residual for fixed variables is zero
    // because these variables are treated as non-existent by the IPM.
    rc_ = model_.c() - zl_ + zu_;
    MultiplyAdd(AI, y_, -1.0, rc_, 'T');
    if (!postprocessed_) {
        for (Int j = 0; j < n+m; j++)
            if (StateOf(j) == State::fixed)
                rc_[j] = 0.0;
    }

    // Bound residuals: rl = lb-x+xl and ru = ub-x-xu. If a variable has no
    // lower or upper bound, then the corresponding bound residual is zero.
    // Also, if the variable is implied (i.e. treated as free by the IPM) or
    // fixed (i.e. treated as non-existent by the IPM), the bound residual is
    // zero.
    // Notice that postprocessing sets xl and xu for fixed/implied variables
    // such that rl and ru are zero; hence the residuals computed here are also
    // correct if the iterate was postprocessed.
    for (Int j = 0; j < n+m; j++) {
        if (has_barrier_lb(j))
            rl_[j] = lb[j] - x_[j] + xl_[j];
        else
            rl_[j] = 0.0;
    }
    for (Int j = 0; j < n+m; j++) {
        if (has_barrier_ub(j))
            ru_[j] = ub[j] - x_[j] - xu_[j];
        else
            ru_[j] = 0.0;
    }

    assert(AllFinite(rb_));
    assert(AllFinite(rc_));
    assert(AllFinite(rl_));
    assert(AllFinite(ru_));

    presidual_ = Infnorm(rb_);
    dresidual_ = Infnorm(rc_);
    presidual_ = std::max(presidual_, Infnorm(rl_));
    presidual_ = std::max(presidual_, Infnorm(ru_));
}

void Iterate::ComputeObjectives() const {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& b = model_.b();
    const Vector& c = model_.c();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    const SparseMatrix& AI = model_.AI();

    if (postprocessed_) {
        // Compute objective values as defined for the LP model.
        offset_ = 0.0;
        pobjective_ = Dot(c, x_);
        dobjective_ = Dot(b, y_);
        for (Int j = 0; j < n+m; j++) {
            if (std::isfinite(lb[j]))
                dobjective_ += lb[j] * zl_[j];
            if (std::isfinite(ub[j]))
                dobjective_ -= ub[j] * zu_[j];
        }
    } else {
        // Compute objective values for the LP that is solved at the very moment
        // (after fixing and implying variables). The offset is such that
        // pobjective_ + offset_ is the primal objective after postprocessing.
        offset_ = 0.0;
        pobjective_ = 0.0;
        for (Int j = 0; j < n+m; j++) {
            if (StateOf(j) != State::fixed)
                pobjective_ += c[j] * x_[j];
            else
                offset_ += c[j] * x_[j];
            if (is_implied(j)) {
                // At the moment, we are solving an LP with the cost coefficient
                // for variable j decreased by zl[j]-zu[j].
                pobjective_ -= (zl_[j]-zu_[j]) * x_[j];
                offset_ += (zl_[j]-zu_[j]) * x_[j];
            }
        }
        dobjective_ = Dot(b, y_);
        for (Int j = 0; j < n+m; j++) {
            if (has_barrier_lb(j))
                dobjective_ += lb[j] * zl_[j];
            if (has_barrier_ub(j))
                dobjective_ -= ub[j] * zu_[j];
            if (StateOf(j) == State::fixed)
                // At the moment, we are solving the LP without variable j,
                // but with the RHS decreased by AI[:,j]*x[j].
                dobjective_ -= x_[j] * DotColumn(AI, j, y_);
        }
    }
}

void Iterate::ComputeComplementarity() const {
    const Int m = model_.rows();
    const Int n = model_.cols();

    complementarity_ = 0.0;
    mu_min_ = INFINITY;
    mu_max_ = 0.0;
    Int num_finite = 0;
    for (Int j = 0; j < n+m; j++) {
        if (has_barrier_lb(j)) {
            complementarity_ += xl_[j] * zl_[j];
            mu_min_ = std::min(mu_min_, xl_[j]*zl_[j]);
            mu_max_ = std::max(mu_max_, xl_[j]*zl_[j]);
            num_finite++;
        }
    }
    for (Int j = 0; j < n+m; j++) {
        if (has_barrier_ub(j)) {
            complementarity_ += xu_[j] * zu_[j];
            mu_min_ = std::min(mu_min_, xu_[j]*zu_[j]);
            mu_max_ = std::max(mu_max_, xu_[j]*zu_[j]);
            num_finite++;
        }
    }
    if (num_finite > 0)
        mu_ = complementarity_ / num_finite;
    else
        mu_ = mu_min_ = 0.0;
}

}  // namespace ipx

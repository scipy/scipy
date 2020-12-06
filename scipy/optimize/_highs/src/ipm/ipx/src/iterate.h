// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_ITERATE_H_
#define IPX_ITERATE_H_

#include <vector>
#include "model.h"

namespace ipx {

// Iterate manages the IPM iterate consisting of primal variables x[n+m],
// xl[n+m], xu[n+m] and dual variables y[m], zl[n+m], zu[n+m], where m, n are
// the # rows and structural columns of the model.
//
// Each x[j] is in one of three "states":
//
// "fixed":     The variable is not moved by the IPM. The IPM works as if b was
//              reduced by x[j]*AI[:,j] and column j of AI would not be present.
// "free":      The variable is treated by the IPM as if lb[j] and ub[j] would
//              be infinite, regardless of their actual values.
// "barrier":   The variable has at least one finite bound and a barrier term is
//              added by the IPM.
//
// Notice that the state refers to the way the variable is treated by the IPM,
// not to its bounds in the model. For example, a variable of fixed state might
// have zero, one or both finite bounds in the model.

class Iterate {
public:
    // Constructs an Iterate object associated with a model. @model must be
    // valid as long as the object is used; no data is copied. The constructor
    // initializes variables and states to be consistent with finite/infinite
    // bounds in the model.
    explicit Iterate(const Model& model);

    // Copies the given values into the Iterate object and sets the state of
    // each variable to
    // - free     if the variable has no finite bound in the model,
    // - barrier  otherwise.
    // The variables must have values that are consistent with the bounds of the
    // model; that is:
    // - If lb[j]==-inf, then   xl[j]==inf and zl[j]==0.0.
    // - If lb[j]> -inf, then 0<xl[j]< inf and zl[j]> 0.0.
    // - If ub[j]==+inf, then   xu[j]==inf and zu[j]==0.0.
    // - If ub[j]< +inf, then 0<xu[j]< inf and zu[j]> 0.0.
    // Otherwise an assertion will fail.
    void Initialize(const Vector& x, const Vector& xl, const Vector& xu,
                    const Vector& y, const Vector& zl, const Vector& zu);

    // Adds sp*[dx;dxl;dxu] to the primal iterate and sd*[dy;dzl;dzu] to the
    // dual iterate (restrictions see below). Each of the pointer arguments can
    // be NULL, in which case its entries are assumed to be 0.0.
    // - x[j] is updated only if StateOf(j) != State::fixed.
    // - xl[j] and zl[j] are updated only if has_barrier_lb(j). The update to
    //                   each of xl[j] and zl[j] is truncated such that
    //                   xl[j] >= kBarrierMin and zl[j] >= kBarrierMin.
    // - xu[j] and zu[j] are updated only if has_barrier_ub(j). The update to
    //                   each of xu[j] and zu[j] is truncated such that
    //                   xu[j] >= kBarrierMin and zu[j] >= kBarrierMin.
    // If a variable is not updated, then the entry in the step direction is not
    // accessed.
    void Update(double sp, const double* dx, const double* dxl,
                const double* dxu, double sd, const double* dy,
                const double* dzl, const double* dzu);

    const Model& model() const { return model_; }

    const Vector& x() const  { return x_; }
    const Vector& xl() const { return xl_; }
    const Vector& xu() const { return xu_; }
    const Vector& y() const  { return y_; }
    const Vector& zl() const { return zl_; }
    const Vector& zu() const { return zu_; }

    double x(Int j) const  { return x_[j]; }
    double xl(Int j) const { return xl_[j]; }
    double xu(Int j) const { return xu_[j]; }
    double y(Int j) const  { return y_[j]; }
    double zl(Int j) const { return zl_[j]; }
    double zu(Int j) const { return zu_[j]; }

    // Returns const rerefences to the residual vectors
    // rb = b-AI*x,
    // rl = lb-x+xl,
    // ru = ub-x-xu,
    // rc = c-AI'y-zl+zu.
    // Warning: The residuals are not evaluated immediately after a change to
    // an Iterate object (such as by Update() or make_fixed()). Instead, the
    // four residual vectors are evaluated when any of them is required, e.g.
    // in calls to rb(), presidual(), etc. Hence if a reference to a residual
    // vector is stored and the Iterate object is changed, you must call e.g.
    // presidual() to ensure that the residuals are updated.
    const Vector& rb() const;
    const Vector& rl() const;
    const Vector& ru() const;
    const Vector& rc() const;

    enum class State { fixed, free, barrier };

    // Returns the state of variable j.
    State StateOf(Int j) const;

    // Returns true if a barrier term exists for the lower or upper bound.
    // The following expression is true for each variable:
    // (StateOf(j)==State::barrier) == (has_barrier_lb(j)||has_barrier_ub(j))
    bool has_barrier_lb(Int j) const;
    bool has_barrier_ub(Int j) const;

    // Returns true if variable j is "implied", i.e. StateOf(j)==State::free
    // but the variable has at least one finite bound in the model.
    bool is_implied(Int j) const;

    // Changes state of variable j to fixed and sets xl[j]=xu[j]=zl[j]=zu[j]=0.
    // If @value is given, sets x[j]=@value.
    void make_fixed(Int j);
    void make_fixed(Int j, double value);

    // Changes the state of variable j to free.
    // make_implied_lb(j) lb[j] must be finite.
    //                    The method leaves zl[j] and zu[j] unchanged. The cost
    //                    coefficient of variable j in the remaining LP becomes
    //                    c[j]-zl[j]+zu[j].
    //                    Postprocess() sets x[j]=lb[j], zu[j]=0 and computes
    //                    zl[j] to satisfy dual feasibility.
    // make_implied_ub(j) ub[j] must be finite.
    //                    The method leaves zl[j] and zu[j] unchanged. The cost
    //                    coefficient of variable j in the remaining LP becomes
    //                    c[j]-zl[j]+zu[j].
    //                    Postprocess() sets x[j]=ub[j], zl[j]=0 and computes
    //                    zu[j] to satisfy dual feasibility.
    // make_implied_eq(j) must have lb[j]==ub[j].
    //                    The method sets zl[j]=0 and zu[j]=0.
    //                    Postprocess sets x[j]=lb[j] and chooses either zl[j]
    //                    or zu[j] to be nonzero depending on the sign of the
    //                    reduced cost.
    void make_implied_lb(Int j);
    void make_implied_ub(Int j);
    void make_implied_eq(Int j);

    // Returns the IPM scaling factor for variable j, which is
    // 0                                 if StateOf(j)==State::fixed,
    // inf                               if StateOf(j)==State::free,
    // 1.0/sqrt(zl[j]/xl[j]+zu[j]/xu[j]) if StateOf(j)==State::barrier.
    // In the latter case the scaling factor is finite and positive.
    double ScalingFactor(Int j) const;

    // Returns the primal/dual objective value in the LP that is currently
    // being solved (after fixing and implying variables).
    double pobjective() const;
    double dobjective() const;

    // Returns the primal/dual objective value obtained after postprocessing.
    double pobjective_after_postproc() const;
    double dobjective_after_postproc() const;

    // Returns the infinity norm of [rb;rl;ru] and rc.
    double presidual() const;
    double dresidual() const;

    // copmlementarity() returns the sum of the pairwise complementarity
    // products xl[j]*zl[j] and xu[j]*zu[j] from all barrier terms. mu()
    // returns the average, mu_min() the minimum and mu_max() the maximum.
    double complementarity() const;
    double mu() const;
    double mu_min() const;
    double mu_max() const;

    // Returns true if the relative primal and dual residuals are <=
    // feasibility_tol().
    bool feasible() const;

    // Returns true if the relative objective gap is <= optimality_tol().
    bool optimal() const;

    // Returns true if the iterate satisfies the IPM termination criterion.
    // If crossover_start() <= 0, the termination criterion is feasible() &&
    // optimal(). If crossover_start() > 0, it is additionally required that
    // the relative residuals which result from dropping the iterate to
    // complementarity are <= crossover_start().
    bool term_crit_reached() const;

    double feasibility_tol() const { return feasibility_tol_; }
    double optimality_tol() const { return optimality_tol_; }
    double crossover_start() const { return crossover_start_; }

    void feasibility_tol(double new_tol) { feasibility_tol_ = new_tol; }
    void optimality_tol(double new_tol) { optimality_tol_ = new_tol; }
    void crossover_start(double new_tol) { crossover_start_ = new_tol; }

    // Substitutes x[j], xl[j], xu[j], zl[j] and zu[j] for fixed and implied
    // variables as defined by their state. After Postprocess() the object is
    // invalid for use as an IPM iterate.
    void Postprocess();

    // Calls Model::EvaluateInteriorSolution().
    void EvaluatePostsolved(Info* info) const;

    // Copies y and constructs x and z from the iterate such that
    // x[j]==lb[j] || x[j]==ub[j] || z[j]==0.0 is true for each j.
    // The method can only be called after Postprocess().
    void DropToComplementarity(Vector& x, Vector& y, Vector& z) const;

private:
    // A (primal or dual) variable that is required to be positive in the IPM is
    // not moved closer to zero than kBarrierMin.
    static constexpr double kBarrierMin = 1e-30;

    void assert_consistency();
    void Evaluate() const;
    void ComputeResiduals() const;
    void ComputeObjectives() const;
    void ComputeComplementarity() const;

    // Computes the maximum primal and dual residual that results from dropping
    // a barrier variable. Dropping means to set either x[j] to a bound or zl[j]
    // and zu[j] to zero, where the decision is made by the ratios or primal to
    // dual variables.
    void ResidualsFromDropping(double* pres, double* dres) const;

    // Internally, the state of each variable is subdivided as follows:
    //
    // State::barrier is represented by StateDetail::
    //
    //     BARRIER_LB   0 < xl < inf, xu = inf, zl > 0, zu = 0.
    //                  The variable is moved by the barrier solver. It has a
    //                  finite lower bound and an infinite upper bound.
    //    
    //     BARRIER_UB   xl = inf, 0 < xu < inf, zl = 0, zu > 0.
    //                  The variable is moved by the barrier solver. It has an
    //                  infinite lower bound and a finite upper bound.
    //    
    //     BARRIER_BOXED 0 < xl < inf, 0 < xu < inf, zl > 0, zu > 0.
    //                   The variable is moved by the barrier solver. It has a
    //                   finite lower bound and a finite upper bound.
    //
    // State::fixed is represented by StateDetail::
    //
    //     FIXED        xl = 0, xu = 0, zl = 0, zu = 0.
    //                  The variable is fixed at its current x value.
    //
    // State::free is represented by StateDetail::
    //
    //     FREE         xl = inf, xu = inf, zl = 0, zu = 0.
    //                  The variable is moved by the barrier solver. It has no
    //                  finite bound.
    //    
    //     IMPLIED_LB   xl = inf, xu = inf, zl >= 0, zu >= 0
    //                  The variable has a finite lb but is treated as free by
    //                  the barrier solver. (The ub may or may not be finite.)
    //                  After termination we set x = lb, zu = 0 and compute zl
    //                  from dual feasibility.
    //                  Evaluate() subtracts (zl-zu) from the cost coefficient
    //                  when computing the primal objective value.
    //    
    //     IMPLIED_UB   xl = inf, xu = inf, zl >= 0, zu >= 0
    //                  The variable has a finite ub but is treated as free by
    //                  the barrier solver. (The lb may or may not be finite.)
    //                  After termination we set x = ub, zl = 0 and compute zu
    //                  from dual feasibility.
    //                  Evaluate() subtracts (zl-zu) from the cost coefficient
    //                  when computing the primal objective value.
    //    
    //     IMPLIED_EQ   xl = inf, xu = inf, zl = 0, zu = 0
    //                  The variable has lb == ub but is treated as free by
    //                  the barrier solver. After termination we set x = lb
    //                  and compute either zl or zu from dual feasibility and
    //                  set the other one to zero.
    //
    enum class StateDetail { BARRIER_LB, BARRIER_UB, BARRIER_BOXED, FREE, FIXED,
                             IMPLIED_LB, IMPLIED_UB, IMPLIED_EQ };

    const Model& model_;
    Vector x_, xl_, xu_, y_, zl_, zu_;
    std::vector<StateDetail> variable_state_;

    // The mutable quantities have a defined value if evaluated_ is true.
    // They are computed by Evaluate(), which is called by any method that
    // accesses any of these quantities.
    mutable Vector rb_, rl_, ru_, rc_;
    mutable double pobjective_{0.0};
    mutable double dobjective_{0.0};
    mutable double presidual_{0.0};
    mutable double dresidual_{0.0};
    mutable double offset_{0.0};
    mutable double complementarity_{0.0};
    mutable double mu_{0.0};
    mutable double mu_min_{0.0};
    mutable double mu_max_{0.0};
    mutable bool evaluated_{false};
    bool postprocessed_{false};

    double feasibility_tol_{1e-6};
    double optimality_tol_{1e-8};
    double crossover_start_{-1.0};
};

inline Iterate::State Iterate::StateOf(Int j) const {
    switch (variable_state_[j]) {
    case StateDetail::FIXED:
        return State::fixed;
    case StateDetail::FREE:
    case StateDetail::IMPLIED_LB:
    case StateDetail::IMPLIED_UB:
    case StateDetail::IMPLIED_EQ:
        return State::free;
    default:
        return State::barrier;
    }
}

inline bool Iterate::has_barrier_lb(Int j) const {
    return variable_state_[j] == StateDetail::BARRIER_LB
        || variable_state_[j] == StateDetail::BARRIER_BOXED;
}

inline bool Iterate::has_barrier_ub(Int j) const{
    return variable_state_[j] == StateDetail::BARRIER_UB
        || variable_state_[j] == StateDetail::BARRIER_BOXED;
}

inline bool Iterate::is_implied(Int j) const {
    return variable_state_[j] == StateDetail::IMPLIED_LB
        || variable_state_[j] == StateDetail::IMPLIED_UB
        || variable_state_[j] == StateDetail::IMPLIED_EQ;
}

}  // namespace ipx

#endif  // IPX_ITERATE_H_

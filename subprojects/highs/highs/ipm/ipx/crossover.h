#ifndef IPX_CROSSOVER_H_
#define IPX_CROSSOVER_H_

// Crossover method
//
// Given a basis and (x,y,z) that satisfies
//
//   lb[j] <= x[j] <= ub[j],                    (1)
//    z[j] <= 0    if x[j] > lb[j],             (2a)
//    z[j] >= 0    if x[j] < ub[j]              (2b)
//
// for all j, the crossover method updates the basis and (x,y,z) such that
// z[j]=0 for all basic variables and x[j]=lb[j] or x[j]=ub[j] for all nonbasic
// variables (or x[j]=0 if the variable is free). Hereby (1) and (2) are
// maintained. The textbook method also keeps Ax and A'y+z unchanged.
//
// Nonbasic variables for which lb[j]<x[j]<ub[j] are called "primal superbasic"
// and basic variables for which z[j]!=0 are called "dual superbasic". The
// crossover algorithm removes superbasic variables in two push phases.
//
// Each iteration of the primal push phase chooses a primal superbasic variable
// jn, moves x[jn] toward a bound (we choose the nearer one if jn has two finite
// bounds and zero if it has none) and updates x[basic] to keep Ax unchanged. If
// jn reaches its bound, then the push is complete. Otherwise a basic variable
// jb reached its bound and blocked the step. In this case a basis update
// exchanges jb by jn.
//
// Each iteration of the dual push phase chooses a dual superbasic variable jb,
// moves z[jb] toward zero and updates z[nonbasic] to keep A'y+z unchanged. If
// jb reaches zero, then the push is complete. Otherwise a nonbasic variable jn
// became zero and blocked the step. In this case a basis update exchanges jb by
// jn.

#include <vector>
#include "ipm/ipx/basis.h"
#include "ipm/ipx/control.h"
#include "ipm/ipx/indexed_vector.h"

namespace ipx {

class Crossover {
public:
    // Constructor stores a reference to the Control object, which must be valid
    // as long as the Crossover object is used.
    Crossover(const Control& control);

    // First runs the dual push phase; if this was successful, then runs the
    // primal push phase.
    //
    // weights: Must either be NULL or an array of size n+m.
    //          If given, then primal pushes are executed in decreasing order of
    //          weights and dual pushes in increasing order.
    //          If NULL, then both phases push variables in increasing order of
    //          index.
    //
    // On return info->status_crossover and info->errflag have been set by
    // PushPrimal() or PushDual().
    //
    void PushAll(Basis* basis, Vector& x, Vector& y, Vector& z,
                 const double* weights, Info* info);

    // Pushes a set of primal variables to one of their bounds.
    //
    // basis:           the basis to be updated; also provides the LP model.
    // x:               size n+m vector, satisfying (1) on entry and return.
    //                  x is updated while A*x (ideally) remains unchanged.
    // variables:       indices of nonbasic variables (without duplicates) to be
    //                  pushed to a bound.
    // fixed_at_bound:  NULL or size n+m array. If fixed_at_bound[j] = true,
    //                  then x[j] must be at a bound on entry and remains at
    //                  this bound. If NULL, then no variable is fixed.
    // info:            on return info->status_crossover is one of
    //                  * IPX_STATUS_optimal      if terminated successfully,
    //                  * IPX_STATUS_time_limit   if interrupted,
    //                  * IPX_STATUS_failed       if failed.
    //                  In the latter case info->errflag is set.
    //
    // If a variable to be pushed has two finite bounds, then the nearer one is
    // chosen. If it has no finite bound, then it is pushed to zero.
    //
    void PushPrimal(Basis* basis, Vector& x, const std::vector<Int>& variables,
                    const bool* fixed_at_bound, Info* info);

    // As above, but with fixed_at_bound defined through z.
    // Variable j is fixed at its bound if z[j]!=0.
    void PushPrimal(Basis* basis, Vector& x, const std::vector<Int>& variables,
                    const Vector& z, Info* info);

    // Pushes a set of dual variables to zero.
    //
    // basis:           the basis to be updated; also provides the LP model.
    // y, z:            size m and n+m vectors that are updated while A'y+z
    //                  (ideally) remains unchanged.
    // variables:       indices of basic variables (without duplicates) to be
    //                  pushed to zero.
    // sign_restrict:   size n+m array. z[j] is restricted to be
    //                  * non-negative iff (sign_restrict[j] & 1) != 0,
    //                  * non-positive iff (sign_restrict[j] & 2) != 0.
    //                  The sign restriction must be satisfied on entry and is
    //                  satisfied on return.
    // info:            on return info->status_crossover is one of
    //                  * IPX_STATUS_optimal      if terminated successfully,
    //                  * IPX_STATUS_time_limit   if interrupted,
    //                  * IPX_STATUS_failed       if failed.
    //                  In the latter case info->errflag is set.
    //
    void PushDual(Basis* basis, Vector& y, Vector& z,
                  const std::vector<Int>& variables,
                  const int sign_restrict[], Info* info);

    // As above, but with the sign restriction on z defined through x.
    // z[j] is restricted to be
    // * non-negative iff x[j] < ub[j],
    // * non-positive iff x[j] > lb[j].
    void PushDual(Basis* basis, Vector& y, Vector& z,
                  const std::vector<Int>& variables,
                  const Vector& x, Info* info);

    // Number of pushes in last call to PushPrimal() and PushDual().
    Int primal_pushes() const { return primal_pushes_; }
    Int dual_pushes() const { return dual_pushes_; }

    // Number of basis updates in last call to PushPrimal() and PushDual().
    Int primal_pivots() const { return primal_pivots_; }
    Int dual_pivots() const { return dual_pivots_; }

    // Runtime of last call to PushPrimal() and PushDual().
    double time_primal() const { return time_primal_; }
    double time_dual() const { return time_dual_; }

private:
    // An entry of the tableau row/column is eligible as pivot only if it is
    // larger than kPivotZeroTol in absolute value.
    static constexpr double kPivotZeroTol = 1e-5;

    // Two-pass ratio tests that allow infeasibilities up to feastol in order
    // to choose a larger pivot.
    Int PrimalRatioTest(const Vector& xbasic, const IndexedVector& ftran,
                        const Vector& lbbasic, const Vector& ubbasic,
                        double step, double feastol, bool* block_at_lb);
    Int DualRatioTest(const Vector& z, const IndexedVector& row,
                      const int sign_restrict[], double step,
                      double feastol);

    const Control& control_;

    Int primal_pushes_{0};
    Int dual_pushes_{0};
    Int primal_pivots_{0};
    Int dual_pivots_{0};
    double time_primal_{0.0};
    double time_dual_{0.0};
};

}  // namespace ipx

#endif  // IPX_CROSSOVER_H_

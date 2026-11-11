#ifndef IPX_CONJUGATE_RESIDUALS_H_
#define IPX_CONJUGATE_RESIDUALS_H_

// Implementation of the (preconditioned) Conjugate Residuals (CR) method for
// symmetric positive definite linear systems. Without preconditioning, the
// method is implemented as in [1, Algorithm 6.20]. The implementation with
// preconditioning is described in [2, Section 6.3].
//
// [1] Y. Saad, "Iterative Methods for Sparse Linear Systems", 2nd (2003)
// [2] L. Schork, "Basis Preconditioning in Interior Point Methods", PhD thesis
//     (2018)

#include "ipm/ipx/control.h"
#include "ipm/ipx/linear_operator.h"

namespace ipx {

class ConjugateResiduals {
public:
    // Constructs a ConjugateReisdual object. The object does not have any
    // memory allocated between calls to Solve().
    // @control for InterruptCheck() and Debug(). No parameters are accessed.
    ConjugateResiduals(const Control& control);

    // Solves C*lhs = rhs. @lhs has initial iterate on entry, solution on
    // return. The method terminates when reaching the accuracy criterion
    //
    //   Infnorm(residual) <= tol               (if resscale == NULL), or
    //   Infnorm(resscale.*residual) <= tol     (if resscale != NULL).
    //
    // In the latter case, @resscale must be an array of dimension @rhs.size().
    // The method also stops after @maxiter iterations. If @maxiter < 0, a
    // maximum of @rhs.size()+100 iterations is performed. (In exact arithmetic
    // the solution would be found after @rhs.size() iterations. It happened on
    // some LP models with m << n, e.g. "rvb-sub" from MIPLIB2010, that the CR
    // method did not reach the termination criterion within m iterations,
    // causing the IPM to fail. Giving the CR method 100 extra iterations
    // resolved the issue on all LP models from our test set where it occurred.)
    //
    // If the @P argument is given, it is used as preconditioner (which
    // approximates inverse(C)) and must be symmetric positive definite.
    //
    void Solve(LinearOperator& C, const Vector& rhs,
               double tol, const double* resscale, Int maxiter, Vector& lhs);
    void Solve(LinearOperator& C, LinearOperator& P, const Vector& rhs,
               double tol, const double* resscale, Int maxiter, Vector& lhs);

    // Returns 0 if the last call to Solve() terminated successfully (i.e.
    // the system was solved to the required accuracy). Otherwise returns
    // IPX_ERROR_cr_iter_limit          if iteration limit was reached
    // IPX_ERROR_cr_matrix_not_posdef   if v'*C*v <= 0 for some vector v
    // IPX_ERROR_cr_precond_not_posdef  if v'*P*v <= 0 for some vector v
    // IPX_ERROR_cr_inf_or_nan          if overflow occurred
    // IPX_ERROR_cr_no_progress         if no progress due to round-off errors
    // IPX_ERROR_user_interrupt         if interrupted by user in control
    // IPX_ERROR_time_interrupt         if interrupted by time limit in control
    Int errflag() const;

    // Returns the # iterations in the last call to Solve().
    Int iter() const;

    // Returns the runtime of the last call to Solve().
    double time() const;

private:
    const Control& control_;
    Int errflag_{0};
    Int iter_{0};
    double time_{0.0};
};

}  // namespace ipx

#endif  // IPX_CONJUGATE_RESIDUALS_H_

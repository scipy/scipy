#ifndef IPX_KKT_SOLVER_BASIS_H_
#define IPX_KKT_SOLVER_BASIS_H_

#include "ipm/ipx/basis.h"
#include "ipm/ipx/control.h"
#include "ipm/ipx/kkt_solver.h"
#include "ipm/ipx/model.h"
#include "ipm/ipx/splitted_normal_matrix.h"

namespace ipx {

// KKTSolverBasis implements a KKT solver using the CR method with basis
// preconditioning. Maintaining a basis matrix allows to remove degenerate
// variables from the optimization process. This technique is applied in
// Factorize() if parameter ipm_drop_primal or ipm_drop_dual is positive.
//
// In the call to Factorize() @iterate must not be NULL.

class KKTSolverBasis : public KKTSolver {
public:
    // Constructor wraps a Basis object into a KKTSolverBasis object. The basis
    // must have been initialized to a starting basis and must be valid as long
    // as the KKT solver is used. Its associated model defines the LP data.
    KKTSolverBasis(const Control& control, Basis& basis);

    Int maxiter() const { return maxiter_; }
    void maxiter(Int new_maxiter) { maxiter_ = new_maxiter; }

private:
    // DropPrimal() and DropDual() do not choose entries as pivots that are <=
    static constexpr double kPivotZeroTol = 1e-7;

    void _Factorize(Iterate* iterate, Info* info) override;
    void _Solve(const Vector& a, const Vector& b, double tol,
                Vector& x, Vector& y, Info* info) override;
    Int _iter() const override { return iter_; }
    Int _basis_changes() const override { return basis_changes_; }
    const Basis* _basis() const override { return &basis_; }

    // Processes basic variables that are close to a bound by either pivoting
    // them out of the basis or making them "implied"; i.e. the dual is fixed at
    // its current value and the variable is treated as free by the interior
    // point solver. After the LP solve finished, the primal is set to its
    // bound.
    void DropPrimal(Iterate* iterate, Info* info);

    // Processes nonbasic variables whose dual is close to zero by either
    // pivoting them into the basis or "fixing" the primal; i.e. the variable is
    // removed from optimization and the primal remains at its current value.
    // The dual is set to zero.
    void DropDual(Iterate* iterate, Info* info);

    const Control& control_;
    const Model& model_;
    Basis& basis_;
    SplittedNormalMatrix splitted_normal_matrix_;
    Vector colscale_;           // interior point column scaling factors
    bool factorized_{false};    // preconditioner prepared?
    Int maxiter_{-1};
    Int iter_{0};
    Int basis_changes_{0};
};

}  // namespace ipx

#endif  // IPX_KKT_SOLVER_BASIS_H_

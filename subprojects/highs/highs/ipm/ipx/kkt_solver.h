#ifndef IPX_KKT_SOLVER_H_
#define IPX_KKT_SOLVER_H_

#include "ipm/ipx/basis.h"
#include "ipm/ipx/iterate.h"

namespace ipx {

// Interface to KKT solver implementations. A KKT solver implements a direct or
// iterative method for solving
//
//   [ G   AI' ] (x) = (a) ,                        (1)
//   [ AI   0  ] (y)   (b)
//
// where AI is the m-by-(n+m) matrix defined by a model and G is a positive
// semidefinite diagonal matrix. The solver may add regularization to G and/or
// the zero block.
//
// An iterative solver must compute an approximate solution of the form
//
//   [ G   AI' ] (x) = (a) + (res)                  (2)
//   [ AI   0  ] (y)   (b)   ( 0 )
//
// that satisfies Infnorm(D*res) <= tol, where D is the diagonal matrix with
// entries D[i,i] = sqrt(1/G[i,i]) if G[i,i] != 0 and D[i,i] = 1 otherwise.

class KKTSolver {
public:
    KKTSolver& operator=(const KKTSolver&) = delete;
    KKTSolver& operator=(KKTSolver&&) = delete;
    virtual ~KKTSolver() {}

    // Factorizes the KKT matrix (direct solver) or prepares preconditioner
    // (iterative solver). The diagonal matrix G is built from @iterate.
    // The implementation is allowed to change variable statuses to eliminate
    // close-to-converged variables from the IPM solve. Some implementations
    // allow @iterate to be NULL, in which case G is assumed to be the identity
    // matrix.
    void Factorize(Iterate* iterate, Info* info);

    // Solves KKT system. If an iterative method is used, @tol is the required
    // tolerance for the residual in (2) as specified above.
    void Solve(const Vector& a, const Vector& b, double tol,
               Vector& x, Vector& y, Info* info);

    // If an iterative method is used, returns the # iterations in all Solve()
    // calls since the last call to Factorize(). A direct solver returns the #
    // iterative refinement steps.
    Int iter() const;

    // If a basis matrix is maintained, returns the # basis changes in the last
    // call to Factorize(). Otherwise returns 0.
    Int basis_changes() const;

    // If a basis matrix is maintained, returns a pointer to it.
    // Otherwise returns NULL.
    const Basis* basis() const;

private:
    virtual void _Factorize(Iterate* iterate, Info* info) = 0;
    virtual void _Solve(const Vector& a, const Vector& b, double tol,
                         Vector& x, Vector& y, Info* info) = 0;
    virtual Int _iter() const = 0;
    virtual Int _basis_changes() const { return 0; }
    virtual const Basis* _basis() const { return nullptr; }
};

}  // namespace ipx

#endif  // IPX_KKT_SOLVER_H_

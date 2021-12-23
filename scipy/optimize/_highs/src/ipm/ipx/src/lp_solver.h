// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_LP_SOLVER_H_
#define IPX_LP_SOLVER_H_

#include <memory>
#include "basis.h"
#include "control.h"
#include "ipm.h"
#include "iterate.h"
#include "model.h"

namespace ipx {

class LpSolver {
public:
    // Solves an LP problem in the form given in the reference documentation.
    // @num_var: number of variables, must be > 0.
    // @obj: size num_var array of objective coefficients.
    // @lb: size num_var array of variable lower bounds, can have -INFINITY.
    // @lb: size num_var array of variable upper bounds, can have +INFINITY.
    // @num_constr: number of constraints, must be >= 0.
    // @Ap, @Ai, @Ax: constraint matrix in CSC format; indices can be unsorted.
    // @rhs: size num_constr array of right-hand side entries.
    // @constr_type: size num_constr array of entries '>', '<' and '='.
    // Returns GetInfo().status.
    Int Solve(Int num_var, const double* obj, const double* lb,
              const double* ub, Int num_constr, const Int* Ap, const Int* Ai,
              const double* Ax, const double* rhs, const char* constr_type);

    // Returns the solver info from the last call to Solve(). See the reference
    // documentation for the meaning of Info values.
    Info GetInfo() const;

    // Returns the final IPM iterate from the last call to Solve() into user
    // arrays. An iterate is available if GetInfo().status_ipm !=
    // IPX_STATUS_not_run. If no iterate is available, the method does nothing.
    // Each of the pointer arguments must either be NULL or an array of
    // appropriate dimension. If NULL, the quantity is not returned.
    // Returns -1 if no IPM iterate was available and 0 otherwise.
    Int GetInteriorSolution(double* x, double* xl, double* xu, double* slack,
                            double* y, double* zl, double* zu) const;

    // Returns the basic solution and basis from the last call to Solve() into
    // user arrays. A basic solution and basis are available if
    // GetInfo().status_crossover == IPX_STATUS_optimal ||
    // GetInfo().status_crossover == IPX_STATUS_imprecise. Otherwise the method
    // does nothing. Each of the pointer arguments must either be NULL or an
    // array of appropriate dimension. If NULL, the quantity is not returned.
    // Returns -1 if no basic solution was available and 0 otherwise.
    Int GetBasicSolution(double* x, double* slack, double* y, double* z,
                         Int* cbasis, Int* vbasis) const;

    // Returns/sets all paramters. Without calling SetParameters(), the solver
    // uses the default values of a Parameters object.
    Parameters GetParameters() const;
    void SetParameters(Parameters new_parameters);

    // Discards the model and solution (if any) but keeps the parameters.
    void ClearModel();

    // -------------------------------------------------------------------------
    // The remaining methods are for debugging.
    // -------------------------------------------------------------------------

    // Returns the current IPM iterate without postsolve. The method does
    // nothing when no iterate is available (i.e. when IPM was not started).
    // @x, @xl, @xu, @zl, @zu: either NULL or size num_cols_solver arrays.
    // @y: either NULL or size num_rows_solver array.
    // Returns -1 if no IPM iterate was available and 0 otherwise.
    Int GetIterate(double* x, double* y, double* zl, double* zu, double* xl,
                   double* xu);

    // Returns the current basis postsolved.
    // - If crossover terminated successfully, this is the basis returned by
    //   GetBasicSolution().
    // - If crossover failed, this is the basis at which failure occured.
    // - If crossover was not called, this is the basis from the IPM
    //   preconditioner.
    // - If no basis is available, the method does nothing.
    // @cbasis: either NULL or size num_constr array.
    // @vbasis: either NULL or size num_var array.
    // Returns -1 if no basis was available and 0 otherwise.
    Int GetBasis(Int* cbasis, Int* vbasis);

    // Returns the constraint matrix from the solver (including slack columns)
    // and the diagonal from the (1,1) block of the KKT matrix corresponding to
    // the current IPM iterate. The method does nothing when no IPM iterate is
    // available (i.e. when IPM was not started).
    // @AIp: either NULL or size num_cols_solver + 1 array.
    // @AIi, @AIx: either NULL or size num_entries_solver arrays.
    // (If any of the three arguments is NULL, the matrix is not returned.)
    // @g: either NULL or size num_cols_solver array.
    // Returns -1 if no IPM iterate was available and 0 otherwise.
    Int GetKKTMatrix(Int* AIp, Int* AIi, double* AIx, double* g);

    // (Efficiently) computes the number of nonzeros per row and column of the
    // symbolic inverse of the basis matrix.
    // @rowcounts, @colcounts: either NULL or size num_rows_solver arrays.
    // Returns -1 if no basis was available and 0 otherwise.
    Int SymbolicInvert(Int* rowcounts, Int* colcounts);

private:
    void InteriorPointSolve();
    void RunIPM();
    void ComputeStartingPoint(IPM& ipm);
    void RunInitialIPM(IPM& ipm);
    void BuildStartingBasis();
    void RunMainIPM(IPM& ipm);
    void RunCrossover();
    void PrintSummary();

    Control control_;
    Info info_;
    Model model_;
    std::unique_ptr<Iterate> iterate_;
    std::unique_ptr<Basis> basis_;

    // Basic solution computed by crossover and basic status of each variable
    // (one of IPX_nonbasic_lb, IPX_nonbasic_ub, IPX_basic, IPX_superbasic).
    // If crossover was not run or failed, then basic_statuses_ is empty.
    Vector x_crossover_, y_crossover_, z_crossover_;
    std::vector<Int> basic_statuses_;
};

}  // namespace ipx

#endif  // IPX_LP_SOLVER_H_

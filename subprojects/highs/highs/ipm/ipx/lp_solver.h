#ifndef IPX_LP_SOLVER_H_
#define IPX_LP_SOLVER_H_

#include <memory>
#include "ipm/ipx/basis.h"
#include "ipm/ipx/control.h"
#include "ipm/ipx/ipm.h"
#include "ipm/ipx/iterate.h"
#include "ipm/ipx/model.h"
#include "lp_data/HighsCallback.h"

namespace ipx {

class LpSolver {
public:
    // Loads an LP model in the form given in the reference documentation.
    // @num_var: number of variables, must be > 0.
    // @obj: size num_var array of objective coefficients.
    // @lb: size num_var array of variable lower bounds, can have -INFINITY.
    // @lb: size num_var array of variable upper bounds, can have +INFINITY.
    // @num_constr: number of constraints, must be >= 0.
    // @Ap, @Ai, @Ax: constraint matrix in CSC format; indices can be unsorted.
    // @rhs: size num_constr array of right-hand side entries.
    // @constr_type: size num_constr array of entries '>', '<' and '='.
    // Returns:
    //  0
    //  IPX_ERROR_argument_null
    //  IPX_ERROR_invalid_dimension
    //  IPX_ERROR_invalid_matrix
    //  IPX_ERROR_invalid_vector
    Int LoadModel(Int num_var, const double offset,
		  const double* obj, const double* lb,
                  const double* ub, Int num_constr, const Int* Ap,
                  const Int* Ai, const double* Ax, const double* rhs,
                  const char* constr_type);

    // Loads a primal-dual point as starting point for the IPM.
    // @x: size num_var array
    // @xl: size num_var array, must satisfy xl[j] >= 0 for all j and
    //      xl[j] == INFINITY iff lb[j] == -INFINITY.
    // @xu: size num_var array, must satisfy xu[j] >= 0 for all j and
    //      xu[j] == INFINITY iff ub[j] == INFINITY.
    // @slack: size num_constr array, must satisfy
    //         slack[i] == 0 if constr_type[i] == '='
    //         slack[i] >= 0 if constr_type[i] == '<'
    //         slack[i] <= 0 if constr_type[i] == '>'
    // @y: size num_constr array, must satisfy
    //     y[i] >= 0 if constr_type[i] == '>'
    //     y[i] <= 0 if constr_type[i] == '<'
    // @zl: size num_var array, must satisfy zl[j] >= 0 for all j and
    //      zl[j] == 0 if lb[j] == -INFINITY
    // @zu: size num_var array, must satisfy zu[j] >= 0 for all j and
    //      zu[j] == 0 if ub[j] == INFINITY
    // When a starting point was loading successfully (return value 0), then
    // the next call to Solve() will start the IPM from that point, except that
    // primal and dual slacks with value 0 are made positive if necessary. The
    // IPM will skip the initial iterations and start directly with basis
    // preconditioning.
    // At the moment loading a starting point is not possible when the model was
    // dualized during preprocessing. See parameters to turn dualization off.
    // Returns:
    // 0                            success
    // IPX_ERROR_argument_null      an argument was NULL
    // IPX_ERROR_invalid_vector     a sign condition was violated
    // IPX_ERROR_not_implemented    the model was dualized during preprocessing
    Int LoadIPMStartingPoint(const double* x, const double* xl,
                             const double* xu, const double* slack,
                             const double* y, const double* zl,
                             const double* zu);

    // Solves the model that is currently loaded in the object.
    // Returns GetInfo().status.
    Int Solve();

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

    // Returns/sets all parameters. Without calling SetParameters(), the solver
    // uses the default values of a Parameters object.
    Parameters GetParameters() const;
    void SetParameters(Parameters new_parameters);

    // Sets the callback in control
    void SetCallback(HighsCallback* callback);

    // Discards the model and solution (if any) but keeps the parameters.
    void ClearModel();

    // Discards the starting point (if any).
    void ClearIPMStartingPoint();

    // Runs crossover for the given starting point. The starting point must be
    // complementary and satisfy bound and sign conditions; i.e. a dual variable
    // must be non-positive (non-negative) when its primal is not at lower
    // (upper) bound. Each of the pointer arguments can be NULL, in which case
    // all elements of the vector are assumed to be zero.
    // Returns:
    // 0                            starting point OK
    // IPX_ERROR_invalid_vector     starting point not complementary or violates
    //                              bound or sign conditions
    Int CrossoverFromStartingPoint(const double* x_start,
                                   const double* slack_start,
                                   const double* y_start,
                                   const double* z_start);

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
    // - If crossover failed, this is the basis at which failure occurred.
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

  void setTimerOffset(const double offset) {this->control_.setTimerOffset(offset); }

private:
    void ClearSolution();
    void InteriorPointSolve();
    void RunIPM();
    void MakeIPMStartingPointValid();
    void ComputeStartingPoint(IPM& ipm);
    void RunInitialIPM(IPM& ipm);
    void BuildStartingBasis();
    void RunMainIPM(IPM& ipm);
    void BuildCrossoverStartingPoint();
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
    // If crossover_weights_ is non-empty at call to RunCrossover(), then it
    // contains model_.cols() + model_.rows() weights that define the order of
    // primal and dual pushes.
    Vector x_crossover_, y_crossover_, z_crossover_;
    Vector crossover_weights_;
    std::vector<Int> basic_statuses_;

    // IPM starting point provided by user (presolved).
    Vector x_start_, xl_start_, xu_start_, y_start_, zl_start_, zu_start_;
};

}  // namespace ipx

#endif  // IPX_LP_SOLVER_H_

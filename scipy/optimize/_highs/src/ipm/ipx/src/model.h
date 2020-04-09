// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_MODEL_H_
#define IPX_MODEL_H_

#include <vector>
#include "control.h"
#include "sparse_matrix.h"

namespace ipx {

// Model provides the interface between an LP model given by the user,
//
//   minimize   obj'x                                                    (1)
//   subject to A*x {=,<,>} rhs, lbuser <= x <= ubuser,
//
// and the computational form used by the solver,
//
//   minimize   c'x
//   subject to AI*x = b,                              (dual: y)
//              x-xl = lb, xl >= 0,                    (dual: zl >= 0)
//              x+xu = ub, xu >= 0.                    (dual: zu >= 0)
//
// The matrix AI has m >= 0 rows and n+m columns, where n > 0 is the number of
// "structural" variables. The last m columns of AI form the identity matrix.
// The last m components of c do not need to be zero (can happen when the model
// was dualized in preprocessing). Entries of -lb and ub can be infinity.
//
// The user model is translated into computational form in two steps:
// (a) scaling, which consists of
//     - applying an automatic scaling algorithm to A (optional), and
//     - "flipping" variables for which lbuser[j] is infinite but ubuser[j] is
//        finite by multiplying the column of A by -1.
// (b) dualization if appropriate
//
// A Model object cannot be modified other than discarding the data and loading
// a new user model.

class Model {
public:
    // Constructs an empty model.
    Model() = default;

    // Initializes a Model object from the form (1).
    // @num_constr: number of rows of A
    // @num_var: number of columns of A
    // @Ap, @Ai, @Ax: matrix A in CSC format, 0-based indexing
    // @rhs: array of size num_constr
    // @constr_type: array of size num_constr with entries '=', '<' or '>'
    // @obj: array of size num_var
    // @lbuser: array of size num_var, entries can be -INFINITY
    // @ubuser: array of size num_var, entries can be +INFINITY
    // If the input is invalid, info->errflag is set to nonzero and the Model
    // object becomes empty.
    void Load(const Control& control, Int num_constr, Int num_var,
              const Int* Ap, const Int* Ai, const double* Ax,
              const double* rhs, const char* constr_type, const double* obj,
              const double* lbuser, const double* ubuser, Info* info);

    // Returns true if the model is empty.
    bool empty() const { return cols() == 0; }

    // Deallocates all memory; the object becomes empty.
    void clear();

    // Returns the number of rows of AI.
    Int rows() const { return num_rows_; }

    // Returns the number of structural columns of AI (i.e. without the
    // rightmost identity matrix).
    Int cols() const { return num_cols_; }

    // Returns the number of columns classified as dense.
    Int num_dense_cols() const { return num_dense_cols_; }

    // Returns true if column j is classified as dense (0 <= j < n+m).
    bool IsDenseColumn(Int j) const {
        return AI_.entries(j) >= nz_dense_;
    }

    // Returns true if the user model was dualized in preprocessing.
    bool dualized() const { return dualized_; }

    // Returns a reference to the matrix AI in CSC and CSR format.
    const SparseMatrix& AI() const { return AI_; }
    const SparseMatrix& AIt() const { return AIt_; }

    // Returns a reference to a model vector.
    const Vector& b() const { return b_; }
    const Vector& c() const { return c_; }
    const Vector& lb() const { return lb_; }
    const Vector& ub() const { return ub_; }

    // Returns an entry of a model vector.
    double b(Int i) const { return b_[i]; }
    double c(Int j) const { return c_[j]; }
    double lb(Int j) const { return lb_[j]; }
    double ub(Int j) const { return ub_[j]; }

    // Returns the infinity norm of [b; lb; ub], ignoring infinite entries.
    double norm_bounds() const { return norm_bounds_; }

    // Returns the infinity norm of c.
    double norm_c() const { return norm_c_; }

    // Transforms point from user model to solver model. Each of the pointer
    // arguments can be NULL, in which case its components are assumed 0.0.
    void PresolveStartingPoint(const double* x_user, const double* slack_user,
                               const double* y_user, const double* z_user,
                               Vector& x_solver, Vector& y_solver,
                               Vector& z_solver) const;

    // Given an IPM iterate, recovers the solution to the user model (see the
    // reference documentation). Each of the pointer arguments can be NULL, in
    // which case the quantity is not returned. The sign conditions on the dual
    // variables and those on the primal slack variables are staisfied if
    // xl_solver, xu_solver, zl_solver and zu_solver are nonnegative.
    void PostsolveInteriorSolution(const Vector& x_solver,
                                   const Vector& xl_solver,
                                   const Vector& xu_solver,
                                   const Vector& y_solver,
                                   const Vector& zl_solver,
                                   const Vector& zu_solver,
                                   double* x_user,
                                   double* xl_user, double* xu_user,
                                   double* slack_user,
                                   double* y_user,
                                   double* zl_user, double* zu_user) const;

    // Evaluates the solution to the user model obtained from postsolving the
    // IPM iterate. The following info members are set:
    // abs_presidual, abs_dresidual, rel_presidual, rel_dresidual,
    // pobjval, dobjval, rel_objgap, complementarity, normx, normy, normx.
    void EvaluateInteriorSolution(const Vector& x_solver,
                                  const Vector& xl_solver,
                                  const Vector& xu_solver,
                                  const Vector& y_solver,
                                  const Vector& zl_solver,
                                  const Vector& zu_solver,
                                  Info* info) const;

    // Given a basic solution to the solver model, recovers the basic solution
    // to the user model. Each of the pointer arguments can be NULL.
    void PostsolveBasicSolution(const Vector& x_solver,
                                const Vector& y_solver,
                                const Vector& z_solver,
                                const std::vector<Int>& basic_status_solver,
                                double* x_user, double* slack_user,
                                double* y_user, double* z_user) const;

    // Evaluates the solution to the user model obtained from postsolving the
    // basic solution from the solver. The following info members are set:
    // primal_infeas, dual_infeas, objval
    void EvaluateBasicSolution(const Vector& x_solver,
                               const Vector& y_solver,
                               const Vector& z_solver,
                               const std::vector<Int>& basic_status_solver,
                               Info* info) const;

    // Given a basic status for each variable in the solver model, recovers the
    // basic statuses for constraints and variables in the user model. Each
    // of the pointer arguments can be NULL.
    void PostsolveBasis(const std::vector<Int>& basic_status_solver,
                        Int* cbasis, Int* vbasis) const;

private:
    // Checks that the input is valid, and if so copies into the members below
    // (see "User model after scaling"). If the input is invalid, info->errflag
    // is set and the object remains unchanged.
    void CopyInput(Int num_constr, Int num_var, const Int* Ap, const Int* Ai,
                   const double* Ax, const double* rhs, const char* constr_type,
                   const double* obj, const double* lbuser,
                   const double* ubuser, Info* info);

    // Scales A_, scaled_obj_, scaled_rhs_, scaled_lbuser_ and scaled_ubuser_
    // according to parameter control.scale(). The scaling factors are stored in
    // colscale_ and rowscale_. If all factors are 1.0 (either because scaling
    // was turned off or because the algorithm did nothing), rowscale_ and
    // colscale_ have size 0.
    // In any case, variables for which lbuser is infinite but ubbuser is finite
    // are "flipped" and their indices are kept in flipped_vars_.
    void ScaleModel(const Control& control);

    // Builds computational form without dualization. In Julia notation:
    // num_rows = nc
    // num_cols = nv
    // AI       = [A eye(nc)]
    // b        = rhs
    // c        = [obj    ; zeros(nc)                      ]
    // lb       = [lbuser ; constr_type_ .== '>' ? -Inf : 0]
    // ub       = [ubuser ; constr_type_ .== '<' ? +Inf : 0]
    // dualized = false
    // Here nc = num_constr and nv = num_var. The data must have been loaded
    // into the class member below ("User model after scaling") before calling
    // this method.
    void LoadPrimal();

    // Builds computational form with dualization. In Julia notation:
    // num_rows = nv
    // num_cols = nc + nb
    // AI       = [A' -eye(nv)[:,jboxed] eye(nv)]
    // b        = obj
    // c        = [-rhs                          ; ubuser[jb]  ; -lbuser     ]
    // lb       = [constr_type .== '>' ? 0 : -Inf; zeros(nb)   ; zeros(nv)   ]
    // ub       = [constr_type .== '<' ? 0 : +Inf; Inf*ones(nb); Inf*ones(nv)]
    // dualized = true
    // Here nc = num_constr, nv = num_var, nb is the number of boxed variables
    // and jboxed are their indices. Every variable with a finite upper bound
    // must have a finite lower bound (this is ensured after scaling). If a
    // variable j of the input LP is a free variable, then the j-th slack
    // variable of the model gets a zero upper bound (i.e. it is fixed at zero)
    // and its objective coefficient is set to zero.
    void LoadDual();

    // Recursively equilibrates A_ in infinity norm using the algorithm from
    // [1]. The scaling factors are truncated to powers of 2. Terminates when
    // the entries of A_ are within the range [0.5,8).
    // [1] P. A. Knight, D. Ruiz, B. Ucar, "A symmetry preserving algorithm for
    //     matrix scaling", SIAM J. Matrix Anal., 35(3), 2014.
    void EquilibrateMatrix();

    // Initializes num_dense_cols_ and nz_dense_. We classify the maximum #
    // columns as "dense" which have more than 40 nonzeros and more than 10
    // times the # nonzeros than any column that is not "dense". If this yields
    // more than 1000 dense columns, then no columns are classified as dense.
    void FindDenseColumns();

    // Prints the coefficient ranges of input data to control.Log(). Must be
    // called after CopyInput() and before ScaleModel().
    void PrintCoefficientRange(const Control& control) const;

    // Prints preprocessing operations to control.Log().
    void PrintPreprocessingLog(const Control& control) const;

    // Writes statistics of input data and preprocessing to @info.
    void WriteInfo(Info* info) const;

    // ScaleBasicSolution() applies the operations from ScaleModel() to a
    // primal-dual point.
    void ScaleBasicSolution(Vector& x, Vector& slack, Vector& y, Vector& z)
        const;

    // ScaleBack*() do the reverse operation of ScaleModel().
    void ScaleBackInteriorSolution(Vector& x, Vector& xl, Vector& xu,
                                   Vector& slack, Vector& y, Vector& zl,
                                   Vector& zu) const;
    void ScaleBackResiduals(Vector& rb, Vector& rc, Vector& rl,
                            Vector& ru) const;
    void ScaleBackBasicSolution(Vector& x, Vector& slack, Vector& y,
                                Vector& z) const;
    void ScaleBackBasis(std::vector<Int>& cbasis,
                        std::vector<Int>& vbasis) const;

    // DualizeBasicSolution() applies the operations of LoadPrimal() or
    // LoadDual() to a primal-dual point.
    void DualizeBasicSolution(const Vector& x_user, const Vector& slack_user,
                              const Vector& y_user, const Vector& z_user,
                              Vector& x_solver, Vector& y_solver,
                              Vector& z_solver) const;

    // DualizeBack*() do the reverse operations of LoadPrimal() or LoadDual().
    // Given the solution from the solver, they recover the solution to the
    // scaled user model.
    void DualizeBackInteriorSolution(const Vector& x_solver,
                                     const Vector& xl_solver,
                                     const Vector& xu_solver,
                                     const Vector& y_solver,
                                     const Vector& zl_solver,
                                     const Vector& zu_solver,
                                     Vector& x_user,
                                     Vector& xl_user,
                                     Vector& xu_user,
                                     Vector& slack_user,
                                     Vector& y_user,
                                     Vector& zl_user,
                                     Vector& zu_user) const;
    void DualizeBackBasicSolution(const Vector& x_solver,
                                  const Vector& y_solver,
                                  const Vector& z_solver,
                                  Vector& x_user,
                                  Vector& slack_user,
                                  Vector& y_user,
                                  Vector& z_user) const;
    void DualizeBackBasis(const std::vector<Int>& basic_status_solver,
                          std::vector<Int>& cbasis_user,
                          std::vector<Int>& vbasis_user) const;

    void CorrectScaledBasicSolution(Vector& x, Vector& slack, Vector& y,
                                    Vector& z,
                                    const std::vector<Int> cbasis,
                                    const std::vector<Int> vbasis) const;

    // Performs lhs += alpha*A*rhs or lhs += alpha*A'rhs, where A is the user
    // matrix after scaling. This matrix is not stored explicitly, but is used
    // implicitly through AI.
    // @trans: 't' or 'T' for multiplication with A'.
    void MultiplyWithScaledMatrix(const Vector& rhs, double alpha, Vector& lhs,
                                  char trans) const;

    // Computational form model.
    bool dualized_{false};        // model was dualized in preprocessing?
    Int num_rows_{0};             // # rows of AI
    Int num_cols_{0};             // # structural columns of AI
    Int num_dense_cols_{0};       // # columns classified as dense
    Int nz_dense_{0};             // minimum # nonzeros in a dense column
    SparseMatrix AI_;             // matrix AI columnwise
    SparseMatrix AIt_;            // matrix AI rowwise
    Vector b_;
    Vector c_;
    Vector lb_;
    Vector ub_;
    double norm_bounds_{0.0};     // infinity norm of [b;lb;ub]
    double norm_c_{0.0};          // infinity norm of c

    // User model after scaling. The data members are first initialized by
    // CopyInput() and the vectors and matrix are then modified by ScaleModel().
    Int num_constr_{0};           // # constraints
    Int num_eqconstr_{0};         // # equality constraints
    Int num_var_{0};              // # variables
    Int num_free_var_{0};         // # free variables
    Int num_entries_{0};          // # entries in input matrix
    std::vector<Int> boxed_vars_; // indices of boxed variables
    std::vector<char> constr_type_;
    double norm_obj_{0.0};        // Infnorm(obj) as given by user
    double norm_rhs_{0.0};        // Infnorm(rhs,lb,ub) as given by user
    Vector scaled_obj_;
    Vector scaled_rhs_;
    Vector scaled_lbuser_;
    Vector scaled_ubuser_;
    SparseMatrix A_;              // is cleared after preprocessing

    // Data from ScaleModel() that is required by ScaleBack*().
    std::vector<Int> flipped_vars_;
    Vector colscale_;
    Vector rowscale_;
};

// Returns the maximum violation of lb <= x <= ub.
double PrimalInfeasibility(const Model& model, const Vector& x);

// Returns the maximum violation of the dual feasibility condition
//   z[j] <= 0 if x[j] > lb[j],
//   z[j] >= 0 if x[j] < ub[j].
// Note that dual feasibility implies complementarity, i.e.
//   x[j] == lb[j] || x[j] == ub[j] || z[j] == 0.
double DualInfeasibility(const Model& model, const Vector& x, const Vector& z);

// Returns the maximum violation of Ax=b.
double PrimalResidual(const Model& model, const Vector& x);

// Returns the maximum violation of A'y+z=c.
double DualResidual(const Model& model, const Vector& y, const Vector& z);

}  // namespace ipx

#endif  // IPX_MODEL_H_

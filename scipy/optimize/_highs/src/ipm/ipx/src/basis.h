// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_BASIS_H_
#define IPX_BASIS_H_

#include <memory>
#include <vector>
#include "control.h"
#include "indexed_vector.h"
#include "lu_update.h"
#include "model.h"
#include "sparse_matrix.h"

namespace ipx {

// A Basis object is associated with a model. It manages an ordered set of m
// column indices such that AI[:,basis] is nonsingular, where AI is the
// m-by-(n+m) matrix of the model. The class provides the usual simplex-type
// linear algebra operations.

class Basis {
public:
    // Constructor initializes to slack basis. The object stores a reference to
    // a model, which must be valid as long as the basis is used. No data from
    // model is copied.
    Basis(const Control& control, const Model& model);

    // A basis object cannot be copied (at the moment) because the LU
    // factorization is owned by a std::unique_ptr.
    Basis(const Basis&) = delete;
    Basis& operator=(const Basis&) = delete;

    // Move is OK.
    Basis(Basis&&) = default;
    Basis& operator=(Basis&&) = default;

    ~Basis() = default;

    // Returns the index of the variable at position p in the basis.
    Int operator[](Int p) const;

    // If variable j is basic, then returns its position in the basis.
    // Otherwise returns -1.
    Int PositionOf(Int j) const;

    // Returns true if variable j is in the basis.
    bool IsBasic(Int j) const;

    // Returns true if variable j is not in the basis.
    bool IsNonbasic(Int j) const;

    // Each variable has a status from one of the following:
    //
    //  NONBASIC_FIXED: variable is not in the basis (and will never enter)
    //  NONBASIC:       variable is not in the basis (but may enter)
    //  BASIC:          variable is in the basis (but may leave)
    //  BASIC_FREE:     variable is in the basis (and will never leave)
    //
    // The notes in () are not enforced or required by the implementation.
    // From the view point of the class the two nonbasic and the two basic
    // statuses behave the same. The four statuses are provided to make writing
    // user code more convenient.
    enum BasicStatus { NONBASIC_FIXED = -2, NONBASIC, BASIC, BASIC_FREE };

    // Returns status of variable j.
    BasicStatus StatusOf(Int j) const;

    // Switches status of variable j to NONBASIC_FIXED. The variable must not be
    // in the basis.
    void FixNonbasicVariable(Int j);

    // Switches status of variable j to BASIC_FREE. The variable must be in the
    // basis.
    void FreeBasicVariable(Int j);

    // Switches status from NONBASIC_FIXED to NONBASIC for all variables.
    void UnfixVariables();

    // Switches status from BASIC_FREE to BASIC for all variables.
    void UnfreeVariables();

    // Sets the basis to the slack basis.
    void SetToSlackBasis();

    // Loads basis and factorizes it.
    // @basic_status: size n+m array with BasicStatus of each variable
    // Returns: IPX_ERROR_invalid_basis if the basis is invalid (basic_status
    //          contains an invalid entry or # basic variables != m). In this
    //          case the old basis is unchanged.
    //          Otherwise the return code from Factorize() is returned and the
    //          old basis has been replaced. If the given basis is singular, it
    //          will be repaired with slack variables by the LU factorization.
    Int Load(const int* basic_status);

    // Factorizes the current basis matrix from scratch. If nonsingular, a
    // stability check is performed afterwards and the factorization is repeated
    // with a tighter pivot tolerance if the LU factors were unstable.
    // Returns: IPX_ERROR_basis_singular if the factorization was singular and
    //          slack columns have been inserted into the basis. This is not an
    //          "error" from the view point of the Basis object, which remains
    //          in a perfectly valid state.
    //          0 otherwise.
    Int Factorize();

    // Returns true if the LU factorization has not been updated.
    bool FactorizationIsFresh() const;

    // Extracts L, U and permutations such that
    //
    //   B[rowperm,colperm] = (L+I)*U,
    //
    // where B is the current basis matrix, L is strictly lower triangular and
    // U is upper triangular. All arguments can be NULL, in which case the
    // quantity is not returned. Can only be called if FactorizationIsFresh()
    // returns true.
    void GetLuFactors(SparseMatrix* L, SparseMatrix* U, Int* rowperm,
                      Int* colperm) const;

    // Solves linear system with dense right-hand side.
    // @rhs: size m right-hand side vector
    // @lhs: size m solution vector
    // @trans: 't' or 'T' for transposed, other character for forward system
    // @rhs and @lhs may refer to the same object.
    void SolveDense(const Vector& rhs, Vector& lhs, char trans) const;

    // Solves linear system in preparation for update.
    // @j:   column index defining the linear system to be solved.
    //       If j is basic then BTRAN is computed. RHS is the unit vector
    //       corresponding to position of j in the basis.
    //       If j is nonbasic then FTRAN is computed. RHS is AI[:,j].
    // @lhs: returns the solution to the linear system if given.
    void SolveForUpdate(Int j, IndexedVector& lhs);
    void SolveForUpdate(Int j);

    // Computes a row of the (simplex) tableau matrix and performs BTRAN in
    // preparation for an update.
    // @jb:    basic variable. When jb is at position p in the basis, then row p
    //         of the tableau matrix is computed.
    // @btran: BTRAN solution
    // @row:   returns the tableau row. Basic variables have value zero.
    // @ignore_fixed: If true, then the tableau row entry of variables with
    //                status NONBASIC_FIXED is set to zero.
    // The method chooses between a sparse-vector*sparse-matrix and a
    // dense-vector*sparse-matrix operation. Accordingly the pattern of row is
    // or is not set up.
    void TableauRow(Int jb, IndexedVector& btran, IndexedVector& row,
                    bool ignore_fixed = false);

    // Exchanges basic variable jb with nonbasic variable jn if the update to
    // the factorization is stable. In detail, the following steps are done:
    //
    // (1)  Updates the LU factorization.
    // (2a) If the LU update was stable, updates the basis.
    // (2b) If the LU update was unstable, keeps the basis unchanged and
    //      refactorizes, possibly after tightening the LU pivot tolerance.
    //
    // Returns:
    // (2a) The new basis might be refactorized for speed or memory reasons.
    //      In this case returns Factorize(). Otherwise returns 0.
    // (2b) If the factorization was not fresh and/or the pivot tolerance was
    //      tightened, returns Factorize(). Otherwise returns
    //      IPX_ERROR_basis_too_ill_conditioned. In this case the Basis object
    //      is in an inconsistent state. (A basis repair could be done here; at
    //      the moment we simply give up.)
    //
    // At least one of the forward or transposed system must have been solved
    // before in preparation for the update.
    //
    // @jb: Basic variable that leaves the basis. Status of jb becomes NONBASIC.
    // @jn: Nonbasic variable that enters the basis. Status of jn becomes BASIC.
    // @tableau_entry: entry in pivot col and pivot row of the tableau matrix.
    // @sys: > 0 if forward system needs to be solved in preparation for update.
    //       < 0 if transposed sys needs to be solved in preparation for update.
    //         0 if both systems have already been solved.
    // @exchanged: true on return if the basis exchange was performed.
    Int ExchangeIfStable(Int jb, Int jn, double tableau_entry, int sys,
                         bool* exchanged);

    // Computes x[basic], y and z[nonbasic] such that Ax=b and A'y+z=c.
    // @x vector of size n+m. On entry the nonbasic components of x must be set.
    // @y vector of size m, undefined on entry.
    // @z vector of size n+m. On entry the basic components of z must be set.
    // On return the remaining components have been computed.
    void ComputeBasicSolution(Vector& x, Vector& y, Vector& z) const;

    // Constructs a (nonsingular) basis. Given a nonnegative weight for each
    // column of AI, columns with larger weight are preferably chosen as basic
    // columns. Columns with infinite weight will always become basic unless the
    // basis matrix would become singular. Columns with zero weight will always
    // become nonbasic unless (for a slack column) they are required to form a
    // nonsingular basis.
    //
    // If parameter crash_basis is nonzero, a crash procedure is used.
    // Otherwise the method starts with the slack basis and pivots columns with
    // zero or infinite weight out of or into the basis one at a time.
    //
    // @colweights: vector of length n+m with nonnegative entries
    void ConstructBasisFromWeights(const double* colweights, Info* info);

    // Estimates the smallest singular value of the basis matrix.
    double MinSingularValue() const;

    // Computes the # structural nonzero entries per row and column of
    // inverse(B).
    // @rowcounts, @colcounts: either NULL or size m array
    void SymbolicInvert(Int* rowcounts, Int* colcounts) const;

    // Computes the structural density of inverse(B).
    double DensityInverse() const;

    // Returns the associated Model.
    const Model& model() const;

    // Returns statistics.
    Int factorizations() const;       // # LU factorizations
    Int updates_total() const;        // total # basis updates
    double frac_ftran_sparse() const; // fraction of FTRAN solutions sparse
    double frac_btran_sparse() const; // fraction of BTRAN solutions sparse
    double time_factorize() const;    // time LU factorizations
    double time_ftran() const;        // time FTRAN, including partial
    double time_btran() const;        // time BTRAN, including partial
    double time_update() const;       // time LU update
    double mean_fill() const;         // geom. mean of LU fill factors
    double max_fill() const;          // max LU fill factor
    
private:
    // Basis repair terminates when the maximum absolute entry in inverse(B)
    // is smaller than kBasisRepairThreshold. At most kMaxBasisRepair repair
    // operations are performed.
    static constexpr double kBasisRepairThreshold = 1e5;
    static constexpr Int kMaxBasisRepair = 200;

    // Adjusts basis_ and map2basis_ after a singular factorization. Must be
    // called exactly once after the factorization. Returns the # slack
    // variables inserted into the basis (0 if the factorization was
    // nonsingular).
    Int AdaptToSingularFactorization();

    // If possible, tightens the pivot tolerance for subsequent LU
    // factorizations. Returns true if the tolerance was tightened.
    bool TightenLuPivotTol();

    // "Crashes" a basis with preference for variables with larger weight.
    // On return the object has been initialized to a basis that is nonsingular
    // in exact arithmetic. The condition number of the basis matrix can be
    // unacceptably high, however.
    void CrashBasis(const double* colweights);

    // Repairs singularities in the basis matrix by replacing basic columns by
    // slack columns. The status of slack variables that enter the basis becomes
    // BASIC and the status of variables that leave the basis becomes NONBASIC.
    // On return info->basis_repairs >= 0 if repaired successfully, < 0 if
    // failed.
    void Repair(Info* info);

    // Factorizes the basis matrix using a strict absolute pivot tolerance (if
    // supported by the LU implementation). Does not perform the stability check
    // that Factorize() does.
    // @num_dropped: if not NULL, returns the # columns that were dropped from
    //               the basis matrix and replaced by unit columns.
    void CrashFactorize(Int* num_dropped);

    // Similar to ExchangeIfStable() but is guaranteed to exchange jb and jn.
    // If refactorization is required (either for speed or because the LU
    // update was unstable) calls CrashFactorize() on the new basis.
    // @num_dropped: is passed to CrashFactorize() if called; otherwise is set
    //               to 0 if not NULL.
    void CrashExchange(Int jb, Int jn, double tableau_entry, int sys,
                       Int* num_dropped);

    // Pivots free variables into the basis if they can replace a nonfree
    // basic variable. A variable is "free" if its weight is infinite. If a free
    // variable cannot be pivoted into the basis, then its column is linearly
    // dependent on columns of free basic variables. info->dependent_cols
    // reports the # such variables.
    void PivotFreeVariablesIntoBasis(const double* colweights, Info* info);

    // Pivots fixed variables out of the basis if they can be replaced by a
    // nonfixed variable. A variable is "fixed" if its weight is zero. If a
    // fixed variable cannot be pivoted out of the basis, then
    // the rows of AI without that variable are linearly dependent.
    // info->dependent_rows reports the # such variables.
    void PivotFixedVariablesOutOfBasis(const double* colweights, Info* info);

    const Control& control_;
    const Model& model_;
    std::vector<Int> basis_;    // m column indices of AI

    // For 0 <= j < n+m, map2basis_[j] is one of the following:
    //  -2:            variable is NONBASIC_FIXED
    //  -1:            variable is NONBASIC
    //   0 <= p < m:   variable is BASIC and at position p in the basis
    //   m <= p < 2*m: variable is BASIC_FREE and at position p-m in the basis
    std::vector<Int> map2basis_;

    mutable std::unique_ptr<LuUpdate> lu_; // LU factorization of basis matrix
    bool factorization_is_fresh_;  // true if LU factorization not updated

    Int num_factorizations_{0};    // # LU factorizations
    Int num_updates_{0};           // # basis updates
    Int num_ftran_{0};             // # FTRAN operations, excluding partial
    Int num_btran_{0};             // # BTRAN operations, excluding partial
    Int num_ftran_sparse_{0};      // # ... with sparse solution
    Int num_btran_sparse_{0};      // # ... with sparse solution
    double time_ftran_{0.0};       // time for FTRAN ops, including partial
    double time_btran_{0.0};       // time for BTRAN ops, including partial
    double time_update_{0.0};      // time for LU updates
    double time_factorize_{0.0};   // time for LU factorizations
    std::vector<double> fill_factors_; // fill factors from LU factorizations
};

#include <cassert>

inline Int Basis::operator[](Int p) const {
    return basis_[p];
}

inline Basis::BasicStatus Basis::StatusOf(Int j) const {
    const Int m = model_.rows();
    const Int p = map2basis_[j];
    assert(p >= -2 && p < 2*m);
    if (p < 0)
        return p == -1 ? NONBASIC : NONBASIC_FIXED;
    else
        return p < m ? BASIC : BASIC_FREE;
}

inline Int Basis::PositionOf(Int j) const {
    const Int m = model_.rows();
    const Int p = map2basis_[j];
    assert(p >= -2 && p < 2*m);
    return p < 0 ? -1 : p < m ? p : p-m;
}

inline bool Basis::IsBasic(Int j) const {
    return StatusOf(j) == BASIC || StatusOf(j) == BASIC_FREE;
}

inline bool Basis::IsNonbasic(Int j) const {
    return StatusOf(j) == NONBASIC || StatusOf(j) == NONBASIC_FIXED;
}

// Returns x[basis] (in Matlab notation).
Vector CopyBasic(const Vector& x, const Basis& basis);

}  // namespace ipx

#endif  // IPX_BASIS_H_

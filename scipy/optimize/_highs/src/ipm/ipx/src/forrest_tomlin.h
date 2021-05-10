// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_FORREST_TOMLIN_H_
#define IPX_FORREST_TOMLIN_H_

#include <memory>
#include <vector>
#include "control.h"
#include "lu_factorization.h"
#include "lu_update.h"

namespace ipx {

// Generic implementation of the Forrest-Tomlin update [1] that can be used with
// any LU factorization. The implementation does not exploit hypersparsity,
// which exludes its use for such problems. For non-hypersparse problems the
// implementation is better suited than BASICLU, however, because it stores L
// and U in compressed form with permuted indices; hence solving triangular
// systems with a dense rhs/lhs accesses memory contiguously. BASICLU could not
// be implemented in that form due to re-permuting U to triangular form in the
// updates. Hypersparsity support could be added to the present implementation
// when required.
//
// [1] J.J.H. Forrest and J.A. Tomlin, "Updated triangular factors of the basis
//     to maintain sparsity in the product form simplex method", Math.
//     Programming 2 (1972)

class ForrestTomlin : public LuUpdate {
public:
    // Constructor takes ownership of the LuFactorization object by a move from
    // @lu.
    ForrestTomlin(const Control& control, Int dim,
                  std::unique_ptr<LuFactorization>& lu);
    ~ForrestTomlin() = default;

private:
    Int _Factorize(const Int* Bbegin, const Int* Bend, const Int* Bi,
                   const double* Bx, bool strict_abs_pivottol) override;
    void _GetFactors(SparseMatrix* L, SparseMatrix* U, Int* rowperm,
                     Int* colperm, std::vector<Int>* dependent_cols) override;
    void _SolveDense(const Vector& rhs, Vector& lhs, char trans) override;
    void _FtranForUpdate(Int nz, const Int* bi, const double* bx) override;
    void _FtranForUpdate(Int nz, const Int* bi, const double* bx,
                         IndexedVector& lhs) override;
    void _BtranForUpdate(Int j) override;
    void _BtranForUpdate(Int j, IndexedVector& lhs) override;
    Int _Update(double pivot) override;
    bool _NeedFreshFactorization() override;
    double _fill_factor() const override;
    double _pivottol() const override;
    void _pivottol(double new_pivottol) override;

    // Maximum # updates before refactorization is required.
    static constexpr Int kMaxUpdates = 5000;

    // Solves a linear system with the basis matrix. On entry @x holds the
    // permuted right-hand side, on return it holds the permuted solution.
    // @x must have dimension at least dim_ + # computed updates; the additional
    // components are used as workspace.
    // @trans: 't' or 'T' for transposed system.
    void SolvePermuted(Vector& x, char trans);

    // Computes the spike column for the FT update from
    //   R_k^{-1} * ... * R_1^{-1} * L^{-1} * b
    // and stores it in compressed form at the end of U. The spike is also
    // returned as a full vector in work_.
    // @nb, @bi, @bx: b in compressed form.
    void ComputeSpike(Int nb, const Int* bi, const double* bx);

    // Computes the partial BTRAN solution r = ep' * U^{-1}, where ep is the
    // p-th unit vector and p the position of column @j in the pivot sequence.
    // The row eta vector -r/r[p] without the unit diagonal entry is stored in
    // compressed form at end of R. r is returned as a full vector in work_.
    void ComputeEta(Int j);

    const Control& control_;
    const Int dim_;
    std::unique_ptr<LuFactorization> lu_;

    std::vector<Int> rowperm_;     // row permutation from factorization
    std::vector<Int> colperm_;     // col permutation from factorization
    std::vector<Int> rowperm_inv_; // inverse permutation of rowperm_
    std::vector<Int> colperm_inv_; // inverse permutation of colperm_
    std::vector<Int> dependent_cols_;
                             // passed trough from _Factorize() to _GetFactors()

    SparseMatrix L_;            // L from factorization
    SparseMatrix U_;            // U from factorization with spike cols appended
    SparseMatrix R_;            // cols of R build row eta matrices from updates

    // replaced_[k] == p if update k replaced position p in pivot sequence.
    // replaced_.size() is the # updates performed.
    std::vector<Int> replaced_;
    Int replace_next_;          // position to be replaced in next update
    bool have_btran_{false};    // true if row eta has been computed
    bool have_ftran_{false};    // true if spike has been computed
    double fill_factor_{0.0};   // fill factor from last factorization
    double pivottol_{0.1};      // LU pivot tolerance for next factorization
    Vector work_;               // size dim_ + kMaxUpdates workspace
};

}  // namespace ipx

#endif  // IPX_FORREST_TOMLIN_H_

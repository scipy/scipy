// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "splitted_normal_matrix.h"
#include <cassert>
#include <cmath>
#include "timer.h"
#include "utils.h"

namespace ipx {

SplittedNormalMatrix::SplittedNormalMatrix(const Model& model) : model_(model) {
    Int m = model_.rows();
    colperm_.resize(m);
    rowperm_inv_.resize(m);
    work_.resize(m);
}

void SplittedNormalMatrix::Prepare(const Basis& basis, const double* colscale) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const SparseMatrix& AI = model_.AI();
    assert(colscale);
    prepared_ = false;
    N_.clear();                 // deallocate old memory

    basis.GetLuFactors(&L_, &U_, rowperm_inv_.data(), colperm_.data());
    rowperm_inv_ = InversePerm(rowperm_inv_);

    // Scale columns of U.
    for (Int k = 0; k < m; k++) {
        Int p = colperm_[k];
        Int j = basis[p];
        // Nothing to do for BASIC_FREE variables.
        if (basis.StatusOf(j) == Basis::BASIC) {
            double d = colscale[j];
            assert(std::isfinite(d) && d > 0.0);
            ScaleColumn(U_, k, d);
        }
    }

    // Build N with permuted row indices.
    std::vector<Int> nonbasic_vars;
    for (Int j = 0; j < n+m; j++)
        if (basis.StatusOf(j) == Basis::NONBASIC)
            nonbasic_vars.push_back(j);
    N_ = CopyColumns(AI, nonbasic_vars);
    PermuteRows(N_, rowperm_inv_);

    // Scale columns of N.
    for (Int k = 0; k < (Int) nonbasic_vars.size(); k++) {
        Int j = nonbasic_vars[k];
        double d = colscale[j];
        assert(std::isfinite(d));
        ScaleColumn(N_, k, d);
    }

    // Build list of free variables.
    free_positions_.clear();
    for (Int k = 0; k < m; k++) {
        Int p = colperm_[k];
        Int j = basis[p];
        if (basis.StatusOf(j) == Basis::BASIC_FREE)
            free_positions_.push_back(k);
    }
    prepared_ = true;
}

const Int* SplittedNormalMatrix::colperm() const {
    return colperm_.data();
}

double SplittedNormalMatrix::time_B() const {
    return time_B_;
}

double SplittedNormalMatrix::time_Bt() const {
    return time_Bt_;
}

double SplittedNormalMatrix::time_NNt() const {
    return time_NNt_;
}

void SplittedNormalMatrix::reset_time() {
    time_B_ = 0.0;
    time_Bt_ = 0.0;
    time_NNt_ = 0.0;
}

void SplittedNormalMatrix::_Apply(const Vector& rhs, Vector& lhs,
                                  double* rhs_dot_lhs) {
    assert(prepared_);
    Timer timer;

    // Compute work = inverse(B') * rhs.
    work_ = rhs;
    timer.Reset();
    BackwardSolve(L_, U_, work_);
    time_Bt_ += timer.Elapsed();

    // Compute lhs = N*N' * work.
    lhs = 0.0;
    timer.Reset();
    AddNormalProduct(N_, nullptr, work_, lhs);
    time_NNt_ += timer.Elapsed();

    // Compute lhs := inverse(B) * lhs.
    timer.Reset();
    ForwardSolve(L_, U_, lhs);
    time_B_ += timer.Elapsed();

    lhs += rhs;
    for (Int i : free_positions_)
        lhs[i] = 0.0;
    if (rhs_dot_lhs)
        *rhs_dot_lhs = Dot(rhs,lhs);
}

}  // namespace ipx

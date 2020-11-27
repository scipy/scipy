// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "forrest_tomlin.h"
#include <algorithm>
#include <cassert>
#include <utility>
#include "utils.h"

namespace ipx {

ForrestTomlin::ForrestTomlin(const Control& control, Int dim,
                             std::unique_ptr<LuFactorization>& lu) :
    control_(control), dim_(dim) {
    work_.resize(dim_ + kMaxUpdates);
    lu_ = std::move(lu);
}

Int ForrestTomlin::_Factorize(const Int* Bbegin, const Int* Bend, const Int* Bi,
                              const double* Bx, bool strict_abs_pivottol) {
    // Reset updates.
    R_.resize(dim_, 0);
    replaced_.clear();
    replace_next_ = -1;
    have_btran_ = false;
    have_ftran_ = false;

    lu_->Factorize(dim_, Bbegin, Bend, Bi, Bx, pivottol_, strict_abs_pivottol,
                   &L_, &U_, &rowperm_, &colperm_, &dependent_cols_);
    rowperm_inv_ = InversePerm(rowperm_);
    colperm_inv_ = InversePerm(colperm_);

    // Compute fill factor.
    Int bnz = 0;
    for (Int j = 0; j < dim_; j++)
        bnz += Bend[j]-Bbegin[j];
    fill_factor_ = 1.0 * (L_.entries()+U_.entries()) / bnz;

    if (control_.Debug(3)) {
        double normLinv = NormestInverse(L_, "lower", 1);
        double normUinv = NormestInverse(U_, "upper", 0);
        control_.Debug(3)
            << " normLinv = " << sci2(normLinv) << ','
            << " normUinv = " << sci2(normUinv) << ','
            << " stability = " << sci2(lu_->stability()) << '\n';
    }

    Int ret = 0;
    if (lu_->stability() > kLuStabilityThreshold)
        ret |= 1;
    if (!dependent_cols_.empty())
        ret |= 2;
    return ret;
}

void ForrestTomlin::_GetFactors(SparseMatrix* L, SparseMatrix* U, Int* rowperm,
                                Int* colperm, std::vector<Int>* dependent_cols){
    if (L) *L = L_;
    if (U) *U = U_;
    if (rowperm)
        std::copy(rowperm_.begin(), rowperm_.end(), rowperm);
    if (colperm)
        std::copy(colperm_.begin(), colperm_.end(), colperm);
    if (dependent_cols)
        *dependent_cols = dependent_cols_;
}

void ForrestTomlin::_SolveDense(const Vector& rhs, Vector& lhs, char trans) {
    if (trans == 't' || trans == 'T') {
        PermuteBack(colperm_, rhs, work_);
        SolvePermuted(work_, 'T');
        Permute(rowperm_, work_, lhs);
    }
    else {
        PermuteBack(rowperm_, rhs, work_);
        SolvePermuted(work_, 'N');
        Permute(colperm_, work_, lhs);
    }
}

void ForrestTomlin::_FtranForUpdate(Int nz, const Int* bi, const double* bx) {
    ComputeSpike(nz, bi, bx);
}

void ForrestTomlin::_FtranForUpdate(Int nb, const Int* bi, const double* bx,
                                   IndexedVector& lhs) {
    ComputeSpike(nb, bi, bx);
    TriangularSolve(U_, work_, 'n', "upper", 0);

    // Move extra variables from updates to replaced positions.
    Int num_updates = replaced_.size();
    for (Int k = num_updates-1; k >= 0; k--)
        work_[replaced_[k]] = work_[dim_+k];

    // Return lhs without pattern.
    for (Int p = 0; p < dim_; p++)
        lhs[colperm_[p]] = work_[p];
    lhs.InvalidatePattern();
}

void ForrestTomlin::_BtranForUpdate(Int j) {
    ComputeEta(j);
}

void ForrestTomlin::_BtranForUpdate(Int j, IndexedVector& lhs) {
    ComputeEta(j);

    // Apply update etas and solve with L'.
    Int num_updates = replaced_.size();
    for (Int k = num_updates-1; k >= 0; k--) {
        ScatterColumn(R_, k, -work_[dim_+k], work_);
        work_[replaced_[k]] = work_[dim_+k];
        work_[dim_+k] = 0.0;
    }
    TriangularSolve(L_, work_, 't', "lower", 1);

    // Return lhs without pattern.
    for (Int p = 0; p < dim_; p++)
        lhs[rowperm_[p]] = work_[p];
    lhs.InvalidatePattern();
}

// Dot product of entries in queues of A1 and A2. Indices must be sorted in
// increasing order.
static double SparseDot(const SparseMatrix& A1, const SparseMatrix& A2) {
    Int q1 = A1.queue_size();
    Int q2 = A2.queue_size();
    Int p1 = 0;
    Int p2 = 0;
    double d = 0.0;
    while (p1 < q1 && p2 < q2) {
        Int i1 = A1.qindex(p1);
        Int i2 = A2.qindex(p2);
        if (i1 == i2) {
            d += A1.qvalue(p1) * A2.qvalue(p2);
            p1++; p2++;
        } else if (i1 < i2) {
            p1++;
        } else {
            p2++;
        }
    }
    return d;
}

Int ForrestTomlin::_Update(double pivot) {
    Int num_updates = replaced_.size();
    assert(have_ftran_);
    assert(have_btran_);

    // Find the position of the entry with row index replace_next_ in spike.
    // If not present, we will have where == qend.
    assert(replace_next_ >= 0);
    Int qend = U_.queue_size();
    Int where = 0;
    while (where < qend && U_.qindex(where) != replace_next_)
        where++;

    // Compute new diagonal entry of U. newdiag1 will be inserted into U;
    // newdiag2 would be the same in exact arithmetic and is for monitoring
    // numerical stability.
    double olddiag = U_.value(U_.end(replace_next_)-1);
    double newdiag1 = pivot * olddiag;
    double newdiag2 = (where == qend ? 0.0 : U_.qvalue(where)) -
        SparseDot(U_, R_);
    double newdiag_err = std::abs(newdiag1-newdiag2);
    double rel_newdiag_err = newdiag_err / std::abs(newdiag1);

    // Put new diagonal entry at end of spike.
    if (where < qend) {
        for (Int l = where; l < qend-1; l++) {
            U_.qindex(l) = U_.qindex(l+1);
            U_.qvalue(l) = U_.qvalue(l+1);
        }
        U_.qindex(qend-1) = dim_+num_updates;
        U_.qvalue(qend-1) = newdiag1;
    } else {
        U_.push_back(dim_+num_updates, newdiag1);
    }

    // Overwrite replaced column by unit column in U.
    Int end = U_.end(replace_next_);
    for (Int l = U_.begin(replace_next_); l < end-1; l++)
        U_.value(l) = 0.0;
    U_.value(end-1) = 1.0;

    // Finish update.
    U_.add_column();
    R_.add_column();
    replaced_.push_back(replace_next_);
    replace_next_ = -1;
    have_btran_ = false;
    have_ftran_ = false;

    if (newdiag1 == 0.0)
        return -1;

    // Print a debugging message if a new eta entry is large.
    double max_eta = 0.0;
    for (Int l = R_.begin(num_updates); l < R_.end(num_updates); l++)
        max_eta = std::max(max_eta, std::abs(R_.value(l)));
    if (max_eta > 1e10)
        control_.Debug(3) << " max eta = " << sci2(max_eta) << '\n';

    // stability check
    if (rel_newdiag_err > kFtDiagErrorTol) {
        control_.Debug(3)
            << " relative error in new diagonal entry of U = "
            << sci2(rel_newdiag_err) << '\n';
        return 1;
    }
    return 0;
}

bool ForrestTomlin::_NeedFreshFactorization() {
    Int num_updates = replaced_.size();
    Int Rnz = R_.entries();        // nnz in acculumated row etas
    Int Lnz = L_.entries() + dim_; // nnz(L) incl. diagonal
    Int Unz = U_.entries();        // nnz(U) incl. zeroed out columns
    Int U0nz = U_.begin(dim_);     // nnz(U) after factorization

    if (num_updates == kMaxUpdates)
        return true;
    if (num_updates < 100)
        return false;
    if (Rnz > 1.0 * Lnz) {
        // control_.Debug()
        //     << " Rnz = " << Rnz << ", Lnz = " << Lnz
        //     << " after " << num_updates << " updates\n";
        return true;
    }
    if (Unz > 1.7 * U0nz) {
        // control_.Debug()
        //     << " Unz = " << Unz << ", U0nz = " << U0nz
        //     << " after " << num_updates << " updates\n";
        return true;
    }
    return false;
}

double ForrestTomlin::_fill_factor() const {
    return fill_factor_;
}

double ForrestTomlin::_pivottol() const {
    return pivottol_;
}

void ForrestTomlin::_pivottol(double new_pivottol) {
    pivottol_ = new_pivottol;
}

void ForrestTomlin::SolvePermuted(Vector& lhs, char trans) {
    Int num_updates = replaced_.size();
    assert(U_.cols() == dim_+num_updates);

    // Require num_updates elements workspace at end of lhs.
    assert((int)lhs.size() >= dim_+num_updates);

    if (trans == 't' || trans == 'T') {
        // Move replaced entries to the end of the pivot sequence and zero out
        // their old position. Because the corresponding columns of U are unit
        // columns now, the replaced positions remain zero when solving with U'.
        // This is crucial because we neglected to zero out the corresponding
        // rows of U in the update, so there are entries in U which actually
        // should not be there.
        for (Int k = 0; k < num_updates; k++) {
            lhs[dim_+k] = lhs[replaced_[k]];
            lhs[replaced_[k]] = 0.0;
        }
        TriangularSolve(U_, lhs, 't', "upper", 0);
        // Solve backwards with row eta matrices (leading scatter operations)
        // and put the entry from the end of the pivot sequence back into the
        // position that its update replaced.
        for (Int k = num_updates-1; k >= 0; k--) {
            ScatterColumn(R_, k, -lhs[dim_+k], lhs);
            lhs[replaced_[k]] = lhs[dim_+k];
            lhs[dim_+k] = 0.0;
        }
        TriangularSolve(L_, lhs, 't', "lower", 1);
    }
    else {
        TriangularSolve(L_, lhs, 'n', "lower", 1);
        // Solve forward with row eta matrices (leading gather operations) and
        // put the newly computed entry at the end of the pivot sequence. We
        // actually would not need to zero out the replaced positions here
        // because they will not be accessed by any subsequent eta matrix.
        for (Int k = 0; k < num_updates; k++) {
            lhs[dim_+k] = lhs[replaced_[k]] - DotColumn(R_, k, lhs);
            lhs[replaced_[k]] = 0.0;
        }
        // The triangular solve with U fills the replaced position with garbage.
        // This does not hurt because by the fact that the corresponding columns
        // of U have no nonzero off-diagonal entries, the computed values are
        // not propagated further. We finally overwrite them with their actual
        // values from the end of the pivot sequence.
        TriangularSolve(U_, lhs, 'n', "upper", 0);
        for (Int k = num_updates-1; k >= 0; k--) {
            lhs[replaced_[k]] = lhs[dim_+k];
            lhs[dim_+k] = 0.0;
        }
    }
}

void ForrestTomlin::ComputeSpike(Int nb, const Int* bi, const double* bx) {
    Int num_updates = replaced_.size();

    // Solve L*lhs=b.
    work_ = 0.0;
    for (Int p = 0; p < nb; p++)
        work_[rowperm_inv_[bi[p]]] = bx[p];
    TriangularSolve(L_, work_, 'n', "lower", 1);

    // Apply update etas.
    for (Int k = 0; k < num_updates; k++) {
        work_[dim_+k] = work_[replaced_[k]] - DotColumn(R_, k, work_);
        work_[replaced_[k]] = 0.0;
    }

    // Store spike in U. Indices are sorted, which is required for the sparse
    // dot product in Update().
    // Int nz = 0;
    U_.clear_queue();
    for (Int p = 0; p < dim_+num_updates; p++) {
        if (work_[p] != 0.0)
            U_.push_back(p, work_[p]);
    }
    have_ftran_ = true;
}

void ForrestTomlin::ComputeEta(Int j) {
    Int num_updates = replaced_.size();
    assert(U_.cols() == dim_+num_updates);

    // Find permuted position of j.
    Int pos = colperm_inv_[j];
    for (Int k = 0; k < num_updates; k++)
        if (replaced_[k] == pos)
            pos = dim_+k;

    // Solve lhs'U=e_pos'. Replaced position must remain zero (by the fact that
    // the corresponding columns of U are unit), so replaced indices cannot
    // occur in the new eta matrix.
    work_ = 0.0;
    work_[pos] = 1.0;
    TriangularSolve(U_, work_, 't', "upper", 0);
    #ifndef NDEBUG
    for (Int k = 0; k < num_updates; k++)
        assert(work_[replaced_[k]] == 0.0);
    #endif

    // Queue eta at end of R. Indices are sorted, which is required for the
    // sparse dot product in Update().
    // Int nz = 0;
    R_.clear_queue();
    double pivot = work_[pos];
    for (Int i = pos+1; i < dim_+num_updates; i++) {
        if (work_[i] != 0.0)
            R_.push_back(i, -work_[i]/pivot);
    }
    have_btran_ = true;
    replace_next_ = pos;
}

}  // namespace ipx

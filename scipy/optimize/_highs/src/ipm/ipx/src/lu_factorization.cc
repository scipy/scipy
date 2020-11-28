// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "lu_factorization.h"
#include <algorithm>
#include <cassert>
#include "utils.h"

namespace ipx {

// Returns the matrix which in exact arithmetic would be L*U.
// @Bbegin, @Bend, @Bi, @Bx: was the input to LU factorization
// @rowperm, @colperm, @dependent_cols: was the output from LU factorization
static SparseMatrix PermutedMatrix(const Int* Bbegin, const Int* Bend,
                                   const Int* Bi, const double* Bx,
                                   const std::vector<Int>& rowperm,
                                   const std::vector<Int>& colperm,
                                   const std::vector<Int>& dependent_cols) {
    Int dim = rowperm.size();
    std::vector<Int> permuted_row = InversePerm(rowperm);
    std::vector<bool> dependent(dim, false);
    for (Int k : dependent_cols)
        dependent[k] = true;

    SparseMatrix B(dim, 0);
    for (Int k = 0; k < dim; k++) {
        if (dependent[k]) {
            B.push_back(k, 1.0);
        } else {
            Int j = colperm[k];
            for (Int p = Bbegin[j]; p < Bend[j]; p++)
                B.push_back(permuted_row[Bi[p]], Bx[p]);
        }
        B.add_column();
    }
    return B;
}

// Given a strict lower triangular matrix L and an upper triangular matrix U,
// chooses rhs entries +/-1 such that (L+I)\rhs becomes large. Computes
// lhs = U\(L+I)\rhs.
static void SolveForward(const SparseMatrix& L, const SparseMatrix& U,
                         Vector& rhs, Vector& lhs) {
    Int dim = rhs.size();
    lhs = 0.0;
    for (Int i = 0; i < dim; i++) {
        rhs[i] = lhs[i] >= 0.0 ? 1.0 : -1.0;
        lhs[i] += rhs[i];
        ScatterColumn(L, i, -lhs[i], lhs);
    }
    TriangularSolve(U, lhs, 'n', "upper", 0);
}

// Given a strict lower triangular matrix L and an upper triangular matrix U,
// chooses rhs entries +/-1 such that U'\rhs becomes large. Computes
// lhs = (L+I)'\U'\rhs.
static void SolveBackward(const SparseMatrix& L, const SparseMatrix& U,
                          Vector& rhs, Vector& lhs) {
    Int dim = rhs.size();
    lhs = 0.0;
    for (Int j = 0; j < dim; j++) {
        lhs[j] -= DotColumn(U, j, lhs);
        rhs[j] = lhs[j] >= 0.0 ? 1.0 : -1.0;
        lhs[j] += rhs[j];
        Int p = U.end(j)-1;
        assert(U.index(p) == j);
        lhs[j] /= U.value(p);
    }
    TriangularSolve(L, lhs, 't', "lower", 1);
}

// Returns a stability measure for the factorization L*U=B. Considering a linear
// system Bx=b, the solution x computed via LU factorization is "close to" the
// exact solution if the scaled residual
//
//   norm(b-Bx) / (norm(b) + norm(B)*norm(x))
//
// is small (see [1, Section 4.7]). We compute the scaled residual in 1-norm for
// two right-hand sides b1 and b1 and associated x1=B\b1 and x2=B'\b2, and use
// the maximum of the two residuals as stability measure for the factorization.
// This is exactly the method implemented in BASICLU.
//
// [1] I.S. Duff, A.M. Erisman, J.K. Reid, "Direct Methods for Sparse Matrices",
//     second edition (2017).
//
// @Bbegin, @Bend, @Bi, @Bx: was the input to LU factorization
// @L, @U, @rowperm, @colperm, @dependent_cols: output from LU factorization
static double StabilityEstimate(const Int* Bbegin, const Int* Bend,
                                const Int* Bi, const double* Bx,
                                const SparseMatrix& L, const SparseMatrix& U,
                                const std::vector<Int>& rowperm,
                                const std::vector<Int>& colperm,
                                const std::vector<Int>& dependent_cols) {
    Int dim = rowperm.size();
    Vector rhs(dim), lhs(dim);

    // Compute 1-norm and infinity-norm of B.
    SparseMatrix B = PermutedMatrix(Bbegin, Bend, Bi, Bx, rowperm, colperm,
                                    dependent_cols);
    double onenorm = Onenorm(B);
    double infnorm = Infnorm(B);

    SolveForward(L, U, rhs, lhs); // builds some rhs and computes lhs = B\rhs
    double norm_ftran = Onenorm(lhs);
    MultiplyAdd(B, lhs, -1.0, rhs, 'N'); // overwrites rhs by residual
    double norm_ftran_res = Onenorm(rhs);

    SolveBackward(L, U, rhs, lhs); // builds some rhs and computes lhs = B'\rhs
    double norm_btran = Onenorm(lhs);
    MultiplyAdd(B, lhs, -1.0, rhs, 'T'); // overwrites rhs by residual
    double norm_btran_res = Onenorm(rhs);

    return std::max(norm_ftran_res / (dim + onenorm*norm_ftran),
                    norm_btran_res / (dim + infnorm*norm_btran));
}

void LuFactorization::Factorize(Int dim, const Int* Bbegin, const Int* Bend,
                                const Int* Bi, const double* Bx,
                                double pivottol, bool strict_abs_pivottol,
                                SparseMatrix* L, SparseMatrix* U,
                                std::vector<Int>* rowperm,
                                std::vector<Int>* colperm,
                                std::vector<Int>* dependent_cols) {
    _Factorize(dim, Bbegin, Bend, Bi, Bx, pivottol, strict_abs_pivottol,
               L, U, rowperm, colperm, dependent_cols);
    stability_ = StabilityEstimate(Bbegin, Bend, Bi, Bx, *L, *U, *rowperm,
                                   *colperm, *dependent_cols);
}

double LuFactorization::stability() const {
    return stability_;
}

}  // namespace ipx

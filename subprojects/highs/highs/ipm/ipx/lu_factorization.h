#ifndef IPX_LU_FACTORIZATION_H_
#define IPX_LU_FACTORIZATION_H_

#include <vector>
#include "ipm/ipx/ipx_internal.h"
#include "ipm/ipx/sparse_matrix.h"

namespace ipx {

// Interface class for LU factorization. Implementations perform the
// factorization and immediately return the factors. They do not preserve any
// memory in their object.

class LuFactorization {
public:
    virtual ~LuFactorization() {}
    LuFactorization& operator=(const LuFactorization&) = delete;
    LuFactorization& operator=(LuFactorization&&) = delete;

    // Computes the factorization
    //
    //   B[rowperm,colperm] = (L+I)*U
    //
    // into a strict lower triangular matrix L and an upper triangular matrix U.
    // If B is singular, the dependent columns are replaced by unit columns in
    // the product (L+I)*U.
    //
    // @dim: dimension of matrix B
    // @Bbegin, @Bend, @Bi, @Bx: matrix B in sparse column format. The row
    //                           indices and nonzero values of column j are in
    //                           Bi[Bbegin[j]..Bend[j]-1] and
    //                           Bx[Bbegin[j]..Bend[j]-1].
    //                           Indices need not be sorted, but there must be
    //                           no duplicates or invalid entries (need not be
    //                           checked by the implementation).
    // @pivottol: tolerance for partial threshold pivoting, from (0,1].
    // @strict_abs_pivottol: If true, then the implementation is asked to use
    //                       kLuDependencyTol as absolute pivot tolerance and to
    //                       remove columns from the active submatrix
    //                       immediately when all entries became smaller than
    //                       the absolute pivot tolerance. Need not be supported
    //                       by the implementation.
    // @L, @U: return the matrix factors with sorted indices. The objects are
    //         resized as necessary.
    // @rowperm, @colperm: return the row and column permutation. The objects
    //                     are resized for dim entries.
    // @dependent_cols: returns the column indices of B[:,colperm] in which
    //                  no pivot was found.
    //
    // Factorize() cannot fail other than out of memory, in which case
    // std::bad_alloc is thrown.
    void Factorize(Int dim, const Int* Bbegin, const Int* Bend,
                   const Int* Bi, const double* Bx, double pivottol,
                   bool strict_abs_pivottol,
                   SparseMatrix* L, SparseMatrix* U,
                   std::vector<Int>* rowperm, std::vector<Int>* colperm,
                   std::vector<Int>* dependent_cols);

    // Returns a stability measure for the factorization computed in the last
    // call to Factorize(). 0.0 means ultimately stable; a value >> machine
    // precision means that the factorization was numerically unstable, i.e.
    // the relative pivot tolerance was too small.
    double stability() const;

private:
    virtual void _Factorize(Int dim, const Int* Bbegin, const Int* Bend,
                            const Int* Bi, const double* Bx, double pivottol,
                            bool strict_abs_pivottol,
                            SparseMatrix* L, SparseMatrix* U,
                            std::vector<Int>* rowperm,
                            std::vector<Int>* colperm,
                            std::vector<Int>* dependent_cols) = 0;

    double stability_{0.0};
};

}  // namespace ipx

#endif  // IPX_LU_FACTORIZATION_H_
